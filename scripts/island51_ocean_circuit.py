#!/usr/bin/env python3
"""Source-aligned marine circuits with distributed original-port terminals.

The coastline graph determines allowed connections independently of responses.
AIS is an observation-density scalar and wind is an undirected speed proxy.
Distributed sources preserve separate port nodes; an island is never used as
an artificial zero-resistance bridge through land.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as ss
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

EARTH_RADIUS_KM = 6371.0088


def load_ais(path):
    with np.load(path) as data:
        values = np.asarray(data["AIS_mean"], dtype=float)
        valid = np.asarray(data["source_valid_mask"], dtype=bool)
        transform = np.asarray(data["raster_transform"], dtype=float)
    if transform.shape != (6,) or values.ndim != 2 or valid.shape != values.shape:
        raise ValueError("AIS arrays are not aligned with a source geographic grid")
    step = 360. / values.shape[1]
    expected = np.array([step, 0., -180., 0., -step, 85.])
    if (not np.allclose(transform, expected, atol=1e-12, rtol=0.) or
            not np.isclose(values.shape[0] * step, 170., atol=1e-12)):
        raise ValueError("AIS grid is not aligned with the fixed WGS84 extent")
    if not np.isfinite(values[valid]).all() or np.any(values[valid] < 0):
        raise ValueError("Invalid observed AIS values")
    values[~valid] = np.nan
    return values, valid


def conductance_fields(ais, valid, wind, sea):
    if not (ais.shape == valid.shape == wind.shape == sea.shape):
        raise ValueError("Source marine fields have different grids")
    domain = np.asarray(sea, bool) & np.asarray(valid, bool)
    if not np.isfinite(wind[domain]).all() or np.any(wind[domain] <= 0):
        raise ValueError("Wind must be finite and positive on the common source domain")
    if not np.isfinite(ais[domain]).all() or np.any(ais[domain] < 0):
        raise ValueError("AIS must be observed and nonnegative on the common source domain")
    return domain, {
        "uniform": np.where(domain, 1., np.nan),
        "ais": np.where(domain, 1. + np.log1p(ais), np.nan),
        "wind": np.where(domain, wind, np.nan),
    }


def edge_conductance(longitude, latitude, row, col, cell_conductance):
    """Harmonic endpoint conductance divided by spherical edge distance.

    A caller with finite-volume cells additionally multiplies by the verified
    shared water-boundary width. The returned distance is in kilometers.
    """
    longitude, latitude = np.asarray(longitude), np.asarray(latitude)
    row, col = np.asarray(row), np.asarray(col)
    cell_conductance = np.asarray(cell_conductance)
    if not np.isfinite(cell_conductance).all() or np.any(cell_conductance <= 0):
        raise ValueError("Electrical conductance must be finite and positive")
    lon1, lon2 = np.deg2rad(longitude[row]), np.deg2rad(longitude[col])
    lat1, lat2 = np.deg2rad(latitude[row]), np.deg2rad(latitude[col])
    h = np.sin((lat2-lat1)/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2
    distance = 2*EARTH_RADIUS_KM*np.arcsin(np.sqrt(np.clip(h, 0, 1)))
    if np.any(distance <= 0) or not np.isfinite(distance).all():
        raise ValueError("Electrical edges must have positive geographic length")
    harmonic = 2. / (1. / cell_conductance[row] + 1. / cell_conductance[col])
    return distance, harmonic / distance


def distributed_effective_resistance(laplacian, weights):
    """Return (w_i-w_j)' L^+ (w_i-w_j), retaining source distributions.

    Solve only the common connected component of terminals. Unused isolated
    marine basins are not falsely connected. Every positive-weight terminal
    must already belong to that same source-defined component.
    """
    lap = ss.csr_matrix(laplacian, dtype=float)
    weights = ss.csr_matrix(weights, dtype=float)
    if lap.shape[0] != lap.shape[1] or weights.shape[1] != lap.shape[0]:
        raise ValueError("Laplacian and distributed terminal dimensions differ")
    if weights.shape[0] < 2:
        raise ValueError("At least two source distributions are required")
    if (not np.isfinite(weights.data).all() or np.any(weights.data < 0) or
            not np.allclose(np.asarray(weights.sum(axis=1)).ravel(), 1., atol=1e-12, rtol=0)):
        raise ValueError("Distributed terminal weights must be nonnegative and sum to one")
    skew = lap - lap.T
    if ((skew.nnz and np.max(np.abs(skew.data)) > 1e-10) or
            not np.isfinite(lap.data).all() or
            not np.allclose(np.asarray(lap.sum(axis=1)).ravel(), 0., atol=1e-8)):
        raise ValueError("Invalid symmetric electrical Laplacian")
    adjacency = -lap.copy()
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    if np.any(adjacency.data <= 0):
        raise ValueError("Laplacian contains nonpositive edge conductance")
    _, labels = connected_components(adjacency, directed=False)
    terminals = np.unique(weights.nonzero()[1])
    groups = np.unique(labels[terminals])
    if len(groups) != 1:
        raise ValueError("Source terminal ports are disconnected; no artificial marine bridges allowed")
    component = np.flatnonzero(labels == groups[0])
    lap = lap[component][:, component].tocsc()
    weights = weights[:, component].tocsr()
    # Ground an arbitrary existing node solely to fix the potential gauge.
    keep = np.arange(1, len(component))
    reduced = lap[keep][:, keep].tocsc()
    n = weights.shape[0]
    reference = ss.vstack([weights[0]] * (n-1), format="csr")
    injection = (weights[1:] - reference).T.tocsr()
    rhs_sparse = injection[keep]
    factor = splu(reduced, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.,
                  options={"SymmetricMode": True})
    gram = np.zeros((n, n))
    max_residual = 0.
    for start in range(0, n-1, 4):
        stop = min(start+4, n-1)
        rhs = rhs_sparse[:, start:stop].toarray()
        volts = factor.solve(rhs)
        if not np.isfinite(volts).all():
            raise ValueError("Circuit solver produced nonfinite potentials")
        residual = reduced @ volts - rhs
        denom = np.maximum(np.max(np.abs(rhs), axis=0), 1e-30)
        max_residual = max(max_residual, float(np.max(np.max(np.abs(residual), axis=0)/denom)))
        gram[1:, start+1:stop+1] = rhs_sparse.T @ volts
    if not np.isfinite(gram).all() or not np.isfinite(max_residual):
        raise ValueError("Circuit solver produced nonfinite energy or residuals")
    asymmetry = float(np.max(np.abs(gram-gram.T)))
    gram = (gram+gram.T)/2
    diagonal = np.diag(gram)
    resistance = diagonal[:, None] + diagonal[None, :] - 2*gram
    scale = max(float(np.max(np.abs(resistance))), 1.)
    if resistance.min() < -1e-9*scale or max_residual > 1e-7:
        raise ValueError("Circuit numerical accuracy gate failed")
    resistance = np.maximum(resistance, 0.)
    np.fill_diagonal(resistance, 0.)
    return resistance, {
        "solver": "scipy sparse grounded LU, batched distributed unit-current right-hand sides",
        "island_terminal_distributions": n,
        "unique_terminal_nodes": len(terminals),
        "solved_component_nodes": len(component),
        "max_linear_system_relative_residual": max_residual,
        "gram_max_asymmetry_before_symmetrization": asymmetry,
        "ports_are_not_short_circuited": True,
    }


def effective_resistance(laplacian, terminals):
    terminals = np.asarray(terminals, dtype=int)
    weights = ss.csr_matrix((np.ones(len(terminals)), (np.arange(len(terminals)), terminals)),
                            shape=(len(terminals), laplacian.shape[0]))
    return distributed_effective_resistance(laplacian, weights)


def resistance_feature(resistance, row, col):
    """Geography-only median normalization removes arbitrary electrical units."""
    upper = resistance[np.triu_indices(len(resistance), 1)]
    positive = upper[upper > 0]
    if not len(positive):
        raise ValueError("All source terminal distributions are electrically identical")
    return np.log1p(resistance[row, col]/np.median(positive))


def port_weights(island_ids, groups, port_to_node, n_nodes):
    rows, cols, values = [], [], []
    seen = set()
    for i, island in enumerate(island_ids):
        ports = groups[island]
        if not ports or len(set(ports)) != len(ports) or seen.intersection(ports):
            raise ValueError("Original single-landmass port groups are empty or duplicated")
        seen.update(ports)
        for port in ports:
            if port not in port_to_node:
                raise ValueError("An original source port is missing from marine attachments")
            node = int(port_to_node[port])
            if node < 0 or node >= n_nodes:
                raise ValueError("Port attachment lies outside source marine graph")
            rows.append(i)
            cols.append(node)
            values.append(1./len(ports))
    # CSR construction sums equal node entries, retaining equal original-port
    # weights instead of assigning equal weights to unique coarse cells.
    return ss.csr_matrix((values, (rows, cols)), shape=(len(island_ids), n_nodes))


def nearest_sea_cell(longitude, latitude, cell_to_node, step, *, west=-180., north=85.,
                     distance_limit_km=25.):
    """Independently attach an actual port to the nearest resolved sea center.

    A bounded geographic window avoids a global point tree. Longitude wraps
    only for a true global grid. No alternative node is chosen for uniqueness
    or to select a desired marine component.
    """
    if not np.isfinite([longitude, latitude, step, distance_limit_km]).all():
        raise ValueError("Nonfinite original marine attachment parameters")
    if step <= 0 or distance_limit_km <= 0 or abs(latitude) > 90:
        raise ValueError("Invalid original marine attachment parameters")
    height, width = cell_to_node.shape
    angular = np.rad2deg(distance_limit_km/EARTH_RADIUS_KM)
    # Bound using the smallest cosine anywhere within the candidate latitude
    # interval; this cannot exclude the nearest feasible cell near the poles.
    cosine = max(np.cos(np.deg2rad(min(abs(latitude)+angular, 89.999))), 1e-5)
    radius_y = int(np.ceil(angular/step)) + 1
    radius_x = int(np.ceil(angular/(step*cosine))) + 1
    y0 = int(np.floor((north-latitude)/step))
    x0 = int(np.floor((longitude-west)/step))
    ys = np.arange(max(0, y0-radius_y), min(height, y0+radius_y+1))
    raw_x = np.arange(x0-radius_x, x0+radius_x+1)
    if np.isclose(width*step, 360.):
        xs = np.unique(raw_x % width)
    else:
        xs = raw_x[(raw_x >= 0) & (raw_x < width)]
    if not len(ys) or not len(xs):
        raise ValueError("No verified marine center within source-port distance bound")
    yy, xx = np.meshgrid(ys, xs, indexing="ij")
    allowed = cell_to_node[yy, xx] >= 0
    yy, xx = yy[allowed], xx[allowed]
    if not len(yy):
        raise ValueError("No verified marine center within source-port distance bound")
    lon = west+(xx+.5)*step
    lat = north-(yy+.5)*step
    dlat = np.deg2rad(lat-latitude)
    dlon = np.deg2rad(lon-longitude)
    h = np.sin(dlat/2)**2 + np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(latitude))*np.sin(dlon/2)**2
    distances = 2*EARTH_RADIUS_KM*np.arcsin(np.sqrt(np.clip(h, 0, 1)))
    k = int(np.argmin(distances))
    if distances[k] > distance_limit_km:
        raise ValueError("Nearest verified marine center exceeds source-port distance bound")
    return {"node": int(cell_to_node[yy[k], xx[k]]),
            "raster_row": int(yy[k]), "raster_col": int(xx[k]),
            "marine_longitude": float(lon[k]), "marine_latitude": float(lat[k]),
            "displacement_km": float(distances[k])}
