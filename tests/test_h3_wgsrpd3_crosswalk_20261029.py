"""Synthetic geographic tests, including coastless-centre and dateline cases."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1] / "scripts"))

import pytest
import h3
from shapely.geometry import box
from shapely.ops import transform

from build_h3_wgsrpd3_crosswalk_20261029 import (
    RES,EXPECTED_H3_CELLS,all_global_h3_cells,build_tree,
    fractional_cell,h3_cell_polygon,projector,
)


def simulate(regions,cell):
    tree,geoms,codes=build_tree(regions)
    project,_=projector()
    return fractional_cell(cell,tree,geoms,codes,project)


def test_global_number_and_resolution_are_frozen():
    assert RES==4
    assert EXPECTED_H3_CELLS==h3.get_num_cells(4)
    assert h3.average_hexagon_area(4,unit='km^2')==pytest.approx(1770.347654491,rel=1e-6)
    assert h3.get_resolution(h3.latlng_to_cell(51.5,-0.13,4))==4


def test_center_within_very_small_source_area_mapped_but_not_100percent_of_hex():
    h=h3.latlng_to_cell(51.5,-.13,4)
    lat,lon=h3.cell_to_latlng(h)
    project,_=projector()
    # A tiny original botanical area centred within the H3 cell: only a small
    # fraction of the coarse grid cell is botanically mapped.
    small=project(box(lon-.02,lat-.02,lon+.02,lat+.02))
    record,overlap=simulate({'TST':(small,'1')},h)
    assert record['center_wgsrpd3']=='TST'
    assert record['dominant_area_wgsrpd3']=='TST'
    assert record['mapped_land_fraction_of_cell']<.02
    assert len(overlap)==1 and overlap[0]['wgsrpd3']=='TST'


def test_island_polygon_without_any_hex_center_is_retained():
    h=h3.latlng_to_cell(40,-2,4)
    polygon,centre=h3_cell_polygon(h)
    lon,lat=polygon.exterior.coords[0]
    project,_=projector()
    # A tiny island touches the rim of a hex whose centre remains offshore.
    island=project(box(lon-.03,lat-.03,lon+.03,lat+.03))
    item=simulate({'ISL':(island,'1')},h)
    assert item is not None
    row,overlap=item
    assert row['coastal_cell_center_outside_tdwg']==1
    assert row['center_wgsrpd3']==''
    assert row['mapped_land_fraction_of_cell']<.05


def test_boundary_not_silently_forced_to_one_botanical_region():
    h=h3.latlng_to_cell(51,-1,4)
    poly,centre=h3_cell_polygon(h)
    xmin,ymin,xmax,ymax=poly.bounds
    mid=centre[0]
    project,_=projector()
    left=project(box(xmin-.2,ymin-.2,mid,ymax+.2))
    right=project(box(mid,ymin-.2,xmax+.2,ymax+.2))
    item=simulate({'AAA':(left,'1'),'BBB':(right,'1')},h)
    assert item is not None
    row,pairs=item
    assert row['multiple_region_boundary']==1
    assert len(pairs)==2
    assert sum(v['cell_fraction'] for v in pairs)==pytest.approx(1,abs=.01)
    assert 0<row['dominant_mapped_land_fraction']<1


def test_antimeridian_cell_is_not_giant_false_polygon():
    h=h3.latlng_to_cell(10,179.95,4)
    poly,(lon,lat)=h3_cell_polygon(h)
    assert (poly.bounds[2]-poly.bounds[0])<6
    project,_=projector()
    # The same tiny source lies just across the +180° seam in source geometry.
    target=-179.98 if lon>0 else 179.98
    patch=project(box(target-.03,lat-.03,target+.03,lat+.03))
    out=simulate({'DAT':(patch,'1')},h)
    if out is not None:
        row,pairs=out
        assert row['mapped_land_fraction_of_cell']<.04
        assert all(0<=p['cell_fraction']<=1.015 for p in pairs)


def test_no_actual_source_polygon_is_not_a_planted_host_presence():
    h=h3.latlng_to_cell(55,30,4)
    project,_=projector()
    remote=project(box(-11,49,-10,50))
    assert simulate({'REM':(remote,'1')},h) is None
