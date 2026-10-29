#!/usr/bin/env python3
"""Build a source-verified GLOBAL H3-r4 ↔ WGSRPD3 polygon-area crosswalk.

Each H3 index is geography only, never proof an accepted host plant is present.
Enumerating *all* H3 r4 cells captures islands/coasts with no H3 centre on land.
Original botanical source polygons are not redistributed; only derived cell IDs,
WGSRPD codes and land overlap area fractions are written.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv,gzip,hashlib,io,json,math
from pathlib import Path

import h3
from pyproj import Transformer
from shapely import make_valid
from shapely.affinity import translate
from shapely.geometry import Point,Polygon,shape
from shapely.ops import transform,unary_union
from shapely.strtree import STRtree

from build_country_area_atoms_20261020 import atom_areas

SCHEMA='chocho_global_h3r4_WGSRPD3_fractional_geometry_bridge_20261029_v01'
RES=4
EXPECTED_REGIONS=369
# Full planet: exactly 2 + 120*7**4, including the 12 pentagons.
EXPECTED_H3_CELLS=2+120*7**4
MIN_INTERSECTION_M2=1.0
MAX_PROJECTED_COVERAGE=1.015  # material errors fail instead of clamping source geometry


def digest(path:Path)->str:
    d=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):d.update(block)
    return d.hexdigest()


def polygons_only(g):
    q=make_valid(g)
    if q.geom_type in ('Polygon','MultiPolygon'):return q
    if q.geom_type=='GeometryCollection':
        parts=[p for x in q.geoms for p in (list(x.geoms) if x.geom_type=='MultiPolygon' else [x])
               if p.geom_type=='Polygon' and p.area>0]
        return unary_union(parts) if parts else Polygon()
    return Polygon()


def projector():
    # Existing chocho country→WGSRPD audit uses this same equal-area projection
    # and force_over=True for original MAG polygon -180.000026° boundary.
    t=Transformer.from_crs('EPSG:4326','EPSG:6933',always_xy=True,force_over=True)
    return lambda g: transform(t.transform,g),t


def load_source(input_geojson,source_summary):
    manifest=json.loads(Path(source_summary).read_text(encoding='utf8'))
    if manifest.get('schema')!='country_wgsrpd_area_crosswalk_v1':
        raise ValueError('Wrong previously source-audited WGSRPD3 geography manifest')
    src_sha=digest(input_geojson)
    if src_sha!=manifest['source']['tdwg_sha256']:
        raise ValueError('Author-only Level3 source geometry does not match audited prior crosswalk SHA256')
    source=json.loads(Path(input_geojson).read_text(encoding='utf8'))
    groups=defaultdict(list)
    levels={}
    duplicate=0
    for f in source['features']:
        p=f['properties'];code=str(p['LEVEL3_COD']);level=str(p['LEVEL1_COD'])
        if code in levels and levels[code]!=level:
            raise ValueError('WGSRPD3 region has incompatible hierarchy: '+code)
        if code in groups:duplicate+=1
        levels[code]=level
        g=polygons_only(shape(f['geometry']))
        if g.is_empty:raise ValueError('Empty original WGSRPD3 region '+code)
        groups[code].append(g)
    if len(groups)!=EXPECTED_REGIONS:
        raise ValueError('Expected fixed 369-region original WGSRPD3, got '+str(len(groups)))
    project,t=projector()
    out={}
    for code,parts in sorted(groups.items()):
        g=polygons_only(unary_union(parts))
        gp=polygons_only(project(g))
        if gp.is_empty or gp.area<100:
            raise ValueError('Original region projection lost '+code)
        out[code]=(gp,levels[code])
    return out,{'tdwg_geojson_sha256':src_sha,'tdwg_source_features':len(source['features']),
                'tdwg_regions':len(groups),'source_duplicate_code_features':duplicate,
                'prior_original_country_crosswalk_sha256':digest(source_summary)}


def build_tree(regions):
    """Include shifted copies to cover dateline cells without flattening a 360° seam."""
    _,trans=projector()
    x0=trans.transform(0.,0.)[0];x360=trans.transform(360.,0.)[0]
    period=x360-x0
    if not(3.0e7<period<4.5e7):raise ValueError('Projected longitude full-turn scale changed')
    geoms=[];codes=[]
    for code,(g,level) in sorted(regions.items()):
        for shift in (-period,0.,period):
            geoms.append(translate(g,xoff=shift));codes.append(code)
    return STRtree(geoms),geoms,codes


def h3_cell_polygon(cell):
    """Unwrap longitude about H3 centre before force_over equal-area projection."""
    clat,clon=h3.cell_to_latlng(cell)
    vertices=[]
    for lat,lon in h3.cell_to_boundary(cell):
        lon2=lon+360*round((clon-lon)/360)
        vertices.append((float(lon2),float(lat)))
    g=Polygon(vertices)
    if not g.is_valid:g=polygons_only(g)
    if g.is_empty or g.area<=0:raise ValueError('Invalid H3 geographical boundary '+cell)
    return g,(clon,clat)


def fractional_cell(cell,tree,geoms,codes,project):
    polygon,centre=h3_cell_polygon(cell)
    projected=polygons_only(project(polygon))
    cell_area=projected.area
    if not math.isfinite(cell_area) or cell_area<=100:
        raise ValueError('Invalid equal-area projected H3 cell '+cell)
    ix=tree.query(projected,predicate='intersects')
    if len(ix)==0:return None
    pairs=defaultdict(list)
    for idx in ix:
        j=int(idx)
        overlap=projected.intersection(geoms[j])
        if overlap.area>=MIN_INTERSECTION_M2:
            pairs[codes[j]].append(overlap)
    if not pairs:return None
    # Several geographic aliases may represent one code at ±180°. Union their
    # overlaps per original code, instead of double-counting shifted polygons.
    region_overlaps={code:polygons_only(unary_union(gg)) for code,gg in pairs.items()}
    positive={code:g for code,g in region_overlaps.items() if g.area>=MIN_INTERSECTION_M2}
    if not positive:return None
    areas={code:float(g.area) for code,g in positive.items()}
    mapped_area=unary_union(list(positive.values())).area
    if mapped_area>MAX_PROJECTED_COVERAGE*cell_area:
        raise ValueError('WGSRPD3 mapped land > 101.5% of cell area; inspect dateline and source topology '+cell)
    cp=Point(*centre)
    cpp=project(cp)
    center_codes=sorted({codes[int(idx)] for idx in tree.query(cpp,predicate='intersects')
                        if geoms[int(idx)].covers(cpp)})
    center=center_codes[0] if len(center_codes)==1 else ''
    rank=sorted(areas,key=lambda code:(-areas[code],code))
    dominant=rank[0]
    dominant_fraction=areas[dominant]/mapped_area if mapped_area else 0
    # An unambiguous coarse region is assigned only when its area dominates
    # mapped land. Raw centre result remains separate from that area assignment.
    majority=dominant if dominant_fraction>.50 and len(rank)==1 or (
        dominant_fraction>.50 and areas[rank[0]]>areas[rank[1]]+1e-5) else ''
    # Geographic memberships can overlap; disjoint membership atoms preserve
    # source union areas and make downstream native vs introduced resources
    # additive WITHOUT silently double-counting border areas.
    if len(positive)==1:
        only=next(iter(positive))
        atoms={(only,):mapped_area}
    else:
        atoms=atom_areas(projected,positive)
        if not math.isclose(sum(atoms.values()),mapped_area,rel_tol=1e-7,abs_tol=10.):
            raise ValueError('Disjoint WGSRPD atoms do not reconstruct the cell union area '+cell)
    atom_rows=[{'h3_r4':cell,'wgsrpd3_members':'|'.join(members),
                'cell_fraction':area/cell_area,
                'mapped_land_fraction':area/mapped_area,
                'atomic_intersection_km2':area/1e6}
               for members,area in sorted(atoms.items()) if area>=MIN_INTERSECTION_M2]
    rows=[{'h3_r4':cell,'wgsrpd3':code,'cell_fraction':areas[code]/cell_area,
           'mapped_land_fraction':areas[code]/mapped_area,'intersection_km2':areas[code]/1e6}
          for code in sorted(areas)]
    return {
       'h3_r4':cell,'center_wgsrpd3':center,
       'center_match_count':len(center_codes),
       'dominant_area_wgsrpd3':majority,
       'dominant_mapped_land_fraction':dominant_fraction,
       'mapped_land_fraction_of_cell':mapped_area/cell_area,
       'cell_area_km2':cell_area/1e6,'mapped_area_km2':mapped_area/1e6,
       'n_intersecting_wgsrpd3':len(rows),
       'multiple_region_boundary':int(len(rows)>1),
       'coastal_cell_center_outside_tdwg':int(len(center_codes)==0),
       'source_region_overlap_excess_fraction':sum(areas.values())/cell_area-mapped_area/cell_area,
     },rows,atom_rows


def all_global_h3_cells(res=RES):
    """Never mask via centres before intersection: retain subcell island areas."""
    for base in sorted(h3.get_res0_cells()):
        for cell in sorted(h3.cell_to_children(base,res)):
            yield cell


def deterministic_csv_gz(path,columns,rows):
    with Path(path).open('wb') as target:
        with gzip.GzipFile(filename='',fileobj=target,mode='wb',compresslevel=7,mtime=0) as z:
            with io.TextIOWrapper(z,encoding='utf-8',newline='') as w:
                writer=csv.DictWriter(w,fieldnames=columns)
                writer.writeheader()
                for row in rows:writer.writerow(row)


def run(args):
    regions,source=load_source(args.tdwg,args.previous_crosswalk)
    tree,geoms,codes=build_tree(regions)
    project,_=projector()
    outdir=Path(args.output_dir);outdir.mkdir(parents=True,exist_ok=True)
    # Stream compressed tables instead of materializing all ~100k cells. Avoid
    # non-reproducible gzip timestamps to permit byte-identical replay.
    cells_path=outdir/'global_h3r4_wgsrpd3_cells_20261029.csv.gz'
    rows_path=outdir/'global_h3r4_wgsrpd3_overlap_20261029.csv.gz'
    atoms_path=outdir/'global_h3r4_wgsrpd3_atoms_20261029.csv.gz'
    cc=['h3_r4','center_wgsrpd3','center_match_count','dominant_area_wgsrpd3',
        'dominant_mapped_land_fraction','mapped_land_fraction_of_cell',
        'cell_area_km2','mapped_area_km2','n_intersecting_wgsrpd3',
        'multiple_region_boundary','coastal_cell_center_outside_tdwg',
        'source_region_overlap_excess_fraction']
    rc=['h3_r4','wgsrpd3','cell_fraction','mapped_land_fraction','intersection_km2']
    ac=['h3_r4','wgsrpd3_members','cell_fraction','mapped_land_fraction','atomic_intersection_km2']
    # Two output streams are created simultaneously to avoid recomputing global
    # intersections. Gzip header contains no timestamp or original filename.
    import contextlib
    counters=Counter(); regions_seen=set();max_excess=0.;area_sum=0.
    with contextlib.ExitStack() as stack:
        a=stack.enter_context(cells_path.open('wb'))
        b=stack.enter_context(rows_path.open('wb'))
        c=stack.enter_context(atoms_path.open('wb'))
        ag=stack.enter_context(gzip.GzipFile(filename='',fileobj=a,mode='wb',compresslevel=7,mtime=0))
        bg=stack.enter_context(gzip.GzipFile(filename='',fileobj=b,mode='wb',compresslevel=7,mtime=0))
        cg=stack.enter_context(gzip.GzipFile(filename='',fileobj=c,mode='wb',compresslevel=7,mtime=0))
        tw=stack.enter_context(io.TextIOWrapper(ag,encoding='utf-8',newline=''))
        rw=stack.enter_context(io.TextIOWrapper(bg,encoding='utf-8',newline=''))
        aw=stack.enter_context(io.TextIOWrapper(cg,encoding='utf-8',newline=''))
        cw=csv.DictWriter(tw,fieldnames=cc);cw.writeheader()
        ow=csv.DictWriter(rw,fieldnames=rc);ow.writeheader()
        atom_writer=csv.DictWriter(aw,fieldnames=ac);atom_writer.writeheader()
        for cell in all_global_h3_cells():
            counters['global_H3_r4_cells_tested']+=1
            mapped=fractional_cell(cell,tree,geoms,codes,project)
            if not mapped:continue
            row,overlap,atoms=mapped
            cw.writerow(row);ow.writerows(overlap);atom_writer.writerows(atoms)
            counters['cells_any_tdwg_land_intersection']+=1
            counters['cell_region_positive_overlap_rows']+=len(overlap)
            counters['source_disjoint_membership_atom_rows']+=len(atoms)
            counters['source_cells_with_actual_region_area_overlap']+=int(row['source_region_overlap_excess_fraction']>1e-6)
            counters['cells_with_multiple_code_atomic_membership']+=int(any('|' in a['wgsrpd3_members'] for a in atoms))
            counters['boundary_mixed_region_cells']+=row['multiple_region_boundary']
            counters['coastal_cells_without_tdwg_centroid']+=row['coastal_cell_center_outside_tdwg']
            counters['cells_center_multiple_wgsrpd3']+=int(row['center_match_count']>1)
            counters['cells_area_majority_unresolved']+=int(not row['dominant_area_wgsrpd3'])
            if row['mapped_land_fraction_of_cell']>.5:
                counters['cells_mapped_over_half_area']+=1
            max_excess=max(max_excess,row['source_region_overlap_excess_fraction'])
            area_sum+=row['mapped_area_km2']
            regions_seen.update(z['wgsrpd3'] for z in overlap)
            if counters['global_H3_r4_cells_tested']%50000==0:
                print(json.dumps({'h3_tested':counters['global_H3_r4_cells_tested'],
                    'mapped_land_cells':counters['cells_any_tdwg_land_intersection']}),flush=True)
    if counters['global_H3_r4_cells_tested']!=h3.get_num_cells(RES) or counters['global_H3_r4_cells_tested']!=EXPECTED_H3_CELLS:
        raise ValueError('Global H3 r4 cell count differs from published H3 index hierarchy')
    if len(regions_seen)!=len(regions):
        raise ValueError('Some original WGSRPD Level3 region has zero H3 intersection; source resolution/geometry incomplete')
    if counters['cells_any_tdwg_land_intersection']<25000:
        raise ValueError('Implausibly low worldwide H3↔botanical land intersections')
    if area_sum<1.0e8 or area_sum>2.0e8:
        raise ValueError('Mapped botanical land area outside global terrestrial bound')
    report={
        'schema':SCHEMA,'status':'SOURCE_VERIFIED_GLOBAL_H3R4_WGSRPD3_AREA_CROSSWALK_COMPLETE',
        'h3_resolution':RES,'all_H3_planet_cells_including_sea_tested':counters['global_H3_r4_cells_tested'],
        'h3_average_hex_area_km2':h3.average_hexagon_area(RES,unit='km^2'),
        'h3_average_hex_edge_km':h3.average_hexagon_edge_length(RES,unit='km'),
        'h3_library_version':h3.versions(),
        'equal_area_crs':'EPSG:6933','geographic_crs':'EPSG:4326',
        'longitude_projection_policy':'force_over=True, H3 hex vertices unwrapped about centre, original WGSRPD polygons shifted ±360° for dateline intersection',
        'original_source':source,
        'coverage':dict(sorted(counters.items())),
        'regions_with_at_least_one_area_intersection':len(regions_seen),
        'sum_mapped_botanical_land_area_km2':area_sum,
        'maximum_source_region_overlap_excess_fraction':max_excess,
        'cell_csv_gzip_sha256':digest(cells_path),
        'cell_region_overlaps_gzip_sha256':digest(rows_path),
        'disjoint_atoms_gzip_sha256':digest(atoms_path),
        'disjoint_atoms_csv_filename':atoms_path.name,
        'disjoint_atoms_union_verification':'The sum of atom areas equals the union of source-region land within each H3; support can be aggregated by selecting any atom containing a specified WCVP region once.',
        'cell_csv_filename':cells_path.name,'cell_region_csv_filename':rows_path.name,
        'source_raw_geometry_not_exported':True,
        'per_h3_land_fraction_is_not_actual_host_plant_presence':True,
        'wgsrpd3_membership_fixed_for_all_years':True,
        'source_region_overlap_not_interpreted_as_independent_resource_area':True,
        'GBIF_detection_nondetection_not_computed':True,
        'r4_coastal_boundary_partials_retained_even_without_centroid':True,
        'no_species_host_status_or_remote_sensing_values_joined_yet':True,
        'future_grid_year_key':['h3_r4','year'],
        'inference':'The table describes geographic membership under WGSRPD3 original polygon generalized limits; does NOT establish local host occurrence, managed status, host introduction year, butterfly occupancy, or future climate response.'
    }
    (outdir/'global_h3r4_wgsrpd3_crosswalk_manifest_20261029.json').write_text(
        json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':report['status'],'coverage':report['coverage'],
                      'source_regions':len(regions),
                      'sum_mapped_km2':round(area_sum,1),
                      'cell_csv_sha256':report['cell_csv_gzip_sha256']}),flush=True)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--tdwg',required=True,type=Path)
    p.add_argument('--previous-crosswalk',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path)
    args=p.parse_args()
    run(args)
