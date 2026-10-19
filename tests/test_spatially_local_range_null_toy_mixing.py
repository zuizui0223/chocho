"""Toy enumeration audit for the local double-edge swap sampler.

Exactly one introduced occurrence per plant and per region yields 3!=6
possible adjacency matchings. This test confirms that repeated symmetric
local swaps can visit and approximately equilibrate across all six *toy*
worlds when distance and original tethers do not truncate the support.
It does not prove mixing in the real, constrained 30,373-edge problem.
"""
from __future__ import annotations

import importlib.util
import random
import sys
from collections import Counter
from itertools import permutations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SOURCE = ROOT / "scripts/analyze_butterfly_spatially_local_plant_range_null.py"
spec = importlib.util.spec_from_file_location("geographic_null_toy", SOURCE)
assert spec and spec.loader
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_symmetric_swaps_sample_six_toy_permutations() -> None:
    hosts = ("a", "b", "c")
    regions = (0, 1, 2)
    locations = {0:(0.,0.),1:(0.,1.),2:(0.,2.)}
    neighbours = m.neighbour_index(locations,500)
    initial = {h:{u} for h,u in zip(hosts,regions)}
    tether = m.origin_tether(initial,neighbours)
    live = {h:set(v) for h,v in initial.items()}
    groups = [[(h,u) for h,u in zip(hosts,regions)]]
    rng = random.Random(722)
    for _ in range(1000):
        m.local_swap(groups,[3],live,{},neighbours,tether,6,rng)
    count = Counter()
    for _ in range(3000):
        m.local_swap(groups,[3],live,{},neighbours,tether,5,rng)
        assignment = tuple(next(iter(live[h])) for h in hosts)
        count[assignment] += 1
    expected = set(permutations(regions))
    assert set(count) == expected
    assert sum(count.values()) == 3000
    assert all(abs(count[s] / 3000 - 1/6) < 0.05 for s in expected)
    assert sum(len(v) for v in live.values()) == 3
    assert all(u in tether[h] for h,v in live.items() for u in v)
