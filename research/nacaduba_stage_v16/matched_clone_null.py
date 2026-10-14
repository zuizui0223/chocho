"""Reproduce the exact 4-clone/6-ramet provenance-label rarity diagnostic.

Source: Kato et al (2012), Silvae Genetica, Table 7, DOI 10.1515/sg-2012-0029.
A group of two single-ramet clones and two double-ramet clones is selected
from the 23 observed clones. This is NOT randomized causal inference, host
flower-bud biomass, or a butterfly survival study.
"""
from pathlib import Path
import csv
import itertools
import json
import hashlib
from collections import Counter

P = Path(__file__).resolve().parents[1] / "nacaduba_stage_v15" / "data"
files = {"clones": P / "clones_kato2012_table4_table5.csv",
         "ramets": P / "ramets_kato2012_table7.csv"}
blobs = {"clones": "47f30f98807f27c3f47f8797731f7e05cc9eca9f",
         "ramets": "ed8425db6cbc6830249824a2edd41fbfc7ec2db3"}
for key, path in files.items():
    content = path.read_bytes()
    actual = hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()
    assert actual == blobs[key], (key, actual)

with files["clones"].open(newline="") as f:
    clones = list(csv.DictReader(f))
with files["ramets"].open(newline="") as f:
    ramets = list(csv.DictReader(f))
assert len(clones) == 23 and len(ramets) == 36
assert sum(int(c["ramets"]) for c in clones) == 39
years = (2008, 2009, 2010)
months = ["Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar"]
parts = ["early","mid","late"]
ramet_by_clone = {c["clone"]: [] for c in clones}
masks = {c["clone"]: {year: 0 for year in years} for c in clones}
for row in ramets:
    c = row["clone"]
    ramet_by_clone[c].append(row)
    for year in years:
        b = row["begin_" + str(year)]
        if not b:
            assert not row["duration_" + str(year)]
            continue
        part, month = b.split()
        start = months.index(month) * 3 + parts.index(part)
        length = int(row["duration_" + str(year)])
        assert start + length <= 36
        masks[c][year] |= ((1 << length) - 1) << start
all_ids = set(masks)
real_group = frozenset(c["clone"] for c in clones if c["seed_source"] == "g")
assert real_group == {"146", "147", "148", "149"}

def exclusive_and_lead(group, removed=None):
    active = all_ids - ({removed} if removed else set())
    s1 = active & set(group)
    s2 = active - s1
    score = []
    leads = []
    for year in years:
        a = 0
        b = 0
        for c in s1: a |= masks[c][year]
        for c in s2: b |= masks[c][year]
        score.append((a & ~b).bit_count())
        def first_bit(v): return (v & -v).bit_length() - 1 if v else None
        ea, eb = first_bit(a), first_bit(b)
        leads.append(max(eb-ea, 0) if ea is not None and eb is not None else 0)
    return score, leads

def all_six_flowered_all_years(group):
    for c in group:
        records = ramet_by_clone[c]
        required = int(next(x["ramets"] for x in clones if x["clone"] == c))
        if len(records) != required:
            return False
        if any(not row["begin_" + str(y)] for row in records for y in years):
            return False
    return True

one = [c["clone"] for c in clones if int(c["ramets"]) == 1]
two = [c["clone"] for c in clones if int(c["ramets"]) == 2]
assert (len(one),len(two)) == (7,16)
rows = []
for singles in itertools.combinations(one,2):
    for doubles in itertools.combinations(two,2):
        selection = frozenset(singles + doubles)
        scores, leads = exclusive_and_lead(selection)
        rows.append((selection, scores, leads, all_six_flowered_all_years(selection)))
assert len(rows) == 2520
obs, lead = exclusive_and_lead(real_group)
assert obs == [9,6,8] and lead == [9,6,8]
scores = [sum(x[1]) for x in rows]
assert max(scores) == sum(obs) == 23 and sum(s>=23 for s in scores) == 1
conditional = [r for r in rows if r[3]]
assert len(conditional) == 15 and sum(sum(r[1])>=23 for r in conditional) == 1
drop_one = {}
for c in sorted(real_group):
    exclusive, lead = exclusive_and_lead(real_group-{c},removed=c)
    drop_one[c] = {"exclusive":exclusive,"lead":lead}
assert {k:sum(v["exclusive"]) for k,v in drop_one.items()} == {
    "146":22,"147":23,"148":21,"149":23
}
summary = {
    "data":"Pinned Kato2012 Acacia auriculiformis clone Table 4/5/7",
    "plant_clones":23,"plant_ramets":39,"seasons":3,
    "real_source_g_clones":sorted(real_group),
    "unique_flowering_tenday_bins_by_season":obs,
    "count_matched_groups":len(rows),
    "group_membership_constraint":"two clones with one ramet each + two clones with two ramets each",
    "matched_group_total_score_mean":sum(scores)/len(scores),
    "matched_group_total_score_median":sorted(scores)[len(scores)//2],
    "matched_group_score_second_highest":sorted(set(scores))[-2],
    "count_groups_at_least_observed":sum(s>=23 for s in scores),
    "finite_label_frequency":1/len(rows),
    "all_year_flowering_conditioned_groups":len(conditional),
    "all_year_flowering_conditioned_ge_observed":1,
    "all_year_conditioned_two_singleton_clones_necessarily_g":True,
    "leave_one_g_clone_out":drop_one,
    "causal_genotype_effect":None,
    "butterfly_breeding_survival_effect":None,
    "claim_boundary":"Exploratory nonrandomized labels, Kato et al already reported earlier flowering of provenance g. One label ranking is not inferential p or an independent ecological experiment."
}
print(json.dumps(summary,indent=2,ensure_ascii=False))
