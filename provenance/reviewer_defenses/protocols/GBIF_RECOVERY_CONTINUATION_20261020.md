# GBIF recovery continuation — 20 October 2026

## Purpose and fixed ecological question

Complete the already specified 239-butterfly macro comparison of sampled
butterfly occurrences with introduced-host resource geography. Keep the
original 135 completed species' observations unchanged. Recover only the 96
transport-incomplete species using the original two 300-record ordinal
windows per time period (2010–2017 and 2018–2025). The remaining eight
taxonomically unresolved names stay unresolved in this continuation.

The primary ecological outcome, climate thresholds, effort gates, seeds and
72 contrasts remain those of the 19 October script. Report three cohorts:
original135 targets and effort; expanded targets with original135 effort;
expanded targets and expanded effort. Do not select the headline by sign.

## Source correction identified before the new outcome

Run 39246353812's four successful *jobs* yielded only one completed source
taxon, zero mapped records and 95 incomplete species. The lone completed
query used the wrong key for Chlosyne lacinia: recovery CSV 1914469 versus
the original source receipt's 1914465. The merge correctly rejected this.
Correct that single transcribed key, and validate every recovery name/key
against the unchanged original source receipt before any new acquisition.
Do not weaken the merge check or resolve a different live taxon to pass it.

Original source artifact 12178748642 (run 39233336804) ZIP SHA256:
`a789fec816772c4eb5f24140e8481fd26b5cac49df61cb3d5a241cf8cdb8d23f`.
Original receipt SHA256:
`18f84371cf2a9a71ca06202fac2246d9aa3df0465bd973d06b4e9e8bfa58e8d2`.
Corrected 96-key CSV SHA256:
`d38e117da0053fe15c6a1d049fee899aabac0eb8290a9be0616a8e3507049c45`.

The provider index is live. Equal query rules and ordinal budgets do not
guarantee equal record identities or index offsets between acquisition dates.
Retain and report the existing per-window count/offset-drift audit.

## Information-size gate clarified before inspecting recovered outcomes

The 19 October written protocol required at least 50 informative strata but
did not define informative precisely. Its implementation used at least ten
strata with K>0, where K is the number of sampled positive candidate regions.
This mismatch must be explicit. The old effect estimates and legacy status
remain reproducible; they are not retrospectively called protocol-qualified.

Add a separate conservative gate counting strata with **0<A<N and 0<K<N**,
where N is the candidate-region count and A the added-host-region count.
These are exactly strata with nonzero conditional hypergeometric variance.
Require **50** such strata, as in the written protocol. Also retain the counts
of all mixed-resource strata, positive strata under the old definition, K=0
strata and K=N strata. Passing this sample-size check is not confirmation of
a causal association, sampling representativeness or biological absence.

The source-identical original135 replay gives 276 mixed-resource strata,
30 with K>0, of which one has K=N, hence **29** with nonzero conditional
variance. Its original +3.6819645 pp estimate and genus interval
[-0.9931785,+9.1044238] pp are unchanged and remain exploratory.

## Boundaries

No additional plant-range null variants, individual published-experiment
reanalysis, inferred colonisation, true absence, larval fitness or recruitment.
Native-host resource geography is not butterfly-native geography. Missing
host records and region-level geography can explain occurrences outside
native-host support. Only derived GBIF regional aggregates, receipts, code
and analysis summaries may enter this branch; keep author-only WCVP/HOSTS/
TDWG source files outside the repository and any public review artifact.
