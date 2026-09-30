#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def rho(block: dict, key: str) -> float | None:
    value = block.get(key, {}).get("rho")
    return None if value is None else float(value)


def pval(block: dict, key: str) -> float | None:
    value = block.get(key, {}).get("p_two_sided")
    return None if value is None else float(value)


def all_expected_sign_loo(loo: dict) -> bool:
    if not loo:
        return False
    for block in loo.values():
        er = rho(block, "effective")
        dr = rho(block, "dominance")
        if er is None or dr is None or not (er > 0 and dr < 0):
            return False
    return True


def classify(result: dict, rule: dict) -> dict:
    gate = bool(result.get("evaluable_gate", {}).get("passed"))
    out = {
        "schema": "chocho_geb_realm_claim_promotion_result_v0.1",
        "rule_schema": rule.get("schema"),
        "realm_result_status": result.get("status"),
        "level": "not_evaluable",
        "checks": {},
    }
    if not gate:
        out["checks"]["evaluability_gate"] = False
        return out

    out["checks"]["evaluability_gate"] = True
    gen = result.get("generality", {})
    pooled = gen.get("pooled_within_realm_rank", {})
    eff = pooled.get("effective", {})
    dom = pooled.get("dominance", {})
    eff_r = eff.get("rho")
    dom_r = dom.get("rho")
    descriptive = (
        eff_r is not None and dom_r is not None
        and float(eff_r) > 0 and float(dom_r) < 0
    )
    out["checks"]["expected_primary_signs"] = descriptive
    if not descriptive:
        return out
    out["level"] = "descriptive_main_text"

    primary_sig = (
        eff.get("p_two_sided") is not None
        and dom.get("p_two_sided") is not None
        and float(eff["p_two_sided"]) < 0.05
        and float(dom["p_two_sided"]) < 0.05
    )
    fam = pooled.get("butterfly_family_stratified_permutation", {})
    fam_sign = (
        rho(fam, "effective") is not None
        and rho(fam, "dominance") is not None
        and rho(fam, "effective") > 0
        and rho(fam, "dominance") < 0
    )
    loo_ok = all_expected_sign_loo(pooled.get("leave_one_realm_out", {}))
    core = pooled.get("core_realm_sensitivity", {})
    core_ok = (
        rho(core, "effective") is not None
        and rho(core, "dominance") is not None
        and rho(core, "effective") > 0
        and rho(core, "dominance") < 0
    )
    cell = pooled.get("cell_weighted_realm_sensitivity", {})
    cell_ok = (
        rho(cell, "effective") is not None
        and rho(cell, "dominance") is not None
        and rho(cell, "effective") > 0
        and rho(cell, "dominance") < 0
    )
    out["checks"].update({
        "primary_architecture_p_lt_0_05_both": primary_sig,
        "family_stratified_expected_signs": fam_sign,
        "leave_one_realm_out_expected_signs": loo_ok,
        "core_realm_expected_signs": core_ok,
        "cell_weighted_expected_signs": cell_ok,
    })
    robust = all((primary_sig, fam_sign, loo_ok, core_ok, cell_ok))
    if not robust:
        return out
    out["level"] = "robust_main_text"

    fam_eff_p = pval(fam, "effective")
    fam_dom_p = pval(fam, "dominance")
    fam_abstract = (
        fam_eff_p is not None and fam_dom_p is not None
        and max(fam_eff_p, fam_dom_p) <= 0.10
        and min(fam_eff_p, fam_dom_p) < 0.05
    )
    boot = gen.get("realm_decoupling", {}).get(
        "architecture_vs_magnitude_bootstrap", {}
    )
    eff_ci = boot.get("effective_minus_magnitude_ci95")
    dom_ci = boot.get("dominance_oriented_minus_magnitude_ci95")
    eff_contrast = (
        isinstance(eff_ci, list) and len(eff_ci) == 2
        and float(eff_ci[0]) > 0
    )
    dom_contrast = (
        isinstance(dom_ci, list) and len(dom_ci) == 2
        and float(dom_ci[0]) > 0
    )
    out["checks"].update({
        "family_stratified_abstract_p_rule": fam_abstract,
        "effective_architecture_minus_magnitude_ci_above_zero": eff_contrast,
        "dominance_architecture_minus_magnitude_ci_above_zero": dom_contrast,
    })
    if all((fam_abstract, eff_contrast, dom_contrast)):
        out["level"] = "abstract_eligible"
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-json", type=Path, required=True)
    ap.add_argument("--promotion-rule", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    result = load(args.result_json)
    rule = load(args.promotion_rule)
    if rule.get("schema") != "chocho_geb_realm_claim_promotion_rule_v0.1":
        raise RuntimeError("unexpected promotion rule")
    if rule.get("status") != "FROZEN_BEFORE_REALM_GENERALITY_RESULT":
        raise RuntimeError("promotion rule is not frozen")
    payload = classify(result, rule)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
