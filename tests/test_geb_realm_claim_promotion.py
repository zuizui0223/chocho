import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "classify_geb_realm_claim.py"
SPEC = importlib.util.spec_from_file_location("geb_realm_claim", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


RULE = {
    "schema": "chocho_geb_realm_claim_promotion_rule_v0.1",
    "status": "FROZEN_BEFORE_REALM_GENERALITY_RESULT",
}


def block(rho, p):
    return {"rho": rho, "p_two_sided": p}


def base_result():
    return {
        "status": "EXPLORATORY_POSTHOC_REALM_GENERALITY_AUDIT_COMPLETE",
        "evaluable_gate": {"passed": True},
        "generality": {
            "pooled_within_realm_rank": {
                "effective": block(0.45, 0.002),
                "dominance": block(-0.42, 0.003),
                "butterfly_family_stratified_permutation": {
                    "effective": block(0.45, 0.02),
                    "dominance": block(-0.42, 0.08),
                },
                "leave_one_realm_out": {
                    "A": {
                        "effective": block(0.43, 0.01),
                        "dominance": block(-0.40, 0.02),
                    },
                    "B": {
                        "effective": block(0.47, 0.01),
                        "dominance": block(-0.44, 0.02),
                    },
                },
                "core_realm_sensitivity": {
                    "effective": block(0.40, 0.01),
                    "dominance": block(-0.38, 0.02),
                },
                "cell_weighted_realm_sensitivity": {
                    "effective": block(0.39, 0.01),
                    "dominance": block(-0.36, 0.02),
                },
            },
            "realm_decoupling": {
                "architecture_vs_magnitude_bootstrap": {
                    "effective_minus_magnitude_ci95": [0.10, 0.50],
                    "dominance_oriented_minus_magnitude_ci95": [0.08, 0.46],
                }
            },
        },
    }


def test_gate_failure_is_not_evaluable():
    result = base_result()
    result["evaluable_gate"]["passed"] = False
    out = mod.classify(result, RULE)
    assert out["level"] == "not_evaluable"


def test_expected_signs_only_reach_descriptive_main_text():
    result = base_result()
    result["generality"]["pooled_within_realm_rank"]["effective"]["p_two_sided"] = 0.2
    out = mod.classify(result, RULE)
    assert out["level"] == "descriptive_main_text"


def test_robust_main_text_requires_all_directional_sensitivities():
    result = base_result()
    result["generality"]["pooled_within_realm_rank"][
        "cell_weighted_realm_sensitivity"
    ]["effective"]["rho"] = -0.05
    out = mod.classify(result, RULE)
    assert out["level"] == "descriptive_main_text"


def test_robust_but_not_abstract_when_bootstrap_contrast_crosses_zero():
    result = base_result()
    result["generality"]["realm_decoupling"][
        "architecture_vs_magnitude_bootstrap"
    ]["effective_minus_magnitude_ci95"] = [-0.02, 0.40]
    out = mod.classify(result, RULE)
    assert out["level"] == "robust_main_text"


def test_abstract_eligible_requires_family_permutation_and_both_contrasts():
    result = base_result()
    out = mod.classify(result, RULE)
    assert out["level"] == "abstract_eligible"
    assert out["checks"]["family_stratified_abstract_p_rule"] is True
    assert out["checks"][
        "effective_architecture_minus_magnitude_ci_above_zero"
    ] is True
    assert out["checks"][
        "dominance_architecture_minus_magnitude_ci_above_zero"
    ] is True
