from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "06_results" / "framework_benchmark_v2"
CRC_OUT = ROOT / "06_results" / "second_scenario_crc_v2" / "framework_execution"
OUT.mkdir(parents=True, exist_ok=True)
CRC_OUT.mkdir(parents=True, exist_ok=True)

DIMS = [
    "microbial_reproducibility",
    "human_metabolite_observed",
    "microbe_metabolite_link",
    "metabolite_target_link",
    "host_reproduced_two_datasets",
    "single_cell_localisation",
    "major_contradiction_resolved",
    "robust_to_analytic_choices",
]
METHODS = ["any_positive", "majority_dimensions", "additive_score_ge_0_5", "noncompensatory_all_complete"]


def decisions(values):
    return {
        "any_positive": int(any(v > 0 for v in values)),
        "majority_dimensions": int(sum(v >= .5 for v in values) >= 4),
        "additive_score_ge_0_5": int(np.mean(values) >= .5),
        "noncompensatory_all_complete": int(all(v == 1 for v in values)),
    }


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


# Complete binary state space. The reference label is specification conformance,
# not biological truth: a mechanism-eligible claim requires all essential evidence.
patterns = []
for bits in itertools.product([0, 1], repeat=len(DIMS)):
    d = decisions(bits)
    patterns.append({"pattern": "".join(map(str, bits)), "n_complete": sum(bits),
                     "reference_specification_eligible": int(all(bits)), **d})
write_csv(OUT / "all_binary_patterns.csv", patterns)

# Every one-bit edit is evaluated in both directions. Unlike the previous
# one-way ablation, these transitions can cross the decision boundary.
transitions = []
for bits in itertools.product([0, 1], repeat=len(DIMS)):
    base = decisions(bits)
    for j, dim in enumerate(DIMS):
        changed = list(bits); changed[j] = 1 - changed[j]
        after = decisions(changed)
        row = {"from_pattern": "".join(map(str, bits)), "to_pattern": "".join(map(str, changed)),
               "dimension": dim, "direction": f"{bits[j]}->{changed[j]}",
               "reference_before": int(all(bits)), "reference_after": int(all(changed))}
        for method in METHODS:
            row[f"{method}_before"] = base[method]
            row[f"{method}_after"] = after[method]
            row[f"{method}_changed"] = int(base[method] != after[method])
        transitions.append(row)
write_csv(OUT / "bidirectional_one_bit_transitions.csv", transitions)

transition_summary = []
for method in METHODS:
    for direction in ("0->1", "1->0"):
        subset = [r for r in transitions if r["direction"] == direction]
        transition_summary.append({"workflow": method, "direction": direction, "transitions": len(subset),
                                   "decision_changes": sum(r[f"{method}_changed"] for r in subset),
                                   "decision_change_rate": sum(r[f"{method}_changed"] for r in subset) / len(subset)})
write_csv(OUT / "bidirectional_transition_summary.csv", transition_summary)

# Boundary-focused Monte Carlo coding-error benchmark. Results quantify a
# conservative-versus-permissive trade-off against the declared specification;
# they are not clinical accuracy or biological validation.
rng = np.random.default_rng(20260819)
reference_cases = [("complete_positive", np.ones(len(DIMS), dtype=int), 1)]
for j, dim in enumerate(DIMS):
    x = np.ones(len(DIMS), dtype=int); x[j] = 0
    reference_cases.append((f"single_missing__{dim}", x, 0))

simulation = []
for error_rate in (.01, .05, .10, .20):
    for case_name, truth, reference in reference_cases:
        draws = np.tile(truth, (10000, 1))
        flips = rng.random(draws.shape) < error_rate
        observed = np.where(flips, 1 - draws, draws)
        for method in METHODS:
            pred = np.array([decisions(row)[method] for row in observed])
            simulation.append({"case": case_name, "reference_specification_eligible": reference,
                               "per_dimension_coding_error": error_rate, "workflow": method,
                               "simulations": len(pred), "promotion_rate": float(pred.mean()),
                               "false_promotion_rate": float(pred.mean()) if reference == 0 else np.nan,
                               "false_stop_rate": float(1 - pred.mean()) if reference == 1 else np.nan})
write_csv(OUT / "coding_error_simulation.csv", simulation)

simulation_summary = []
for error_rate in (.01, .05, .10, .20):
    for method in METHODS:
        sub = [r for r in simulation if r["per_dimension_coding_error"] == error_rate and r["workflow"] == method]
        fp = [r["false_promotion_rate"] for r in sub if r["reference_specification_eligible"] == 0]
        fs = [r["false_stop_rate"] for r in sub if r["reference_specification_eligible"] == 1]
        simulation_summary.append({"per_dimension_coding_error": error_rate, "workflow": method,
                                   "mean_false_promotion_near_boundary": float(np.mean(fp)),
                                   "false_stop_complete_case": float(fs[0]),
                                   "interpretation": "specification-conformance trade-off; not clinical accuracy"})
write_csv(OUT / "coding_error_summary.csv", simulation_summary)

# Partial-code boundary cases explicitly test uncertainty rather than silently
# converting PARTIAL to PASS or FAIL.
partial_rows = []
for j, dim in enumerate(DIMS):
    values = np.ones(len(DIMS)); values[j] = .5
    partial_rows.append({"case": f"single_partial__{dim}", "partial_dimension": dim,
                         "reference_claim_class": "uncertain-not-eligible", **decisions(values)})
write_csv(OUT / "partial_code_boundary_cases.csv", partial_rows)

# Updated CRC execution. Four species pass leakage-free internal holdout;
# whole-atlas single-cell support is PARTIAL and cannot satisfy localisation.
crc_cases = [
    {"candidate_family": "Oral-anaerobe inflammatory-response case", "microbial_reproducibility": 1.0,
     "human_metabolite_observed": 0.0, "microbe_metabolite_link": 0.0, "metabolite_target_link": 0.0,
     "host_reproduced_two_datasets": 0.0, "single_cell_localisation": 0.5,
     "major_contradiction_resolved": 0.0, "robust_to_analytic_choices": 1.0},
    {"candidate_family": "Fusobacterium-autophagy case", "microbial_reproducibility": 1.0,
     "human_metabolite_observed": 0.0, "microbe_metabolite_link": 0.0, "metabolite_target_link": 0.0,
     "host_reproduced_two_datasets": 1.0, "single_cell_localisation": 0.5,
     "major_contradiction_resolved": 1.0, "robust_to_analytic_choices": 1.0},
    {"candidate_family": "Oral-anaerobe EMT-barrier case", "microbial_reproducibility": 1.0,
     "human_metabolite_observed": 0.0, "microbe_metabolite_link": 0.0, "metabolite_target_link": 0.0,
     "host_reproduced_two_datasets": 0.0, "single_cell_localisation": 0.5,
     "major_contradiction_resolved": 0.0, "robust_to_analytic_choices": 1.0},
]
crc_rows = []
for case in crc_cases:
    vals = [case[d] for d in DIMS]
    crc_rows.append({**case, "complete_dimensions": sum(v == 1 for v in vals),
                     "partial_dimensions": sum(v == .5 for v in vals), **decisions(vals),
                     "claim_class": "mechanism-eligible" if all(v == 1 for v in vals) else "NO-GO"})
write_csv(CRC_OUT / "crc_candidate_gate_matrix_v2.csv", crc_rows)

comparison = [
    {"dimension": "Primary aim", "present_framework": "claim calibration under missing or contradictory multi-omics evidence",
     "Ascandari_2026": "causal identification and precision-intervention workflow", "difference": "different decision target"},
    {"dimension": "Core operations", "present_framework": "non-compensatory gates and bounded claim classes",
     "Ascandari_2026": "DAG, Mendelian randomisation, double machine learning, mediation and reversibility",
     "difference": "evidence eligibility versus causal estimation"},
    {"dimension": "Missing evidence", "present_framework": "explicit MISSING or PARTIAL state blocks mechanism eligibility",
     "Ascandari_2026": "assumption guardrails constrain progression", "difference": "present framework exposes missing-link states as outputs"},
    {"dimension": "Disease demonstration", "present_framework": "HCC plus CRC deployment stress test",
     "Ascandari_2026": "CRC public-data operational demonstration", "difference": "cross-disease deployment but not validated portability"},
    {"dimension": "Validation status", "present_framework": "specification-conformance and coding-error trade-off only",
     "Ascandari_2026": "operational causal-workflow demonstration", "difference": "neither supplies universal biological truth"},
]
write_csv(OUT / "nearest_framework_comparison.csv", comparison)

manifest = {
    "benchmark_version": "v2 boundary and coding-error stress test",
    "previous_floor_effect_removed": True,
    "bidirectional_transitions": len(transitions),
    "coding_error_simulations": sum(r["simulations"] for r in simulation),
    "external_truth_available": False,
    "inter_rater_agreement": "UNMEASURED; requires independent human coders",
    "permitted_conclusion": "framework executability, decision-boundary behaviour and specification-conformance trade-offs",
    "forbidden_conclusion": "validated clinical accuracy, biological truth or universal cross-disease portability",
    "crc_complete_axes": sum(r["noncompensatory_all_complete"] for r in crc_rows),
}
(OUT / "benchmark_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
(CRC_OUT / "crc_execution_manifest_v2.json").write_text(json.dumps({
    "execution_status": "COMPLETE", "complete_axes": 0,
    "single_cell_localisation_state": "PARTIAL",
    "portability_claim": "cross-disease deployment stress test; not external validation",
}, indent=2), encoding="utf-8")
print(json.dumps(manifest, indent=2))
