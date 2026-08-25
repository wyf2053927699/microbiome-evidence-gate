from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "06_results" / "second_scenario_crc" / "framework_portability"
OUT.mkdir(parents=True, exist_ok=True)

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

# These are transparent framework test cases assembled from the frozen CRC
# component families. They are not asserted to be complete causal axes.
CASES = [
    {
        "candidate_family": "Oral-anaerobe inflammatory-response test case",
        "microbial_reproducibility": 1.0,
        "human_metabolite_observed": 0.0,
        "microbe_metabolite_link": 0.0,
        "metabolite_target_link": 0.0,
        "host_reproduced_two_datasets": 0.0,
        "single_cell_localisation": 0.75,
        "major_contradiction_resolved": 0.0,
        "robust_to_analytic_choices": 0.75,
        "evidence_note": "Five species reproduced; TLR4/NF-kB bulk replication failed; single-cell signal did not survive BH correction.",
    },
    {
        "candidate_family": "Fusobacterium-autophagy stress test case",
        "microbial_reproducibility": 1.0,
        "human_metabolite_observed": 0.0,
        "microbe_metabolite_link": 0.0,
        "metabolite_target_link": 0.0,
        "host_reproduced_two_datasets": 1.0,
        "single_cell_localisation": 0.0,
        "major_contradiction_resolved": 1.0,
        "robust_to_analytic_choices": 0.75,
        "evidence_note": "F. nucleatum reproduced and autophagy/stress bulk program reproduced in two cohorts; metabolite bridge and single-cell localisation were absent.",
    },
    {
        "candidate_family": "Oral-anaerobe EMT-barrier test case",
        "microbial_reproducibility": 1.0,
        "human_metabolite_observed": 0.0,
        "microbe_metabolite_link": 0.0,
        "metabolite_target_link": 0.0,
        "host_reproduced_two_datasets": 0.0,
        "single_cell_localisation": 1.0,
        "major_contradiction_resolved": 0.0,
        "robust_to_analytic_choices": 0.5,
        "evidence_note": "Donor-level single-cell localisation passed, but the two bulk datasets had significant opposite directions.",
    },
]


def decisions(values):
    return {
        "any_positive": int(any(v > 0 for v in values)),
        "majority_dimensions": int(sum(v >= 0.5 for v in values) >= 4),
        "additive_score_ge_0_5": int(sum(values) / len(values) >= 0.5),
        "noncompensatory_all_complete": int(all(v == 1 for v in values)),
    }


rows = []
for case in CASES:
    vals = [case[d] for d in DIMS]
    rows.append({**case, "complete_dimensions": sum(v == 1 for v in vals),
                 "additive_score": sum(vals) / len(vals), **decisions(vals),
                 "claim_class": "eligible" if all(v == 1 for v in vals) else "NO-GO"})

with (OUT / "crc_candidate_gate_matrix.csv").open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

# Disease-grounded ablation: remove each dimension in turn from each observed
# test case and quantify permissive false promotions relative to the unchanged
# non-compensatory decision.
ablations = []
for case in rows:
    base = [case[d] for d in DIMS]
    for dim in DIMS:
        vals = list(base)
        vals[DIMS.index(dim)] = 0.0
        ablations.append({"candidate_family": case["candidate_family"], "ablated_dimension": dim,
                          **decisions(vals), "additive_score": sum(vals) / len(vals)})
with (OUT / "crc_rule_ablation.csv").open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(ablations[0]))
    w.writeheader(); w.writerows(ablations)

summary = []
for method in ["any_positive", "majority_dimensions", "additive_score_ge_0_5", "noncompensatory_all_complete"]:
    summary.append({"workflow": method,
                    "observed_test_cases_promoted": sum(r[method] for r in rows),
                    "observed_test_cases_total": len(rows),
                    "ablation_scenarios_promoted": sum(r[method] for r in ablations),
                    "ablation_scenarios_total": len(ablations),
                    "interpretation": "framework-behaviour comparison; not clinical accuracy"})
with (OUT / "crc_rule_comparison_summary.csv").open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0]))
    w.writeheader(); w.writerows(summary)

primary_positive = (
    all(r["noncompensatory_all_complete"] == 0 for r in rows)
    and any(r["any_positive"] != r["noncompensatory_all_complete"] or r["majority_dimensions"] != r["noncompensatory_all_complete"] or r["additive_score_ge_0_5"] != r["noncompensatory_all_complete"] for r in rows)
    and all(r["noncompensatory_all_complete"] == 0 for r in ablations)
)
manifest = {
    "framework_portability_primary_outcome": "POSITIVE" if primary_positive else "NOT POSITIVE",
    "basis": "unchanged rules executed; observed CRC evidence produced auditable claim classes; permissive and non-compensatory workflows disagreed; non-compensatory status remained stable under every single-dimension ablation",
    "biological_component_positive": True,
    "complete_axis_positive": False,
    "causal_or_therapeutic_claim_permitted": False,
}
(OUT / "crc_framework_portability_outcome.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps(manifest, indent=2))
