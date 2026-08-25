from __future__ import annotations

import csv
import itertools
import json
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "06_results" / "axis_prioritisation" / "axis_dimension_scores.csv"
SPEC = ROOT / "00_protocol" / "EVIDENCE_GATE_OPERATIONAL_SPEC_v2.json"
OUT = ROOT / "06_results" / "framework_benchmark"
SEED = 20260818


def read_rows(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows):
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def additive_score(values):
    return sum(values) / len(values)


def decisions(values):
    return {
        "any_positive": int(any(v > 0 for v in values)),
        "majority_dimensions": int(sum(v >= 0.5 for v in values) >= 4),
        "additive_score_ge_0_5": int(additive_score(values) >= 0.5),
        "noncompensatory_all_complete": int(all(v == 1 for v in values)),
    }


with SPEC.open(encoding="utf-8") as handle:
    spec = json.load(handle)
dims = [item["id"] for item in spec["dimensions"]]
observed = read_rows(INPUT)

# Observed-candidate comparison. This is a status comparison, not a truth-labelled
# accuracy benchmark.
observed_rows = []
for row in observed:
    values = [float(row[d]) for d in dims]
    dec = decisions(values)
    observed_rows.append({
        "candidate_family": row["candidate_family"],
        "complete_dimensions": sum(v == 1 for v in values),
        "dimensions_at_least_partial": sum(v >= 0.5 for v in values),
        "additive_score": f"{additive_score(values):.6f}",
        **dec,
        "observed_gate_status": "eligible" if dec["noncompensatory_all_complete"] else "NO-GO",
        "interpretation": "workflow disagreement only; biological truth is not asserted",
    })
write_rows(OUT / "benchmark_observed_workflow_comparison.csv", observed_rows)

# Structural negative controls: the complete vector and every single-gate broken
# vector. A single broken essential link is known by construction and therefore
# permits a false-promotion stress test without claiming biological ground truth.
controls = []
for failed in [None, *dims]:
    values = [1.0] * len(dims)
    if failed is not None:
        values[dims.index(failed)] = 0.0
    dec = decisions(values)
    expected = int(failed is None)
    controls.append({
        "scenario": "all_complete" if failed is None else f"broken_{failed}",
        "known_broken_dimension": "none" if failed is None else failed,
        "expected_axis_eligibility": expected,
        **dec,
    })
write_rows(OUT / "benchmark_structural_negative_controls.csv", controls)

summary = []
methods = list(decisions([0] * len(dims)))
broken = [row for row in controls if row["known_broken_dimension"] != "none"]
complete = [row for row in controls if row["known_broken_dimension"] == "none"]
for method in methods:
    false_promotions = sum(int(row[method]) for row in broken)
    summary.append({
        "workflow": method,
        "broken_scenarios": len(broken),
        "false_promotions": false_promotions,
        "structural_false_promotion_rate": f"{false_promotions / len(broken):.6f}",
        "complete_scenario_retained": int(complete[0][method]),
        "scope": "synthetic structural stress test; not clinical accuracy",
    })
write_rows(OUT / "benchmark_structural_summary.csv", summary)

# Exhaustive binary perturbation across all 2^8 completeness patterns.
patterns = []
for bits in itertools.product([0.0, 1.0], repeat=len(dims)):
    dec = decisions(bits)
    patterns.append({
        "pattern": "".join(str(int(v)) for v in bits),
        "n_complete": int(sum(bits)),
        **dec,
    })
write_rows(OUT / "benchmark_all_binary_patterns.csv", patterns)

# Partial-code perturbation: vary every non-binary observed code within its
# adjacent admissible levels. This assesses rank/status sensitivity but never
# changes the all-complete promotion rule.
random.seed(SEED)
perturbed = []
levels = [0.0, 0.25, 0.5, 0.75, 1.0]
for row in observed:
    base = [float(row[d]) for d in dims]
    for iteration in range(1, 1001):
        vals = []
        for value in base:
            if value in (0.0, 1.0):
                vals.append(value)
            else:
                idx = levels.index(value)
                vals.append(random.choice(levels[max(0, idx - 1): min(len(levels), idx + 2)]))
        dec = decisions(vals)
        perturbed.append({
            "candidate_family": row["candidate_family"],
            "iteration": iteration,
            "additive_score": f"{additive_score(vals):.6f}",
            **dec,
        })
write_rows(OUT / "benchmark_partial_code_perturbations.csv", perturbed)

perturb_summary = []
for candidate in sorted({row["candidate_family"] for row in perturbed}):
    subset = [row for row in perturbed if row["candidate_family"] == candidate]
    scores = [float(row["additive_score"]) for row in subset]
    perturb_summary.append({
        "candidate_family": candidate,
        "iterations": len(subset),
        "additive_score_min": f"{min(scores):.6f}",
        "additive_score_max": f"{max(scores):.6f}",
        "additive_promotions": sum(int(row["additive_score_ge_0_5"]) for row in subset),
        "noncompensatory_promotions": sum(int(row["noncompensatory_all_complete"]) for row in subset),
        "interpretation": "coding-sensitivity audit; not an inferential confidence interval",
    })
write_rows(OUT / "benchmark_partial_code_summary.csv", perturb_summary)

with (OUT / "README.md").open("w", encoding="utf-8") as handle:
    handle.write("# Evidence-gate comparative stress test\n\n")
    handle.write("This directory compares permissive evidence-selection rules with the exact non-compensatory rule. ")
    handle.write("Synthetic broken-link controls provide structural truth only: an axis with one essential link set to zero must not be promoted. ")
    handle.write("The outputs do not estimate clinical accuracy, causality, or the biological truth of an HCC mechanism.\n")

print(f"Wrote framework benchmark outputs to {OUT}")
