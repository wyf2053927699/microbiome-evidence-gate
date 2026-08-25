from __future__ import annotations

import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "06_results" / "single_cell" / "paired_programs" / "prespecified_program_scores.csv"
OUT = ROOT / "06_results" / "single_cell" / "paired_programs" / "prespecified_program_exact_signflip.csv"

# Two-sided 0.975 Student-t quantiles for the paired sample sizes present here.
T975 = {5: 2.776445, 6: 2.570582, 7: 2.446912, 8: 2.364624, 9: 2.306004, 10: 2.262157}


with INPUT.open(encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.DictReader(handle))

values = defaultdict(dict)
meta = {}
for row in rows:
    key = (row["celltype"], row["program"], row["patient"])
    values[key][row["site"].lower()] = float(row["module_score"])
    meta[(row["celltype"], row["program"])] = row["genes_detected"]

tests = []
for celltype, program in sorted(meta):
    diffs = []
    donors = []
    for c, p, donor in sorted(values):
        if (c, p) != (celltype, program):
            continue
        sites = values[(c, p, donor)]
        if "tumor" in sites and "normal" in sites:
            donors.append(donor)
            diffs.append(sites["tumor"] - sites["normal"])
    n = len(diffs)
    observed = sum(diffs) / n
    permuted = []
    for signs in itertools.product((-1.0, 1.0), repeat=n):
        permuted.append(sum(s * d for s, d in zip(signs, diffs)) / n)
    p_exact = sum(abs(x) >= abs(observed) - 1e-15 for x in permuted) / len(permuted)
    sd = math.sqrt(sum((d - observed) ** 2 for d in diffs) / (n - 1))
    half = T975[n] * sd / math.sqrt(n)
    tests.append({
        "celltype": celltype,
        "program": program,
        "genes_detected": meta[(celltype, program)],
        "paired_n": n,
        "donors": ";".join(donors),
        "mean_tumor_minus_normal": f"{observed:.12g}",
        "mean_difference_ci95_low": f"{observed - half:.12g}",
        "mean_difference_ci95_high": f"{observed + half:.12g}",
        "exact_signflip_p_two_sided": f"{p_exact:.12g}",
        "permutations": len(permuted),
    })

# Benjamini-Hochberg over the same global 22-test programme family.
order = sorted(range(len(tests)), key=lambda i: float(tests[i]["exact_signflip_p_two_sided"]))
q = [1.0] * len(tests)
running = 1.0
for rank_rev, idx in reversed(list(enumerate(order, start=1))):
    raw = float(tests[idx]["exact_signflip_p_two_sided"])
    running = min(running, raw * len(tests) / rank_rev)
    q[idx] = min(1.0, running)
for idx, row in enumerate(tests):
    row["exact_signflip_q_global_bh"] = f"{q[idx]:.12g}"
    row["exact_signflip_global_fdr_pass"] = str(q[idx] < 0.05).lower()

with OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(tests[0]))
    writer.writeheader()
    writer.writerows(tests)

print(f"tests={len(tests)}")
print(f"global_exact_signflip_fdr_passes={sum(r['exact_signflip_global_fdr_pass']=='true' for r in tests)}")
print(OUT)
