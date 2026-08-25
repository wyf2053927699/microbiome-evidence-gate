from __future__ import annotations

import csv
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "01_raw_data" / "second_scenario_candidates" / "CRC_single_cell"
OUT = ROOT / "06_results" / "second_scenario_crc_v2" / "host_validation"
OUT.mkdir(parents=True, exist_ok=True)

PROGRAMS = {
    "TLR4_NFKB_inflammatory_response": ["TLR4", "MYD88", "IRAK1", "TRAF6", "NFKB1", "RELA", "NFKBIA", "TNF", "IL6", "CXCL8", "PTGS2", "BIRC3"],
    "epithelial_EMT_barrier_response": ["CDH1", "OCLN", "TJP1", "CLDN1", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "MMP9"],
    "autophagy_stress_response": ["ATG5", "ATG7", "BECN1", "MAP1LC3B", "SQSTM1", "ULK1"],
}


def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


annotation = pd.read_csv(RAW / "GSE200997_cell_annotation.csv.gz")
annotation = annotation.rename(columns={"Unnamed: 0": "cell"}).set_index("cell")

matrix_path = RAW / "GSE200997_raw_UMI_count_matrix.csv.gz"
with gzip.open(matrix_path, "rt", newline="") as fh:
    header = next(csv.reader(fh))[1:]

cell_meta = annotation.reindex(header)
if cell_meta[["samples", "Condition"]].isna().any().any():
    raise RuntimeError("Raw-matrix cells do not map completely to GEO annotation")
cell_meta["patient"] = cell_meta["samples"].str.extract(r"cac(\d+)")[0]
cell_meta["unit"] = cell_meta["patient"] + "|" + cell_meta["Condition"]
units = sorted(cell_meta["unit"].unique())
unit_index = {u: i for i, u in enumerate(units)}
codes = cell_meta["unit"].map(unit_index).to_numpy(dtype=int)
cell_counts = np.bincount(codes, minlength=len(units)).astype(int)

needed = {g for members in PROGRAMS.values() for g in members}
captured = {}
library_totals = np.zeros(len(units), dtype=float)
genes_scanned = 0
with gzip.open(matrix_path, "rt", newline="") as fh:
    next(fh)
    for line in fh:
        comma = line.find(",")
        if comma < 0:
            continue
        gene = line[:comma].strip('"')
        values = np.fromstring(line[comma + 1 :], sep=",", dtype=float)
        if len(values) != len(header):
            raise RuntimeError(f"Malformed row {gene}: {len(values)} versus {len(header)} cells")
        aggregated = np.bincount(codes, weights=values, minlength=len(units))
        library_totals += aggregated
        if gene in needed:
            captured[gene] = aggregated
        genes_scanned += 1

coverage = {name: [g for g in members if g in captured] for name, members in PROGRAMS.items()}
if any(len(v) < 2 for v in coverage.values()):
    raise RuntimeError(f"Insufficient programme coverage: {coverage}")

genes = sorted(captured)
counts = pd.DataFrame({g: captured[g] for g in genes}, index=units).T
logcpm = np.log1p(counts.divide(library_totals, axis=1) * 1e6)
z = logcpm.sub(logcpm.mean(axis=1), axis=0).div(logcpm.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
scores = pd.DataFrame({name: z.loc[members].mean(axis=0) for name, members in coverage.items()})
scores.index.name = "unit"
scores = scores.reset_index()
scores[["patient", "Condition"]] = scores["unit"].str.split("|", expand=True)
scores["cell_count"] = cell_counts
scores["library_UMI"] = library_totals

test_rows = []
donor_diffs = []
for program in PROGRAMS:
    paired = scores.pivot(index="patient", columns="Condition", values=program).dropna()
    diff = paired["Tumor"] - paired["Normal"]
    observed = abs(float(diff.mean()))
    permuted = [abs(float(np.mean(diff.to_numpy() * np.asarray(signs)))) for signs in itertools.product([-1, 1], repeat=len(diff))]
    p = sum(x >= observed - 1e-12 for x in permuted) / len(permuted)
    test_rows.append({"program": program, "paired_donors": len(diff), "mean_tumor_minus_normal": float(diff.mean()),
                      "median_tumor_minus_normal": float(diff.median()), "exact_signflip_p": p,
                      "minimum_attainable_two_sided_p": 2 / (2 ** len(diff))})
    donor_diffs.extend({"program": program, "patient": patient, "tumor_minus_normal": value} for patient, value in diff.items())
tests = pd.DataFrame(test_rows)
tests["exact_signflip_q_three_programs"] = bh(tests["exact_signflip_p"])

scores.to_csv(OUT / "single_cell_whole_atlas_pseudobulk_scores.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(donor_diffs).to_csv(OUT / "single_cell_whole_atlas_donor_differences.csv", index=False, encoding="utf-8-sig")
tests.to_csv(OUT / "single_cell_whole_atlas_exact_signflip.csv", index=False, encoding="utf-8-sig")
summary = {
    "dataset": "GSE200997",
    "genes_scanned": genes_scanned,
    "cells": len(header),
    "annotated_cell_type_field_available": False,
    "aggregation": "raw UMI summed by donor and condition across the whole atlas",
    "normalization": "log1p CPM followed by gene-wise z score and frozen programme mean",
    "independent_unit": "paired donor",
    "paired_donors": int(tests["paired_donors"].min()),
    "programs_q_lt_0_05": int((tests["exact_signflip_q_three_programs"] < .05).sum()),
    "gate_state": "PARTIAL",
    "claim_class": "donor-level whole-atlas host support; not cell-type or spatial localisation",
    "composition_boundary": "cell-type composition cannot be separated because the supplied annotation has no cell-type field",
    "program_coverage": coverage,
}
(OUT / "single_cell_whole_atlas_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
print(tests.to_string(index=False))
