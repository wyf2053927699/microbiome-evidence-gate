from __future__ import annotations

import csv
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "01_raw_data" / "second_scenario_candidates"
OUT = ROOT / "06_results" / "second_scenario_crc" / "host_validation"
OUT.mkdir(parents=True, exist_ok=True)

# Frozen mechanistic programs. These are literature-defined host response classes,
# not genes selected from either bulk or single-cell outcome matrix.
PROGRAMS = {
    "TLR4_NFKB_inflammatory_response": ["TLR4", "MYD88", "IRAK1", "TRAF6", "NFKB1", "RELA", "NFKBIA", "TNF", "IL6", "CXCL8", "PTGS2", "BIRC3"],
    "epithelial_EMT_barrier_response": ["CDH1", "OCLN", "TJP1", "CLDN1", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "MMP9"],
    "autophagy_stress_response": ["ATG5", "ATG7", "BECN1", "MAP1LC3B", "SQSTM1", "ULK1"],
}
HOUSEKEEPING = ["RPLP0", "RPL3", "RPL4", "RPL5", "RPL7", "RPL8", "RPL9", "RPL10", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL23", "RPL24", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35"]


def bh(vals):
    return multipletests(np.asarray(vals, float), method="fdr_bh")[1]


def parse_annotation(path):
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                break
        reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
        mapping = {}
        for row in reader:
            symbol = (row.get("Gene symbol") or "").split("///")[0].strip()
            if symbol:
                mapping[row["ID"]] = symbol
    return mapping


def parse_series(path):
    meta = {}
    table_lines = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                key, *vals = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
                meta.setdefault(key, []).append(vals)
    data = list(csv.reader(table_lines, delimiter="\t"))
    samples = [x.strip('"') for x in data[0][1:]]
    probes = [row[0].strip('"') for row in data[1:] if row]
    values = np.asarray([[float(x) for x in row[1:]] for row in data[1:] if row], dtype=float)
    expr = pd.DataFrame(values, index=probes, columns=samples)
    sample_meta = pd.DataFrame(index=samples)
    titles = meta.get("!Sample_title", [[]])[0]
    sources = meta.get("!Sample_source_name_ch1", [[]])[0]
    sample_meta["title"] = titles
    sample_meta["source"] = sources
    return expr, sample_meta


def gene_matrix(expr, mapping):
    genes = pd.Series([mapping.get(x) for x in expr.index], index=expr.index)
    keep = genes.notna()
    x = expr.loc[keep].copy()
    x["gene"] = genes[keep].values
    return x.groupby("gene").median(numeric_only=True)


def program_scores(genes):
    z = genes.sub(genes.mean(axis=1), axis=0).div(genes.std(axis=1).replace(0, np.nan), axis=0)
    scores = {}
    coverage = {}
    for name, members in PROGRAMS.items():
        present = [g for g in members if g in z.index]
        scores[name] = z.loc[present].mean(axis=0)
        coverage[name] = present
    return pd.DataFrame(scores), coverage


bulk_rows = []
bulk_coverage = {}
bulk_dir = RAW / "CRC_host_bulk"
for accession, platform in [("GSE74602", "GPL6104"), ("GSE39582", "GPL570")]:
    expr, sm = parse_series(bulk_dir / f"{accession}_series_matrix.txt.gz")
    genes = gene_matrix(expr, parse_annotation(bulk_dir / f"{platform}.annot.gz"))
    scores, coverage = program_scores(genes)
    bulk_coverage[accession] = coverage
    pvals = []
    temp = []
    if accession == "GSE74602":
        sm["condition"] = np.where(sm["source"].str.contains("Normal", case=False), "normal", "tumor")
        sm["patient"] = sm["source"].str.extract(r"([TN]\d+)")[0].str[1:]
        for program in PROGRAMS:
            frame = pd.DataFrame({"score": scores[program], "condition": sm["condition"], "patient": sm["patient"]})
            paired = frame.pivot(index="patient", columns="condition", values="score").dropna()
            diff = paired["tumor"] - paired["normal"]
            stat, p = wilcoxon(diff, alternative="two-sided", zero_method="wilcox")
            temp.append({"dataset": accession, "program": program, "design": "paired tumor-normal", "n_tumor": len(diff), "n_normal": len(diff), "effect": float(diff.mean()), "p": p})
            pvals.append(p)
    else:
        sm["condition"] = np.where(sm["source"].str.contains("non.tum|normal", case=False, regex=True), "normal", "tumor")
        for program in PROGRAMS:
            tumor = scores.loc[sm.index[sm["condition"] == "tumor"], program]
            normal = scores.loc[sm.index[sm["condition"] == "normal"], program]
            u, p = mannwhitneyu(tumor, normal, alternative="two-sided")
            effect = 2 * u / (len(tumor) * len(normal)) - 1
            temp.append({"dataset": accession, "program": program, "design": "independent tumor-normal", "n_tumor": len(tumor), "n_normal": len(normal), "effect": float(effect), "p": p})
            pvals.append(p)
    qvals = bh(pvals)
    for row, q in zip(temp, qvals):
        row["q_within_dataset"] = q
        bulk_rows.append(row)

bulk = pd.DataFrame(bulk_rows)
wide = bulk.pivot(index="program", columns="dataset", values=["effect", "q_within_dataset"])
replicated = []
for program in PROGRAMS:
    a = bulk[(bulk.dataset == "GSE74602") & (bulk.program == program)].iloc[0]
    b = bulk[(bulk.dataset == "GSE39582") & (bulk.program == program)].iloc[0]
    replicated.append({"program": program, "GSE74602_effect": a.effect, "GSE74602_q": a.q_within_dataset,
                       "GSE39582_effect": b.effect, "GSE39582_q": b.q_within_dataset,
                       "replicated_two_bulk_datasets": bool(a.q_within_dataset < 0.05 and b.q_within_dataset < 0.05 and np.sign(a.effect) == np.sign(b.effect))})


# Single-cell localisation: extract only frozen program and housekeeping genes
# from the dense GEO CSV, normalize program counts to a ribosomal housekeeping
# denominator, and perform donor-level exact sign-flip tests in the seven paired donors.
sc_dir = RAW / "CRC_single_cell"
annotation = pd.read_csv(sc_dir / "GSE200997_cell_annotation.csv.gz")
annotation = annotation.rename(columns={"Unnamed: 0": "cell"}).set_index("cell")
needed = set(HOUSEKEEPING)
for genes in PROGRAMS.values():
    needed.update(genes)
captured = {}
with gzip.open(sc_dir / "GSE200997_raw_UMI_count_matrix.csv.gz", "rt", newline="") as f:
    reader = csv.reader(f)
    header = next(reader)[1:]
    for row in reader:
        gene = row[0]
        if gene in needed:
            captured[gene] = np.fromiter((float(x) for x in row[1:]), dtype=float, count=len(header))

cells = pd.Index(header)
hk_present = [g for g in HOUSEKEEPING if g in captured]
denom = np.sum([captured[g] for g in hk_present], axis=0) + 1.0
sc_scores = {}
sc_coverage = {}
for program, genes in PROGRAMS.items():
    present = [g for g in genes if g in captured]
    sc_coverage[program] = present
    normalized = [np.log1p(1000.0 * captured[g] / denom) for g in present]
    sc_scores[program] = np.mean(normalized, axis=0)
score_df = pd.DataFrame(sc_scores, index=cells).join(annotation[["samples", "Condition"]], how="inner")
score_df["patient"] = score_df["samples"].str.extract(r"cac(\d+)")[0]
pseudo = score_df.groupby(["patient", "Condition"])[list(PROGRAMS)].mean().reset_index()
sc_rows = []
for program in PROGRAMS:
    paired = pseudo.pivot(index="patient", columns="Condition", values=program).dropna()
    diff = paired["Tumor"] - paired["Normal"]
    obs = abs(diff.mean())
    perm = []
    for signs in itertools.product([-1, 1], repeat=len(diff)):
        perm.append(abs(np.mean(diff.to_numpy() * np.asarray(signs))))
    p = (sum(x >= obs - 1e-12 for x in perm)) / len(perm)
    sc_rows.append({"program": program, "paired_donors": len(diff), "mean_tumor_minus_normal": float(diff.mean()), "exact_signflip_p": p})
sc = pd.DataFrame(sc_rows)
sc["exact_signflip_q"] = bh(sc["exact_signflip_p"])

bulk.to_csv(OUT / "bulk_program_tests.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(replicated).to_csv(OUT / "bulk_program_replication.csv", index=False, encoding="utf-8-sig")
pseudo.to_csv(OUT / "single_cell_donor_pseudobulk_scores.csv", index=False, encoding="utf-8-sig")
sc.to_csv(OUT / "single_cell_program_exact_signflip.csv", index=False, encoding="utf-8-sig")
(OUT / "program_coverage.json").write_text(json.dumps({"programs": PROGRAMS, "bulk": bulk_coverage, "single_cell": sc_coverage, "single_cell_normalizer": hk_present}, indent=2), encoding="utf-8")

summary = {"bulk_replicated_programs": int(pd.DataFrame(replicated)["replicated_two_bulk_datasets"].sum()),
           "single_cell_programs_q_lt_0_05": int((sc["exact_signflip_q"] < 0.05).sum()),
           "single_cell_paired_donors": int(sc["paired_donors"].min()),
           "single_cell_boundary": "donor-level tumor-versus-normal localisation; no cell-type-specific claim because the selected GEO annotation lacks a complete cell-type field"}
(OUT / "host_validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
