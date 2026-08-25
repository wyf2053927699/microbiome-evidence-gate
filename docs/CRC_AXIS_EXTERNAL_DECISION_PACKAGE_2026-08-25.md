# CRC external microbiome-metabolite-host decision package

Date: 2026-08-25  
Status: independent analysis package; manuscript and approved figures remain frozen.

## Decision

**Conditional GO** for a triangulated CRC taxon-associated metabolite-host response axis, centred on *Peptostreptococcus stomatis*--glycocholate--GCA-responsive SOX14/ZDHHC9/PD-L1 programme.

This is not a same-participant causal chain and must not be described as proof that *P. stomatis* produces glycocholate. Glycocholate is a host-conjugated primary bile acid; the microbial result is an adjusted association and may reflect ecological context, deconjugation, host physiology, or residual confounding.

## Data provenance

- Yachida et al. paired metagenome-metabolome workbook: 347 intersecting participants, matching the published cohort; SHA-256 `A43766D21A09127563C9F319A10800629004FC55D36C1DAB90EA05B87A6B7FB3`.
- GSE241076 official GEO counts: 204,489 bytes; SHA-256 `FB32351613F0EEE2AF159BC35B5B34FEF5C6D7A678071D7E9EC99C87074CBDB6`.
- GSE241076 official GEO FPKM: 803,306 bytes; SHA-256 `F749A6ECB0E8297704B53C89A92C6F3EFC1761EBFA5DDCA0FE51060B9868C002`.

## Frozen primary analysis

- Deterministic discovery/held-out split: 207/140 participants.
- Four locked taxa, 252 discovery-filtered metabolites, 1,008 tests.
- Clinical adjustment: disease group, age, sex, BMI, smoking, alcohol and tumour location where available.
- Discovery BH q < 0.10 within each locked taxon; held-out candidate-family BH q < 0.10 plus directional agreement.
- Outcome: 223 discovery candidates and 96 held-out reproduced associations. The breadth triggered additional attribution audits and was not treated as 96 mechanisms.

## Robustness and attribution audits

- Alternative mOTU profiler, CRC-only models, presence/absence exposure and direction checks reduced neither the broad shared gradient sufficiently nor taxon attribution ambiguity: 72 pairs passed the initial combined sensitivity rule.
- Joint four-taxon models then controlled oral-taxon co-occurrence. Six associations retained q < 0.10 in the complete cohort and CRC-only subset, had the same direction as discovery/validation, were stronger in CRC, and were not supported in healthy-only models.
- These six are post-validation prioritisation results, not a replacement confirmatory gate.

Top candidates:

1. *P. stomatis*--beta-Ala-Lys: joint CRC beta 0.634, q 0.035; limited host-mechanism tractability.
2. *F. nucleatum*--3-aminopropane-1,2-diol: joint CRC beta 0.457, q 0.051; strong taxon attribution but weak public host-programme support.
3. *P. stomatis*--glycocholate: joint complete-cohort beta 0.281, q 0.020; joint CRC beta 0.485, q 0.051; healthy beta 0.055, q 0.879. This is the best three-layer candidate because a direct GCA perturbation transcriptome is publicly available.
4. *P. stomatis*--alanine: joint CRC beta 0.570, q 0.051.
5. *P. stomatis*--2-hydroxybutyrate: joint CRC beta 0.434, q 0.081.
6. *P. stomatis*--pyridoxamine 5'-phosphate: joint CRC beta -0.438, q 0.091.

## External host perturbation layer

GSE241076 contains MC38 CRC cells treated with GCA versus control (n = 3 per group). After log2(CPM + 0.5) normalisation and filtering to 12,655 expressed genes:

- SOX14: log2FC 8.811, BH q 0.0147.
- ZDHHC9: log2FC 5.310, BH q 0.0127.
- CD274/PD-L1: log2FC 7.100, BH q 0.0246.
- The locked three-gene programme increased by 1.822 standardised-score units. Its exact two-sided label-permutation p is 0.10; with 3 versus 3 samples, 0.10 is the smallest attainable two-sided exact value.

Quality checks:

- Library sizes: 18.06--20.53 million reads.
- Median within-group logCPM correlation: 0.991.
- Median between-group correlation: 0.966.
- PC1 explained 64.9% of variance.
- Counts and FPKM effect estimates were concordant for all three mechanism genes.

Condition and any unrecorded batch remain inseparable in this six-sample experiment. The public article supplies orthogonal wet-lab support, but this reanalysis alone is not independent biological replication of that article.

## Claim classes

Allowed:

- reproducible, covariate-adjusted association between *P. stomatis* abundance and faecal glycocholate in the Yachida cohort;
- association survives alternative profiling, CRC-only, presence/absence and joint-taxon attribution audits;
- GCA perturbation is linked to a SOX14/ZDHHC9/CD274 host-response programme in an external experimental dataset;
- the combined result is a cross-dataset, triangulated candidate axis.

Not allowed:

- *P. stomatis* produces glycocholate;
- the taxon causes the host programme;
- the three layers were measured in the same participants;
- the axis is diagnostic, prognostic, therapeutic or causal;
- the other 71 robust pairs are independent mechanisms.

## Single-cell and further-dataset decision

GSE241415 can localise GCA-associated immune changes in a mouse tumour experiment, but it cannot close the microbial-production gap. It should be added only if cell-type localisation is required for the final figure contract. The original route reserved the Gut paired cohort for a negative Yachida result; Yachida is positive, so that fallback is not mandatory. An additional independent paired cohort would strengthen portability but is not required to decide this axis.

## Approval gate

Recommended approval: incorporate the axis as **triangulated external validation**, not as a complete causal microbial-metabolite-host mechanism. On approval, manuscript, figures, legends, tables, supplements, cover letter and scope statement can be synchronised in one controlled revision. Until approval, none of those frozen files should change.

## Audit locations

- Primary: `06_results/crc_axis_external_yachida2019/`
- Robustness: `06_results/crc_axis_external_yachida2019/robustness/`
- Joint attribution: `06_results/crc_axis_external_yachida2019/independence_audit/`
- Host programme: `06_results/crc_axis_external_yachida2019/host_program_GSE241076/`
- Freeze manifest: `09_audit/crc_axis_external_freeze_2026-08-25.csv`
- Freeze recheck: 45 files checked; 0 changed or missing.
