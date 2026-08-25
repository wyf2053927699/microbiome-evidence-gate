# Microbiome Evidence-Gate Framework

Version 1.1.0 supporting the candidate manuscript **“Evidence-gated microbiome data integration identifies a triangulated colorectal cancer taxon-glycocholate-host response axis.”**

This repository separates three claims that should not be conflated:

1. **Biological component evidence**: layer-specific results in public HCC/cirrhosis and CRC datasets.
2. **Framework behaviour**: exhaustive state-space, one-bit transition and coding-error benchmarks.
3. **External biological transportability**: bounded, not universal. The CRC extension combines a paired participant-level taxon-metabolite analysis with an external GCA perturbation transcriptome, but does not establish microbial glycocholate production, same-participant mediation or causality.

## Repository structure

- `verified_release/`: path-independent, clean-room-tested runner for the leakage-controlled CRC analysis, whole-atlas donor pseudobulk and framework boundary benchmark.
- `analysis_scripts/`: analysis and figure-generation scripts used during manuscript development.
- `derived_data/`: author-generated summary tables and figure-source data; no third-party raw data are redistributed.
- `metadata/`: dataset registry and release file manifest.
- `docs/`: licences, third-party data restrictions and claim boundaries.

## Version 1.1.0 extension

This release adds a deterministic discovery/held-out analysis of 347 paired metagenome-metabolome profiles, alternative-profiler and disease-stratified sensitivity analyses, joint adjustment for four locked CRC taxa, and an external GSE241076 glycocholic-acid host-response audit. Source scripts for Figure 7, Figure 8 and Supplementary Figure 4 are included. The author-generated outputs are under `derived_data/crc_axis_v1_1/`.

The permitted positive claim is a **triangulated candidate axis** linking *Peptostreptococcus stomatis* abundance, glycocholate and a GCA-responsive SOX14-ZDHHC9-CD274 programme. It is not evidence that *P. stomatis* produces glycocholate and is not a same-participant causal mechanism.

The five Python analyses use `MICROBIOME_GATE_ROOT` when set and otherwise treat the repository root as the project root. Third-party inputs must be placed under `01_raw_data/crc_axis_external/` as described in `docs/THIRD_PARTY_DATA.md`; they are not included in this repository.

## Verified reproduction

Requirements: PowerShell 7, Python 3.12 or later, and the pinned packages in `verified_release/requirements.txt`.

1. Download the three source files listed in `verified_release/INPUT_MANIFEST.csv` from their original repositories or publisher pages.
2. Place them at the exact relative paths specified in that manifest.
3. Run:

```powershell
pwsh -File ./verified_release/run_release.ps1
```

The runner verifies every input SHA-256 before analysis and writes `RUN_SHA256.csv` for generated outputs. The frozen clean-room audit reproduced 23 of 23 canonical output hashes.

## Data availability and redistribution

Public third-party source data are not included. Their accession numbers and roles are recorded in `metadata/DATASET_REGISTRY.csv`. Users must obtain those files from the original repositories and comply with the original terms. Author-generated derived tables in `derived_data/` are released under CC BY 4.0; code is released under the MIT License.

## Citation

Citation metadata are provided in `CITATION.cff`. Version 1.0.0 is archived at https://doi.org/10.5281/zenodo.22089531. Zenodo will assign a separate version-specific DOI to v1.1.0 after this GitHub release is archived.

## Status

- `v1.0.0`: released on 25 August 2026 and permanently archived at https://doi.org/10.5281/zenodo.22089531.
- `v1.1.0`: adds the paired CRC taxon-metabolite analysis, external GCA host response, updated benchmark visualisation and bounded positive claim.

