# Microbiome Evidence-Gate Framework

Version 1.0.0 release candidate supporting the manuscript **“An executable evidence-gate framework for claim calibration across microbiome multi-omics case studies.”**

This repository separates three claims that should not be conflated:

1. **Biological component evidence**: layer-specific results in public HCC/cirrhosis and CRC datasets.
2. **Framework behaviour**: exhaustive state-space, one-bit transition and coding-error benchmarks.
3. **External biological transportability**: not established by this release. The CRC analysis is a second disease deployment with an internal held-out split, not an independent population validation of universal accuracy.

## Repository structure

- `verified_release/`: path-independent, clean-room-tested runner for the leakage-controlled CRC analysis, whole-atlas donor pseudobulk and framework boundary benchmark.
- `analysis_scripts/`: analysis and figure-generation scripts used during manuscript development.
- `derived_data/`: author-generated summary tables and figure-source data; no third-party raw data are redistributed.
- `metadata/`: dataset registry and release file manifest.
- `docs/`: licences, third-party data restrictions and claim boundaries.

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

Citation metadata are provided in `CITATION.cff`. The version-specific permanent archive is https://doi.org/10.5281/zenodo.22089531.

## Status

`v1.0.0`: released on 25 August 2026 and permanently archived at https://doi.org/10.5281/zenodo.22089531.

