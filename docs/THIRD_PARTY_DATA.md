# Third-party data acquisition and redistribution

This repository intentionally excludes third-party raw matrices and publisher supplements. Users must download them from the accessions or source publications listed in `metadata/DATASET_REGISTRY.csv` and comply with the corresponding repository or publisher terms.

The verified runner additionally requires exact filenames and SHA-256 values listed in `verified_release/INPUT_MANIFEST.csv`. Checksums identify the analysed versions but do not grant redistribution rights.

No literature PDF, GutMGene database export, credential, token, personal identifier or local absolute path is part of the release.

## Additional inputs for v1.1.0

- `01_raw_data/crc_axis_external/yachida_2019/Yachida_Supplementary_Tables_1-15.xlsx`: obtain from the supplementary material of Yachida et al. (2019). It is used for the paired taxon-metabolite analyses and is not redistributed.
- `01_raw_data/crc_axis_external/GSE241076/GSE241076_counts.txt.gz`: obtain from NCBI GEO accession GSE241076.
- `01_raw_data/crc_axis_external/GSE241076/GSE241076_FPKM.txt.gz`: obtain from NCBI GEO accession GSE241076.

The repository contains only author-generated analysis outputs from these resources.

