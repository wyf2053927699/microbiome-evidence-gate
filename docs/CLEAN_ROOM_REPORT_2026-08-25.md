# GitHub v1.0.0 clean-room report

Date: 2026-08-25

Repository: https://github.com/wyf2053927699/microbiome-evidence-gate

The `verified_release` directory was copied to a new audit directory containing no prior outputs. A new Python 3.12 virtual environment was created and populated only from the pinned `requirements.txt`. The three public source inputs were staged afterward and verified against the frozen SHA-256 input manifest.

Results:

- CRC leakage-controlled discovery and held-out validation: PASS
- CRC raw-count whole-atlas donor pseudobulk: PASS
- Framework boundary benchmark: PASS
- Output files generated: 23
- Exact equality against canonical output hashes: 23/23
- Personal absolute paths detected in the publication tree: 0
- credential, token or private-key patterns detected: 0
- Third-party raw source files included in the publication tree: 0

Code is designated MIT. Author-generated derived data are designated CC BY 4.0. The v1.0.0 DOI remains pending until the reviewed GitHub Release is created and archived by Zenodo.
