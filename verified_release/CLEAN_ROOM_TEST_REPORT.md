# Clean-room execution report

Date: 2026-08-19

The release candidate was copied to a new directory containing no prior results. The three raw public inputs were staged only after the code copy. `run_release.ps1` verified all input SHA-256 values, then executed scripts 83, 84 and 85 using the pinned project environment.

Results:

- CRC leakage-free analysis: PASS
- CRC raw-count whole-atlas donor pseudobulk: PASS
- Framework boundary benchmark: PASS
- Output files generated: 23
- Hash equality against the canonical project outputs: 23/23
- Hidden dependence on a pre-existing result directory: not detected

The clean-room output manifest is `CLEAN_ROOM_TEST_SHA256.csv`. Public repository URL, archival DOI and licence remain author-controlled.
