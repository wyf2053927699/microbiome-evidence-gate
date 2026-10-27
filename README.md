# Microbiome Evidence-Gate Framework

## October revision v1.2.0

The October revision code and analysis-ready inputs are supplied in [microbiome-evidence-gate-v1.2.0.zip](microbiome-evidence-gate-v1.2.0.zip). This standalone package is the reference for the revised manuscript results. The earlier directories preserve the August workflow.

Verify [SHA256SUMS.txt](SHA256SUMS.txt), extract the ZIP, install its requirements.txt and the listed R dependencies, and run from the extracted package root:

```sh
python -B -X utf8 run_revision.py --rscript Rscript
```

The package includes corrected numeric covariate encoding, rank-adjusted residual-correlation inference, rerun selection and BH procedures, unified-partition sensitivity with global BH, deposited-count TMM/limma-voom analysis, independent-rating checks, target-panel provenance, and figure source workflows. Its README identifies inputs, commands, dependencies and output roles. INDEPENDENT_REPRODUCTION_20261027.json records successful comparisons of ten regenerated key tables and independent-rating/global-BH verification.

See [RELEASE_NOTES_v1.2.0.md](RELEASE_NOTES_v1.2.0.md) for the version history.

Version-specific archive for v1.2.0: https://doi.org/10.5281/zenodo.23682888 .

## Historical workflow and archives

The original August workflow remains under the immutable v1.1.0 tag and its [version-specific Zenodo archive](https://doi.org/10.5281/zenodo.22092200). The [version-series DOI](https://doi.org/10.5281/zenodo.22089530) identifies the release series. Version v1.0.0 remains at https://doi.org/10.5281/zenodo.22089531. The August archived version predates the October corrections and is not the reference for the revised analyses.

Historical directories include analysis_scripts, derived_data, metadata, docs, and verified_release. Their original runners and outputs retain their historical scope. Use the standalone v1.2.0 package for the October revision.

## Data and licensing

Code is licensed under MIT; author-generated derived data are licensed under CC BY 4.0. Public analysis-ready source inputs in the October package are attributed in its README and retain the original source reuse conditions. Third-party material is excluded from the author-generated data licence. Original upstream processing is documented in the archived workflow and source repositories; the October package reproduces revision analyses from supplied source or analysis-ready tables.

