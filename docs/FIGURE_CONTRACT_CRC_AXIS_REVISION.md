# Figure contract: CRC positive-axis revision

## Figure 7

Core conclusion: A leakage-free paired-cohort analysis narrows broad CRC taxon-metabolite covariation to a *Peptostreptococcus stomatis*-associated glycocholate signal that connects, across datasets, to a GCA-responsive SOX14-ZDHHC9-CD274 host programme.

Figure archetype: asymmetric mixed-modality figure.  
Target journal/output: Frontiers in Microbiology, double-column, 183 x 198 mm, editable SVG/PDF plus 600-dpi TIFF and PNG preview.  
Backend: R only.

Panel map:

- A: Yachida paired-cohort design, deterministic discovery/held-out split and locked taxa.
- B: discovery-versus-validation effects for all reproduced taxon-metabolite pairs, highlighting the six joint-taxon CRC-preferential candidates.
- C: joint four-taxon CRC-only coefficients and 95% confidence intervals for the six prioritised candidates.
- D: all-participant, CRC-only and healthy-only joint estimates for *P. stomatis*-glycocholate.
- E: GSE241076 sample PCA after GCA perturbation.
- F: sample-level expression of SOX14, ZDHHC9 and CD274 under GCA versus control.
- G: three-layer evidence chain with explicit claim-class labels on both links.

Evidence hierarchy:

- Hero evidence: panels C, D and F.
- Validation evidence: panels B and E.
- Design and interpretation boundary: panels A and G.

Statistics needed: discovery and locked-family BH q values; joint-model beta, standard error and BH q; gene-level Welch test with BH correction; exact programme permutation result in legend.

Reviewer risk: the microbial-metabolite link is observational; glycocholate is host-conjugated; cohorts are not participant-linked across all three layers; GSE241076 has n=3 per condition and treatment is inseparable from any unrecorded batch.

## Figure 8

Core conclusion: The same executable framework stops incomplete HCC mechanisms yet retains a bounded, triangulated CRC axis, and its non-compensatory decision boundary behaves predictably under exhaustive and coding-error stress tests.

Figure archetype: quantitative grid with a gate-matrix hero panel.  
Target journal/output: Frontiers in Microbiology, double-column, 183 x 158 mm, editable SVG/PDF plus 600-dpi TIFF and PNG preview.  
Backend: R only.

Panel map:

- A: cross-disease component/gate matrix, separating HCC STOP from CRC triangulated-candidate output.
- B: promotion frequency over all 256 binary evidence states for four decision rules.
- C: decision changes over 2,048 bidirectional one-bit transitions.
- D: false-promotion and false-stop behaviour across coding-error probabilities.
- E: nearest-rule/ablation comparison and operational interpretation.

Reviewer risk: benchmarks verify implementation and boundary behaviour, not biological truth, clinical accuracy or universal superiority.

## Supplementary Figure 4

Core conclusion: The CRC candidate survives prespecified and post-validation attribution audits, while the host perturbation dataset passes basic technical consistency checks.

Figure archetype: quantitative grid.  
Target journal/output: Frontiers supplementary figure, 183 x 205 mm.  
Backend: R only.

Panel map:

- A: candidate attrition from 1,008 tests to six joint-taxon CRC-preferential candidates.
- B: MetaPhlAn-versus-mOTU effect concordance for exact taxon matches.
- C: locked-taxon rank-correlation matrix.
- D: effect-direction heatmap across discovery, validation, mOTU, CRC-only, healthy-only and presence/absence audits.
- E: GSE241076 sample-correlation matrix.
- F: counts-versus-FPKM mechanism-gene effect concordance and library-size audit.

Reviewer risk: the 72-pair intermediate set is not a mechanism set; *Parvimonas micra* in mOTU is a near-neighbour sensitivity, not exact validation of MetaPhlAn *Parvimonas* unclassified.

