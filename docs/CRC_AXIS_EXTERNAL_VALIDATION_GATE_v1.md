# CRC external microbiome-metabolite-host axis gate v1

## Frozen inputs

- Locked CRC species: *Peptostreptococcus stomatis*, *Parvimonas* unclassified, *Gemella morbillorum*, and *Fusobacterium nucleatum*.
- No substitute species is accepted as an exact validation. Genus-level or near-neighbour matches are sensitivity evidence only.
- Primary discovery cohort: KRCA paired tumour microbiome-metabolome samples.
- Held-out validation cohort: CRC paired tumour microbiome-metabolome samples.
- Primary microbial exposure: CLR abundance of the exact locked species.
- Metabolite family: all 55 measured metabolite variables; no outcome-guided subset.
- Covariates where estimable: age, sex, BMI, tumour indicator, and tumour side.
- Discovery gate: BH q < 0.10 across 55 metabolites.
- Validation gate: direction concordance and BH q < 0.10 across locked discovery candidates.
- Missing, discordant, or non-estimable evidence cannot be compensated by literature or database edges.

## Host-layer promotion rule

The host layer is tested only after a measured metabolite passes the discovery-validation bridge. A literature-supported receptor, enzyme, pathway, enrichment result, or single-cell localization cannot rescue a failed measured-metabolite bridge.

## Claim classes

- COMPLETE COMPUTATIONAL TRIANGULATION: species, measured metabolite, and host programme each pass their prespecified external gates.
- PARTIAL AXIS: one adjacent bridge passes but the complete three-layer chain does not.
- COMPONENT ONLY: a layer-specific result passes without an adjacent bridge.
- NO-GO: no measured metabolite passes the leakage-free bridge; host-axis promotion stops.

## Provenance boundary

Zenodo record 10.5281/zenodo.7326674 contains processed tumour microbiome and UPLC-MS metabolomics matrices. Association does not establish mediation, microbial production of a metabolite, or causality.
