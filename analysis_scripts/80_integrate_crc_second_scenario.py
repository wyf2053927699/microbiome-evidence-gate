from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "08_manuscript" / "FULL_MANUSCRIPT_v2_7PLUS1.md"


def replace_between(text, start, end, replacement):
    pattern = re.escape(start) + r".*?(?=" + re.escape(end) + r")"
    out, n = re.subn(pattern, replacement.rstrip() + "\n\n", text, flags=re.S)
    if n != 1:
        raise RuntimeError(f"Expected one replacement for {start!r}; found {n}")
    return out


text = P.read_text(encoding="utf-8")
text = text.replace(
    "# An executable evidence-gate framework for claim calibration in microbiome multi-omics: an HCC case study",
    "# An executable evidence-gate framework for claim calibration and cross-disease portability in microbiome multi-omics",
    1,
)
text = re.sub(r"\*\*Main-text word count:\*\* \d+", "**Main-text word count:** [UPDATE AFTER FINAL TYPESETTING]", text, count=1)

abstract = """# Abstract

Cross-cohort microbiome multi-omics can generate persuasive biological networks even when their contributing layers do not jointly validate a mechanism. We developed an executable, non-compensatory evidence-gate framework that maps dataset roles, independent units, missing evidence, contradictions and replication status to bounded claim classes. Framework behaviour and biological components were first evaluated in a hepatocellular carcinoma (HCC) case study spanning nine public microbiome, metabolomic, bulk-transcriptomic and single-cell datasets. The same frozen architecture was then transported, without outcome-guided rule changes, to an independently screened colorectal cancer (CRC) scenario. In HCC, no microbial feature passed the complete discovery-validation-external sequence, no metabolite reproduced across plasma and liver, and all four proposed microbiota-metabolite-host families failed at least one essential gate, although three broad hepatic programmes reproduced directionally. In CRC, 277 participants with complete species, microbial-orthologue and faecal-metabolite matrices were split deterministically into 165 discovery and 112 validation participants. Five CRC-enriched species reproduced in the locked validation set and were retained in eight of nine threshold-sensitivity settings, whereas zero of 5,124 microbial orthologues and zero of 299 metabolites met the same two-stage rule. An autophagy/stress programme reproduced directionally in two independent bulk-tissue cohorts; a donor-aware single-cell EMT/barrier signal passed global correction but contradicted the bulk direction. Consequently, no complete CRC microbiota-metabolite-host axis was promoted. Across observed CRC cases and 24 single-dimension ablations, permissive rules promoted partial evidence paths, whereas the unchanged non-compensatory rule promoted none. The cross-disease execution therefore provides positive external evidence for operational portability and rule differentiation, while preserving a negative biological mechanism boundary. This framework separates component validation, framework-behaviour validation and external transportability, and supplies auditable inputs, stopping rules and claim classes for reproducibility-first microbiome data science.

**Keywords:** hepatocellular carcinoma; colorectal cancer; gut microbiome; metabolomics; multi-omics integration; reproducibility; evidence gating; single-cell RNA sequencing"""
text = replace_between(text, "# Abstract", "# Contribution to the field", abstract)

contribution = """# Contribution to the field

Microbiome studies commonly connect taxa, metabolites and host pathways measured in separate cohorts, creating a visually coherent network without showing that the complete chain is reproducible in people. We convert that interpretive problem into an executable data-science task. Dataset roles, independent units, multiplicity families, missing inputs, contradictions and claim classes are frozen before integration; an essential failed or missing link cannot be rescued by a high score elsewhere. The HCC application shows how the framework preserves qualified layer-specific findings while rejecting an unsupported gut-liver mechanism. More importantly, an outcome-blind screen selected CRC as a genuinely independent second disease scenario, and the unchanged rules were executed on participant-level microbial species, orthologues and metabolites, two bulk host cohorts and paired-donor single-cell data. CRC yielded five locked-validation microbial species and independently supported host components, yet no validated metabolite bridge. Permissive comparators promoted partial paths, whereas the non-compensatory framework remained stable through 24 ablations. The study therefore contributes both positive biological components and positive evidence of cross-disease operational portability, while demonstrating why portability of a decision framework must not be confused with validation of a causal mechanism."""
text = replace_between(text, "# Contribution to the field", "# Introduction", contribution)

crc_methods = """## Prospective second-disease screening and CRC transportability analysis

Before inspecting any candidate outcome matrix, we froze ten non-compensatory entry gates for a second disease scenario: public and lawful access; participant-level microbial abundance; microbial functional features; measured metabolites; resolvable participant identifiers; a disease-versus-reference contrast; sufficient metadata for the independent biological unit; an independent bulk host-expression cohort; a donor-resolved single-cell host dataset; and stable files that could be checksum-locked. CRC, inflammatory bowel disease, metabolic dysfunction-associated steatotic liver disease and type 2 diabetes were compared only on access, mapping, layer completeness and reproducibility. CRC scored 98/100 and was selected because every entry gate passed. HMP2/IBDMDB was retained as a documented rejected candidate because the required merged and versioned three-layer matrices were not all locked; no HMP2 outcome was inspected.

The CRC microbial analysis used the participant-level supplementary workbook of Yachida et al. (2019). We intersected MetaPhlAn2 species abundances, KEGG orthologue abundances and measured faecal metabolites, excluded multiple-polyp and post-surgical participants from the primary contrast, and retained 277 complete participants: 127 healthy controls and 150 CRC cases. A deterministic stage-stratified split assigned 165 participants to discovery and 112 to locked validation. Features present in at least 10% of participants were tested with two-sided Mann-Whitney tests separately by family. Benjamini-Hochberg correction was applied within species, orthologue and metabolite families. Reproduction required discovery q<0.05, validation P<0.05 in the same direction and full-cohort q<0.05. Nine prespecified sensitivity settings crossed prevalence thresholds of 5%, 10% and 20% with discovery FDR thresholds of 0.025, 0.05 and 0.10. Cross-layer links could be evaluated only between features that independently reproduced under these rules.

Host-component validation used three fixed programmes: TLR4/NF-kB inflammatory response, epithelial EMT/barrier response and autophagy/stress response. Programme scores were calculated in 30 paired tumour-normal samples from GSE74602 and independently in GSE39582; false-discovery correction was applied within each dataset, and replication required adjusted significance and direction agreement. GSE200997 was analysed as paired-donor single-cell pseudobulk after housekeeping normalisation. Exact two-sided sign-flip tests were applied across seven eligible paired donors and corrected across the three programmes. Because the supplied annotation did not support a complete frozen cell-type mapping, this analysis localised programmes to the paired single-cell dataset but did not make cell-type-specific claims.

For the transportability test, the unchanged eight-dimensional gate architecture was applied to three CRC evidence cases: an oral-anaerobe inflammatory path, a *Fusobacterium*-autophagy path and an oral-anaerobe EMT/barrier path. Any-positive, majority-dimension, additive-score and non-compensatory rules were compared on the same matrix. We then performed 24 single-dimension ablations, one for every case-by-dimension combination. The primary transportability outcome was defined prospectively as complete execution on an independent disease, production of bounded claim classes, observable differentiation from permissive rules and stability of the non-compensatory decision. It was not defined by obtaining a favourable biological direction and did not require promotion of a complete mechanism axis."""
text = replace_between(text, "HMP2/IBDMDB (Lloyd-Price et al., 2019) was prospectively assessed as a candidate second scenario under a", "## Comparative stress test of the evidence-gating framework", crc_methods + "\n\n## Comparative stress test of the evidence-gating framework")

crc_results = """## Cross-disease execution in CRC demonstrated operational portability and retained missing bridges

Outcome-blind screening selected CRC before the participant-level outcome matrices were inspected (Figure 7A). Of 277 participants with complete species, microbial-orthologue and faecal-metabolite matrices, 165 entered discovery and 112 entered locked validation. Five species were enriched in CRC and satisfied the complete two-stage reproduction rule: *Peptostreptococcus stomatis*, *Gemella morbillorum*, *Fusobacterium nucleatum*, *Parvimonas* unclassified and *Solobacterium moorei* (Figure 7B). All five were retained in eight of nine prevalence-by-FDR sensitivity settings. In contrast, zero of 5,124 tested microbial orthologues and zero of 299 tested metabolites reproduced under the same rule. No reproduced species-metabolite link was therefore eligible for promotion.

Host evidence was component-positive but internally bounded (Figure 7C,D). The autophagy/stress programme was lower in tumour tissue in both bulk cohorts (GSE74602 effect -0.283, BH q=0.0075; GSE39582 effect -0.612, q=1.72x10^-5). In paired-donor single-cell pseudobulk, the EMT/barrier programme was higher in tumour (mean paired difference 0.040; exact sign-flip BH q=0.0469). The corresponding bulk directions conflicted: -0.387 in GSE74602 (q=0.00046) and +0.443 in GSE39582 (q=0.00154). The TLR4/NF-kB programme did not reproduce. These results support CRC microbial and host components, not a participant-linked or causal axis.

The unchanged eight-dimensional architecture executed completely in CRC and returned bounded claim classes for all three prespecified evidence cases (Figure 7E). The any-positive rule promoted 3/3 partial cases, the majority rule promoted 1/3 and the additive rule promoted 0/3; the non-compensatory rule promoted 0/3 because every case lacked at least one essential bridge. Its decisions were unchanged in all 24 single-dimension ablations. The prespecified primary transportability outcome was therefore positive: the framework was executable in an independent disease, differentiated permissive rules and remained stable under ablation. The biological mechanism outcome remained negative because no complete microbiota-metabolite-host axis passed. This distinction is central: external transportability validates operational behaviour across diseases, not causal truth.

## Counterfactual audit of claim migration"""
text = replace_between(text, "## Framework benchmarking retained negative evidence and bounded claim classes", "To determine whether the gate architecture changed interpretation rather than", crc_results + "\n\nTo determine whether the gate architecture changed interpretation rather than")

text = text.replace(
    "It should be interpreted as a rigorously executed HCC case study of a\nreusable decision architecture, not as validation of portability to every\ndisease or as a new statistical estimator.",
    "The CRC execution extends this contribution from an HCC-only case study to a cross-disease operational validation. It demonstrates portability of the specified decision process, not universal portability, clinical accuracy or a new statistical estimator.",
)

old_lim = """External transportability also remains untested. HMP2/IBDMDB passed the
prospective metadata, biological-unit and catalogue-level layer checks, but the
pilot was stopped before outcome analysis because all three required versioned
matrices and participant mappings were not checksum-locked. This is a
framework-behaviour example of a MISSING/STOP decision, not validation in IBD
and not evidence that the framework generalises across diseases. Reopening the
pilot requires locking all prespecified matrices, identifiers, participant
mappings, checksums and multiplicity families before inspecting results."""
new_lim = """Cross-disease transportability was tested in one independently screened CRC scenario. The positive outcome is operational: unchanged gates could be executed, produced rule differentiation and remained stable under ablation. It does not establish universal disease portability or biological accuracy. The CRC microbial, metabolite, bulk-host and single-cell cohorts were not participant linked, and the absence of a reproduced metabolite bridge prevented a complete mechanism claim. HMP2/IBDMDB remains a documented MISSING/STOP candidate rather than an IBD result because the three required versioned matrices were not all checksum-locked before outcome analysis."""
if old_lim not in text:
    raise RuntimeError("Old limitation paragraph not found")
text = text.replace(old_lim, new_lim)

conclusion = """## Conclusion

An executable non-compensatory framework converted heterogeneous microbiome multi-omics into auditable claim classes in HCC and retained its prespecified behaviour in an independent CRC scenario. The CRC execution produced five locked-validation microbial species and reproducible host components, but no validated metabolite bridge and no complete microbiota-metabolite-host axis. The primary cross-disease transportability outcome was nevertheless positive because unchanged gates executed completely, differentiated permissive rules and remained stable across 24 ablations. The central contribution is therefore neither a universal biomarker nor a causal mechanism: it is a reproducibility-first method that preserves valid positive components while preventing incomplete evidence chains from being promoted beyond their support."""
text = replace_between(text, "## Conclusion", "# References", conclusion + "\n\n# References")

text = text.replace(
    "| Candidate HMP2 second scenario | Resource, metadata and three-layer catalogue feasibility | Feasibility gates passed | All versioned matrices, mappings and checksums locked before outcome inspection | MISSING/STOP; no IBD result or portability claim |",
    "| CRC external transportability scenario | Outcome-blind entry gates and complete participant-level three-layer matrices | 277 participants; 5 reproduced species; 1 replicated bulk programme; 1 single-cell programme | Unchanged eight-dimensional gates, rule comparison and 24 ablations | Positive operational portability; 0 complete mechanism axes |\n| Rejected HMP2 candidate | Resource, metadata and three-layer catalogue feasibility | Catalogue feasibility only | All versioned matrices, mappings and checksums locked before outcome inspection | MISSING/STOP; no IBD result |",
)
text = text.replace(
    "The constructed controls test decision-rule behaviour under deliberately incomplete evidence paths. They do not estimate biological accuracy, clinical performance or cross-disease portability.",
    "The constructed controls test decision-rule behaviour under deliberately incomplete evidence paths. CRC provides a separate external execution test; neither analysis estimates clinical accuracy or causal truth.",
)

fig7 = """## Figure 7. Outcome-blind CRC execution validates cross-disease framework portability while blocking an incomplete mechanism

**(A)**, Four candidate diseases were screened before outcome inspection using access, mapping, layer-completeness and reproducibility criteria. CRC passed all non-compensatory entry gates (98/100); 277 complete participants were assigned to a deterministic stage-stratified discovery set (n=165) and locked validation set (n=112). **(B)**, Rank-biserial CRC-versus-healthy effects for five species that passed discovery FDR, same-direction validation and full-cohort FDR. All five persisted in eight of nine threshold-sensitivity settings. **(C)**, Component audit: five of 217 species, zero of 5,124 microbial orthologues and zero of 299 metabolites reproduced; one of three bulk host programmes and one of three paired-donor single-cell programmes passed their component rules. A validated metabolite bridge is explicitly missing. **(D)**, Host-programme effects in paired GSE74602, GSE39582 and seven paired GSE200997 donors. Autophagy/stress reproduced in both bulk cohorts, whereas the single-cell EMT/barrier signal conflicted with the bulk directions. **(E)**, The unchanged eight-dimensional gate matrix applied to three CRC cases. Permissive rules promote partial evidence, whereas the non-compensatory rule promotes 0/3 and remains stable across 24 single-dimension ablations. This is positive operational transportability, not validation of a causal CRC mechanism."""
text = replace_between(text, "## Figure 7.", "## Supplementary Figure 1.", fig7 + "\n\n## Supplementary Figure 1.")

alt7 = """## Figure 7

Five-panel cross-disease validation figure. Outcome-blind screening selects CRC and locks 277 complete participants into discovery and validation sets. Five CRC-enriched species reproduce, whereas no microbial orthologue or metabolite does. Independent bulk and paired-donor single-cell datasets provide positive but partly contradictory host components. The unchanged non-compensatory gate rejects all three incomplete CRC evidence paths and is stable across 24 ablations, demonstrating operational portability without claiming a validated causal mechanism."""
text = replace_between(text, "## Figure 7\n", "## Supplementary Figure 1", alt7 + "\n\n## Supplementary Figure 1")

refs = """
Khaliq, A. M., Erdogan, C., Kurt, Z., et al. (2022). Refining colorectal cancer classification and clinical stratification through a single-cell atlas. *Genome Biol.* 23, 113. doi: 10.1186/s13059-022-02677-z

Marisa, L., de Reyniès, A., Duval, A., et al. (2013). Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. *PLoS Med.* 10, e1001453. doi: 10.1371/journal.pmed.1001453

Pek, M., Yatim, S. M. J. M., Chen, Y., et al. (2017). Oncogenic KRAS-associated gene signature defines co-targeting of CDK4/6 and MEK as a viable therapeutic strategy in colorectal cancer. *Oncogene* 36, 4975-4986. doi: 10.1038/onc.2017.120

Yachida, S., Mizutani, S., Shiroma, H., et al. (2019). Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the gut microbiota in colorectal cancer. *Nat. Med.* 25, 968-976. doi: 10.1038/s41591-019-0458-7
"""
text = text.replace("\n# Tables\n", refs + "\n# Tables\n", 1)

P.write_text(text, encoding="utf-8")
print(P)
