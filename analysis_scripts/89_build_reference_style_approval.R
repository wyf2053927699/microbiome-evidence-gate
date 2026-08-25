options(stringsAsFactors = FALSE)

root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out <- file.path(root, "07_figures", "reference_style_approval_2026-08-19")
scratch <- file.path(out, "_scratch")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)

font <- "Times New Roman"
if (.Platform$OS.type == "windows") windowsFonts(`Times New Roman` = windowsFont("Times New Roman"))

save_candidate <- function(draw, stem, height_mm, pointsize = 12) {
  base <- file.path(out, stem)
  width_mm <- 183
  png(paste0(base, ".png"), width = width_mm, height = height_mm, units = "mm",
      res = 300, pointsize = pointsize, family = font, bg = "white")
  draw(); dev.off()
  tiff(paste0(base, ".tiff"), width = width_mm, height = height_mm, units = "mm",
       res = 600, compression = "lzw", pointsize = pointsize, family = font, bg = "white")
  draw(); dev.off()
  cairo_pdf(paste0(base, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
            pointsize = pointsize, family = font)
  draw(); dev.off()
  svg(paste0(base, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4,
      pointsize = pointsize, family = font)
  draw(); dev.off()
}

# Figure 1: preserve the frozen information architecture, but use the project typeface.
old_out <- Sys.getenv("HCC_OUTPUT_DIR", unset = NA_character_)
Sys.setenv(HCC_OUTPUT_DIR = scratch)
source(file.path(root, "05_scripts", "65_redesign_figure1_information_architecture.R"), local = FALSE)
font <- "Times New Roman"
save_candidate(draw, "Figure_1_APPROVAL", 145, 11.5)

# Figures 2-5: preserve every value, while correcting the over-tall Q4 canvases.
source(file.path(root, "05_scripts", "68_redesign_figures2_to6_q4.R"), local = FALSE)
font <- "Times New Roman"

# Restore visible, aligned panel titles. The legend carries the detailed wording;
# titles here are intentionally concise, following the two reference papers.
pt <- function(letter, title) {
  compact <- c(
    "Genus effects across cohorts"="Cohort effects",
    "Small-k pooled estimates"="Small-k inference",
    "Directional concordance"="Concordance",
    "Promotion sequence"="Promotion gate",
    "Plasma promotions across four models"="Plasma candidates",
    "Paired tissue promotions and identity tier"="Liver tissue",
    "Identification level"="Identity tier",
    "Exact cross-specimen audit"="Cross-specimen gate",
    "Project evidence-state matrix"="Evidence states",
    "Prior-edge provenance by species"="Prior provenance",
    "Weight sensitivity cannot pass the gate"="Weight sensitivity",
    "Provenance gate"="Provenance gate",
    "Cross-dataset host-program direction and effect"="Host concordance",
    "50 patient-matched TCGA pairs"="Matched pairs",
    "Direct external validation"="External validation",
    "Contradiction is retained"="Retained contradiction"
  )
  if (title %in% names(compact)) title <- unname(compact[title])
  mtext(letter, 3, 1.52, adj = 0, font = 2, cex = .92, col = C["ink"], family = font)
  mtext(title, 3, 1.52, adj = .15, font = 2, cex = .72, col = C["ink"], family = font)
}

bump_top_margin <- function(fun) {
  txt <- paste(deparse(body(fun), width.cutoff = 500), collapse = "\n")
  txt <- gsub(", 1.4,", ", 2.35,", txt, fixed = TRUE)
  body(fun) <- parse(text = txt)[[1]]
  fun
}
draw2 <- bump_top_margin(draw2)
draw3 <- bump_top_margin(draw3)
draw4 <- bump_top_margin(draw4)
draw5 <- bump_top_margin(draw5)

compact_decision_text <- function(fun) {
  txt <- paste(deparse(body(fun), width.cutoff = 500), collapse = "\n")
  txt <- gsub("3 programmes\\nthree-dataset concordant", "3 programmes\\nconcordant in 3 datasets", txt, fixed = TRUE)
  txt <- gsub("Respiratory chain\\ndirectional contradiction", "Respiratory chain\\ncontradictory direction", txt, fixed = TRUE)
  txt <- gsub("cex = 0.67", "cex = 0.52", txt, fixed = TRUE)
  txt <- gsub("cex = 0.64", "cex = 0.50", txt, fixed = TRUE)
  body(fun) <- parse(text = txt)[[1]]
  fun
}
draw5 <- compact_decision_text(draw5)

# Short display labels keep Figure 5's matched-pair x axis inside its panel.
short_program <- c(
  cell_division = "Cell",
  respiratory_electron_transport_chain = "Resp.",
  ribosome_biogenesis = "Ribo.",
  type_I_interferon_response = "IFN-I"
)
tc$program <- unname(short_program[as.character(tc$program)])
tc$program <- factor(tc$program, levels = unname(short_program))
rownames(m5) <- unname(short_program[rownames(m5)])
gse$program <- ifelse(as.character(gse$program) %in% names(short_program),
                      unname(short_program[as.character(gse$program)]), as.character(gse$program))

# A 2 x 2 reference-paper grid eliminates the unused lower-half space in Q4.
body(draw2)[[3]] <- quote(layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.35, 1), heights = c(1.18, .82)))
body(draw3)[[3]] <- quote(layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.25, 1), heights = c(1.05, .95)))
body(draw4)[[3]] <- quote(layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.30, 1), heights = c(1.05, .95)))
body(draw5)[[3]] <- quote(layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.30, 1), heights = c(1.05, .95)))

save_candidate(draw2, "Figure_2_APPROVAL", 135, 12.5)
save_candidate(draw3, "Figure_3_APPROVAL", 140, 12.5)
save_candidate(draw4, "Figure_4_APPROVAL", 140, 12.5)
save_candidate(draw5, "Figure_5_APPROVAL", 140, 12.5)

# Current Figure 6 and Figure 7 are already the mature visual anchors. Their
# authoritative vector/raster exports are copied into the approval set without
# changing analytical content.
copy_set <- function(src_base, stem) {
  for (ext in c("png", "tiff", "pdf", "svg")) {
    src <- paste0(src_base, ".", ext)
    if (file.exists(src)) file.copy(src, file.path(out, paste0(stem, ".", ext)), overwrite = TRUE)
  }
}
copy_set(file.path(root, "07_figures", "single_cell_atlas_redesign", "Figure6_single_cell_atlas"), "Figure_6_APPROVAL")
copy_set(file.path(root, "07_figures", "second_scenario_crc_v2", "Figure7_CRC_leakage_free_framework_benchmark"), "Figure_7_APPROVAL")

# Supplementary Figures 1-2 use the approved large-type implementation.
Sys.setenv(HCC_OUTPUT_DIR = out)
source(file.path(root, "05_scripts", "63_redraw_supplementary_figures_q3_large_type.R"), local = FALSE)
copy_set(file.path(out, "SupplementaryFigure1_microbiome_sensitivity_Q3"), "Supplementary_Figure_1_APPROVAL")
copy_set(file.path(out, "SupplementaryFigure2_screening_corrections_Q3"), "Supplementary_Figure_2_APPROVAL")
copy_set(file.path(root, "07_figures", "single_cell_atlas_redesign", "Supplementary_Figure3_single_cell_atlas_QC"), "Supplementary_Figure_3_APPROVAL")

# Remove only intermediate Q3 duplicate exports from the approval directory.
unlink(file.path(out, "SupplementaryFigure1_microbiome_sensitivity_Q3.*"))
unlink(file.path(out, "SupplementaryFigure2_screening_corrections_Q3.*"))
if (!is.na(old_out)) Sys.setenv(HCC_OUTPUT_DIR = old_out) else Sys.unsetenv("HCC_OUTPUT_DIR")

cat(normalizePath(out, winslash = "/"), "\n")
