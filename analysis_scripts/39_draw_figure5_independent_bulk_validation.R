options(stringsAsFactors = FALSE)
root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
outdir <- Sys.getenv("HCC_OUTPUT_DIR", unset = "")
if (!nzchar(outdir)) outdir <- file.path(root, "07_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "05_scripts", "figure_theme_times_q2.R"))

cross <- read.csv(file.path(root, "06_results", "integration", "figure5_host_program_cross_dataset_source.csv"), check.names = FALSE)
ledger <- read.csv(file.path(root, "06_results", "integration", "figure5_host_program_contradiction_ledger.csv"), check.names = FALSE)
tcga <- read.csv(file.path(root, "06_results", "tcga_lihc_host_program_scores.csv"), check.names = FALSE)
gse <- read.csv(file.path(root, "06_results", "GSE63898", "host_programs", "gse63898_program_validation.csv"), check.names = FALSE)

programs <- c("cell_division", "ribosome_biogenesis", "type_I_interferon_response", "respiratory_electron_transport_chain")
plabs <- c("Cell division", "Ribosome biogenesis", "Type I IFN response", "Respiratory electron\ntransport")
names(plabs) <- programs
dsets <- c("PRJEB54571", "TCGA-LIHC", "GSE63898")
dcols <- c("PRJEB54571" = Q2_COL[["blue"]], "TCGA-LIHC" = "#669B97", "GSE63898" = Q2_COL[["teal"]])
up_col <- Q2_COL[["blue"]]; down_col <- Q2_COL[["amber"]]; conflict_col <- Q2_COL[["red"]]; neutral <- "#F2F2F2"

# Patient-level TCGA paired differences, retaining only patients with exactly one tumour and one normal score per programme.
paired_rows <- list()
for (p in programs) {
  z <- tcga[tcga$module == p, c("patient_id", "group", "score")]
  wide <- reshape(z, idvar = "patient_id", timevar = "group", direction = "wide")
  wide <- wide[is.finite(wide$score.Tumor) & is.finite(wide$score.Normal), ]
  paired_rows[[p]] <- data.frame(patient_id = wide$patient_id, program = p,
                                 paired_difference = wide$score.Tumor - wide$score.Normal)
}
paired <- do.call(rbind, paired_rows); rownames(paired) <- NULL
write.csv(paired, file.path(outdir, "Figure5_source_data_TCGA_paired_differences.csv"), row.names = FALSE)

panel_label <- function(lab) mtext(lab, side = 3, line = 0.16, at = par("usr")[1],
                                   adj = 0, font = 2, cex = 0.88, family = Q2_FONT)

draw_figure <- function() {
  layout(matrix(c(1, 1, 2, 3, 4, 5), nrow = 3, byrow = TRUE), heights = c(0.55, 1.05, 1.05), widths = c(1.1, 1))
  par(family = Q2_FONT, fg = Q2_COL[["ink"]], col.axis = Q2_COL[["ink"]], col.lab = Q2_COL[["ink"]])

  # a: study and contrast registry
  par(mar = c(1.2, 1.1, 2.0, 1.0))
  plot(NA, xlim = c(0.5, 3.5), ylim = c(0.05, 1.0), axes = FALSE, xlab = "", ylab = "")
  card_titles <- c("Discovery", "Paired external validation", "Direct external validation")
  card_body <- c("PRJEB54571 RNA-seq\n17 HCC vs 5 cirrhosis\n22 participants\nGSEA NES",
                 "TCGA-LIHC\nTumour vs matched adjacent\n50 patient pairs\nMedian paired difference",
                 "GSE63898\n228 HCC vs 168 cirrhosis\n396 participants\nPlate-adjusted beta")
  for (i in 1:3) {
    rect(i - 0.43, 0.14, i + 0.43, 0.90, col = adjustcolor(dcols[dsets[i]], 0.10), border = dcols[dsets[i]], lwd = 1.3)
    text(i, 0.78, card_titles[i], font = 2, cex = 0.72, col = dcols[dsets[i]])
    text(i, 0.48, card_body[i], cex = 0.65)
    if (i < 3) arrows(i + 0.45, 0.52, i + 0.53, 0.52, length = 0.08, lwd = 1.2, col = "#888888")
  }
  q2_panel("A", "Three independent human bulk-liver contrasts", title_adj = .105)

  # b: dataset-specific effect matrix
  par(mar = c(5.8, 8.0, 2.0, 1.0))
  plot(NA, xlim = c(0.5, 3.5), ylim = c(0.5, 4.5), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = 1:3, labels = c("PRJEB\nNES", "TCGA\npaired median", "GSE63898\nadjusted beta"), tick = FALSE, cex.axis = 0.68)
  axis(2, at = 4:1, labels = plabs[programs], las = 2, tick = FALSE, cex.axis = 0.68)
  for (i in seq_along(programs)) {
    for (j in seq_along(dsets)) {
      z <- cross[cross$program == programs[i] & cross$dataset == dsets[j], ]
      yy <- 5 - i
      cc <- if (z$effect > 0) up_col else down_col
      rect(j - 0.48, yy - 0.46, j + 0.48, yy + 0.46, col = adjustcolor(cc, 0.25), border = "white")
      arrows(j, yy - 0.08, j, yy + 0.20 * sign(z$effect), length = 0.10, lwd = 2, col = cc)
      text(j, yy - 0.24, sprintf("%.3g", z$effect), cex = 0.63, font = 2)
    }
  }
  box(col = "#BBBBBB")
  q2_panel("B", "Effects remain on dataset-specific scales", title_adj = .105)
  mtext("Blue = higher in HCC; amber = lower in HCC", side = 1, line = 4.6, cex = 0.58, col = "#666666")

  # c: patient-level TCGA paired differences
  par(mar = c(6.3, 4.5, 2.0, 1.0))
  ylim <- range(paired$paired_difference, finite = TRUE)
  plot(NA, xlim = c(0.5, 4.5), ylim = ylim, axes = FALSE, xlab = "", ylab = "Tumour - matched adjacent score")
  axis(1, at = 1:4, labels = c("Cell\ndivision", "Ribosome\nbiogenesis", "Type I\nIFN", "Respiratory\nchain"), tick = FALSE, cex.axis = 0.68)
  axis(2, las = 1, cex.axis = 0.65)
  abline(h = 0, col = "#999999", lty = 2)
  set.seed(20260817)
  for (i in seq_along(programs)) {
    v <- paired$paired_difference[paired$program == programs[i]]
    points(jitter(rep(i, length(v)), amount = 0.11), v, pch = 16, cex = 0.44,
           col = adjustcolor(if (i == 4) conflict_col else dcols["TCGA-LIHC"], 0.48))
    segments(i - 0.22, median(v), i + 0.22, median(v), lwd = 2.4,
             col = if (i == 4) conflict_col else dcols["TCGA-LIHC"])
  }
  q2_panel("C", "Paired TCGA-LIHC validation (n=50 pairs)", title_adj = .105)
  mtext("Dots are patient-level paired differences; bars are medians", side = 1, line = 5.0, cex = 0.58, col = "#666666")

  # d: GSE63898 direct HCC-versus-cirrhosis forest
  par(mar = c(4.2, 8.0, 2.0, 1.0))
  gp <- gse[gse$program %in% programs, ]
  gp <- gp[match(programs, gp$program), ]
  xlim <- range(c(gp$plate_adjusted_ci_low, gp$plate_adjusted_ci_high, 0), finite = TRUE)
  plot(NA, xlim = xlim, ylim = c(0.5, 4.5), axes = FALSE, xlab = "Plate-adjusted programme-score beta", ylab = "")
  axis(1, cex.axis = 0.65); axis(2, at = 4:1, labels = plabs[programs], las = 2, tick = FALSE, cex.axis = 0.68)
  abline(v = 0, col = "#999999", lty = 2)
  for (i in 1:4) {
    yy <- 5 - i; cc <- if (i == 4) conflict_col else dcols["GSE63898"]
    segments(gp$plate_adjusted_ci_low[i], yy, gp$plate_adjusted_ci_high[i], yy, lwd = 2, col = cc)
    points(gp$plate_adjusted_beta[i], yy, pch = 21, bg = cc, col = "white", cex = 1.2)
    text(gp$plate_adjusted_ci_high[i], yy + 0.22, paste0("q=", format(gp$plate_adjusted_q_within_family[i], digits = 2, scientific = TRUE)), adj = 1, cex = 0.55, col = "#555555")
  }
  q2_panel("D", "Direct HCC-versus-cirrhosis validation", title_adj = .105)
  mtext("228 HCC and 168 cirrhotic tissues; 95% CI", side = 3, line = -0.75,
        adj = 1, cex = 0.58, col = "#666666")

  # e: replication decision / contradiction ledger
  par(mar = c(4.3, 6.0, 2.0, 1.0))
  plot(NA, xlim = c(0.5, 3.5), ylim = c(0.1, 4.6), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = 1:3, labels = c("Discovery", "Paired\nexternal", "Direct\nexternal"), tick = FALSE, cex.axis = 0.67)
  axis(2, at = 4:1, labels = c("Cell division", "Ribosome", "Type I IFN", "Respiratory"), las = 2, tick = FALSE, cex.axis = 0.68)
  for (i in seq_along(programs)) for (j in seq_along(dsets)) {
    z <- cross[cross$program == programs[i] & cross$dataset == dsets[j], ]
    yy <- 5 - i; cc <- if (z$effect > 0) up_col else down_col
    points(j, yy, pch = if (z$effect > 0) 24 else 25, bg = cc, col = cc, cex = 1.25)
  }
  rect(0.65, 0.14, 3.35, 0.66, col = "#FCE9E9", border = conflict_col, lwd = 1.1)
  text(2, 0.47, "3/4 strict replications", font = 2, cex = 0.72, col = "#2E6F5B")
  text(2, 0.25, "Respiratory programme: significant direction conflict", cex = 0.56, col = "#9E2F2F")
  q2_panel("E", "Direction-aware replication decision", title_adj = .105)

}

base <- file.path(outdir, "Figure5_independent_bulk_host_validation")
q2_save_base(draw_figure, base, width_mm = 180, height_mm = 135)

cat("Figure 5 exports written to", outdir, "\n")
cat("TCGA paired rows:", nrow(paired), "\n")
