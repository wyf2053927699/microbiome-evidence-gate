options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
outdir <- Sys.getenv("HCC_OUTPUT_DIR", unset = "")
if (!nzchar(outdir)) outdir <- file.path(root, "07_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(file.path(root, "06_results", "axis_prioritisation", "axis_dimension_scores.csv"), check.names = FALSE)
sens <- read.csv(file.path(root, "06_results", "axis_prioritisation", "figure4_weight_sensitivity.csv"), check.names = FALSE)
edges <- read.csv(file.path(root, "06_results", "axis_prioritisation", "edge_level_provenance.csv"), check.names = FALSE)

dims <- c("microbial_reproducibility", "human_metabolite_observed",
          "microbe_metabolite_link", "metabolite_target_link",
          "host_reproduced_two_datasets", "single_cell_localisation",
          "major_contradiction_resolved", "robust_to_analytic_choices")
dim_labs <- c("Microbial\nreproduction", "Human\nmetabolite", "Microbe-\nmetabolite", "Metabolite-\ntarget",
              "Bulk host\nreproduction", "Single-cell\nlocalisation", "Contradiction\nresolved", "Analytic\nrobustness")
short <- c("Bile acid", "Aromatic/AHR", "Indole/AHR", "Choline/GPC")
names(short) <- scores$candidate_family
cols <- c("Bile acid" = "#4477AA", "Aromatic/AHR" = "#CC6677", "Indole/AHR" = "#AA4499", "Choline/GPC" = "#228833")

panel_label <- function(x, y, lab) mtext(lab, side = 3, line = 0.65, at = 0, adj = 0,
                                         font = 2, cex = 1.05)

draw_figure <- function() {
  layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE), widths = c(1.22, 1), heights = c(1.06, 1))
  par(family = "Arial", fg = "#222222", col.axis = "#222222", col.lab = "#222222")

  # a: measured-evidence heat map
  par(mar = c(7.2, 8.4, 2.1, 1.0))
  mat <- as.matrix(scores[, dims])
  mat <- mat[nrow(mat):1, , drop = FALSE]
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat), zlim = c(0, 1),
        col = colorRampPalette(c("#F2F2F2", "#F3D89A", "#3A7D6B"))(101),
        axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq_along(dim_labs), labels = dim_labs, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.5)
  axis(2, at = seq_len(nrow(mat)), labels = rev(short[scores$candidate_family]), las = 2, tick = FALSE, cex.axis = 0.72)
  for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
    text(j, i, sprintf("%.2g", mat[i, j]), cex = 0.62,
         col = if (mat[i, j] >= 0.7) "white" else "#333333", font = 2)
  }
  box(col = "#BBBBBB")
  title("Project-measured evidence gates", adj = 0.06, cex.main = 0.92, font.main = 2, line = 0.6)
  mtext("0 = absent/fail; 0.25-0.75 = partial; 1 = complete", side = 3, adj = 1, cex = 0.58, line = 0.7, col = "#666666")
  panel_label(par("usr")[1], par("usr")[4] + 0.38, "A")

  # b: provenance-coded external edge inventory
  par(mar = c(4.5, 6.8, 2.1, 1.0))
  fam <- unique(scores$candidate_family)
  edge_counts <- table(factor(edges$candidate_family, levels = fam))
  y <- rev(seq_along(fam))
  plot(NA, xlim = c(0.7, 3.35), ylim = c(0.45, 4.55), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = 1:3, labels = c("Microbial\nfunction", "Metabolite", "Host\nprogramme"), tick = FALSE, cex.axis = 0.7)
  axis(2, at = y, labels = short[fam], las = 2, tick = FALSE, cex.axis = 0.7)
  abline(v = 1:3, col = "#E5E5E5", lty = 3)
  for (i in seq_along(fam)) {
    yy <- y[i]
    segments(1, yy, 3, yy, lty = 2, lwd = 1.4, col = adjustcolor(cols[short[fam[i]]], 0.7))
    points(1:3, rep(yy, 3), pch = 21, bg = "white", col = cols[short[fam[i]]], lwd = 1.3, cex = 1.15)
    text(2, yy + 0.22, paste0("external edges n=", edge_counts[i]), cex = 0.62, col = "#444444")
    text(2, yy - 0.22, "project-eligible edges n=0", cex = 0.60, col = "#B23A3A", font = 2)
  }
  legend("bottom", inset = -0.26, xpd = NA, horiz = TRUE, bty = "n", cex = 0.62,
         lty = c(2, 1), lwd = c(1.4, 2), col = c("#777777", "#B23A3A"),
         legend = c("external plausibility", "project-validated edge (none)"))
  title("Edge provenance does not close project links", adj = 0.06, cex.main = 0.92, font.main = 2, line = 0.6)
  panel_label(par("usr")[1], par("usr")[4] + 0.35, "B")

  # c: weighting sensitivity
  par(mar = c(6.1, 4.4, 2.1, 1.0))
  schemes <- c("equal", "measured_layer_emphasis", "bridge_emphasis", "validation_emphasis")
  scheme_labs <- c("Equal", "Measured\nlayer", "Bridge", "Validation")
  plot(NA, xlim = c(0.85, 4.15), ylim = c(0, 0.48), axes = FALSE, xlab = "", ylab = "Weighted score")
  axis(1, at = 1:4, labels = scheme_labs, tick = FALSE, cex.axis = 0.7)
  axis(2, las = 1, cex.axis = 0.68)
  abline(h = seq(0, 0.5, 0.1), col = "#EEEEEE", lty = 1)
  for (f in unique(sens$candidate_family)) {
    z <- sens[sens$candidate_family == f, ]
    z <- z[match(schemes, z$weighting_scheme), ]
    lab <- short[f]
    lines(1:4, z$weighted_score, type = "o", pch = 16, lwd = 1.6, cex = 0.75, col = cols[lab])
  }
  legend("topright", legend = names(cols), col = cols, lty = 1, pch = 16, bty = "n", cex = 0.62)
  title("Rank stability across prespecified weights", adj = 0.06, cex.main = 0.92, font.main = 2, line = 0.6)
  mtext("Rank order is stable, but scores are not promotion probabilities", side = 1, line = 4.5, cex = 0.58, col = "#666666")
  panel_label(par("usr")[1], par("usr")[4] + 0.038, "C")

  # d: essential-gate decision matrix and bounded conclusion
  par(mar = c(3.0, 6.8, 2.1, 1.0))
  plot(NA, xlim = c(0.5, 8.5), ylim = c(0.1, 5.3), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = 1:8, labels = 1:8, tick = FALSE, cex.axis = 0.68)
  axis(2, at = 4:1, labels = short[scores$candidate_family], las = 2, tick = FALSE, cex.axis = 0.7)
  vals <- as.matrix(scores[, dims])
  for (i in 1:4) for (j in 1:8) {
    status_col <- if (vals[i, j] >= 1) "#3A7D6B" else if (vals[i, j] > 0) "#D7A84B" else "#C95A5A"
    points(j, 5 - i, pch = 21, bg = status_col, col = "white", lwd = 0.7, cex = 1.2)
  }
  rect(0.58, 0.15, 8.42, 0.82, col = "#FCE9E9", border = "#C95A5A", lwd = 1.2)
  text(4.5, 0.56, "NO VALIDATED MICROBIOTA-METABOLITE-HOST AXIS", col = "#9E2F2F", font = 2, cex = 0.82)
  text(4.5, 0.29, "All eight gates are required; weighting cannot rescue a failed gate", col = "#555555", cex = 0.60)
  legend("top", horiz = TRUE, bty = "n", inset = -0.03, xpd = NA, cex = 0.62,
         pt.bg = c("#3A7D6B", "#D7A84B", "#C95A5A"), pch = 21, pt.cex = 1.1,
         legend = c("complete", "partial", "failed"))
  title("Essential-gate decision", adj = 0.06, cex.main = 0.92, font.main = 2, line = 0.6)
  panel_label(par("usr")[1], par("usr")[4] + 0.42, "D")

}

base <- file.path(outdir, "Figure4_evidence_aware_axis_prioritisation")
width_mm <- 180
height_mm <- 137
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in,
                   system_fonts = list(Arial = "Arial"), pointsize = 8)
} else {
  svg(paste0(base, ".svg"), width = width_in, height = height_in, family = "Arial", pointsize = 8)
}
par(oma = c(0.2, 0.2, 0.4, 0.2)); draw_figure(); dev.off()

cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 8)
par(oma = c(0.2, 0.2, 0.4, 0.2)); draw_figure(); dev.off()

tiff(paste0(base, ".tiff"), width = width_mm, height = height_mm, units = "mm", res = 600,
     compression = "lzw", pointsize = 8, family = "Arial")
par(oma = c(0.2, 0.2, 0.4, 0.2)); draw_figure(); dev.off()

png(paste0(base, ".png"), width = width_mm, height = height_mm, units = "mm", res = 600,
    pointsize = 8, family = "Arial")
par(oma = c(0.2, 0.2, 0.4, 0.2)); draw_figure(); dev.off()

cat("Figure 4 exports written to", outdir, "\n")
