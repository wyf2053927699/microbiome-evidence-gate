options(stringsAsFactors = FALSE)

library(ggplot2)
library(patchwork)
library(ggrepel)
library(scales)
library(svglite)
library(ragg)

root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
atlas_dir <- file.path(root, "06_results", "single_cell", "atlas")
paired_dir <- file.path(root, "06_results", "single_cell", "paired_programs")
out_dir <- Sys.getenv("HCC_FIGURE_OUTPUT_DIR", unset = file.path(root, "07_figures", "single_cell_atlas_redesign"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Times New Roman"
cell_cols <- c(
  B = "#4477AA", Endothelial = "#66A182", Fibroblast = "#C58B45",
  Hepatocyte = "#D65F5F", Myeloid = "#7A6FAC", `T/NK` = "#2A9D8F"
)
site_cols <- c(Normal = "#4C78A8", Tumor = "#D55E5E", Lymph = "#999999", PVTT = "#8C6BB1")
ink <- "#18242B"; muted <- "#66757C"; rule <- "#D9E1E4"

theme_pub <- function(base_size = 8.4) {
  theme_classic(base_size = base_size, base_family = "Times New Roman") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = ink),
      axis.ticks = element_line(linewidth = 0.35, colour = ink),
      axis.text = element_text(colour = ink, size = 7.1),
      axis.title = element_text(colour = ink, face = "bold", size = 8.0),
      legend.title = element_text(face = "bold"),
      legend.key.height = unit(4.1, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = ink),
      plot.title = element_text(face = "bold", colour = ink),
      plot.tag = element_text(family = font_family, face = "bold", size = 10.8, colour = ink),
      plot.margin = margin(5, 6, 5, 6)
    )
}
theme_set(theme_pub())

save_pub <- function(plot, stem, width_mm = 180, height_mm = 225, dpi = 600) {
  base <- file.path(out_dir, stem)
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite(paste0(base, ".svg"), width = w, height = h, system_fonts = list(serif = font_family))
  print(plot); dev.off()
  cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = font_family)
  print(plot); dev.off()
  agg_tiff(paste0(base, ".tiff"), width = width_mm, height = height_mm,
           units = "mm", res = dpi, compression = "lzw", background = "white")
  print(plot); dev.off()
  agg_png(paste0(base, ".png"), width = width_mm, height = height_mm,
          units = "mm", res = 300, background = "white")
  print(plot); dev.off()
}

# Descriptive cell-atlas layer. No cell-level P values are computed.
pca <- read.csv(gzfile(file.path(atlas_dir, "gse149614_cell_pca_metadata.csv.gz")), check.names = FALSE)
pc_cols <- grep("^PC[0-9]+$", names(pca), value = TRUE)
umap_file <- file.path(atlas_dir, "gse149614_umap_metadata.csv.gz")
if (file.exists(umap_file)) {
  saved_umap <- read.csv(gzfile(umap_file), check.names = FALSE)
  if (nrow(saved_umap) != nrow(pca) || !identical(saved_umap$Cell, pca$Cell)) {
    stop("Saved UMAP does not match the PCA cell manifest")
  }
  pca$UMAP1 <- saved_umap$UMAP1
  pca$UMAP2 <- saved_umap$UMAP2
} else {
  if (!requireNamespace("uwot", quietly = TRUE)) stop("uwot is required only when the saved UMAP is absent")
  set.seed(149614)
  embedding <- uwot::umap(
    as.matrix(pca[, pc_cols]), n_neighbors = 30, min_dist = 0.30,
    metric = "cosine", n_components = 2, n_threads = 1,
    ret_model = FALSE, verbose = TRUE, init = "spectral"
  )
  pca$UMAP1 <- embedding[, 1]
  pca$UMAP2 <- embedding[, 2]
  write.csv(pca[, c("Cell", "sample", "patient", "site", "celltype", "UMAP1", "UMAP2")],
            umap_file, row.names = FALSE)
}

centres <- aggregate(cbind(UMAP1, UMAP2) ~ celltype, pca, median)
p_a <- ggplot(pca, aes(UMAP1, UMAP2, colour = celltype)) +
  geom_point(size = 0.17, alpha = 0.52, stroke = 0) +
  geom_label_repel(data = centres, aes(label = celltype), colour = ink,
                   family = font_family, size = 3.1, fontface = "bold",
                   label.size = 0.2, label.padding = unit(1.3, "mm"),
                   box.padding = 0.35, point.padding = 0.25, seed = 149614,
                   min.segment.length = 0, show.legend = FALSE) +
  scale_colour_manual(values = cell_cols) +
  coord_equal() + labs(x = "UMAP 1", y = "UMAP 2", colour = "Cell class") +
  theme(legend.position = "none")

p_b <- ggplot(pca, aes(UMAP1, UMAP2, colour = site)) +
  geom_point(size = 0.15, alpha = 0.38, stroke = 0) +
  scale_colour_manual(values = site_cols) + coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(size = 7.4), legend.key.width = unit(4.2, "mm"))

# Annotation-validation layer.
markers <- read.csv(file.path(atlas_dir, "gse149614_marker_dotplot.csv"))
marker_order <- unique(markers$gene)
markers$gene <- factor(markers$gene, levels = rev(marker_order))
markers$observed_celltype <- factor(markers$observed_celltype, levels = names(cell_cols))
markers$scaled_mean <- ave(markers$mean_log1p_cpm, markers$gene,
                           FUN = function(x) if (sd(x) == 0) 0 else as.numeric(scale(x)))
p_c <- ggplot(markers, aes(observed_celltype, gene)) +
  geom_point(aes(size = pct_detected, colour = scaled_mean)) +
  scale_size_continuous(name = "% detected", range = c(0.8, 6.2), breaks = c(25, 50, 75)) +
  scale_colour_gradient2(name = "Scaled mean", low = "#3B6FB6", mid = "#F7F7F7",
                         high = "#C94C4C", midpoint = 0, limits = c(-2.2, 2.2), oob = squish) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(title.position = "top", nrow = 1),
         colour = guide_colourbar(title.position = "top", barwidth = unit(24, "mm"),
                                  barheight = unit(2.8, "mm"))) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 6.3),
        axis.text.y = element_text(size = 6.3),
        panel.grid = element_line(colour = "#EDF1F2", linewidth = 0.3),
        legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = 7.4), legend.text = element_text(size = 7.0),
        legend.margin = margin(0, 0, 0, 0))

# Donor-level composition; individual patients, not cells, are the units.
comp <- read.csv(file.path(atlas_dir, "gse149614_donor_celltype_composition.csv"))
comp <- comp[comp$site %in% c("Normal", "Tumor"), ]
wide <- reshape(comp[, c("patient", "site", "celltype", "proportion")],
                idvar = c("patient", "celltype"), timevar = "site", direction = "wide")
composition_rows_before_pairing <- nrow(wide)
wide <- wide[!is.na(wide$proportion.Normal) & !is.na(wide$proportion.Tumor), ]
composition_rows_after_pairing <- nrow(wide)
wide$delta <- 100 * (wide$proportion.Tumor - wide$proportion.Normal)
sum_comp <- do.call(rbind, lapply(split(wide, wide$celltype), function(z) {
  data.frame(celltype = z$celltype[1], mean = mean(z$delta),
             low = mean(z$delta) - qt(.975, df = nrow(z) - 1) * sd(z$delta) / sqrt(nrow(z)),
             high = mean(z$delta) + qt(.975, df = nrow(z) - 1) * sd(z$delta) / sqrt(nrow(z)),
             n = nrow(z))
}))
sum_comp$celltype <- factor(sum_comp$celltype, levels = names(cell_cols))
wide$celltype <- factor(wide$celltype, levels = names(cell_cols))
p_d <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, colour = rule) +
  geom_jitter(data = wide, aes(celltype, delta), width = 0.10, height = 0,
              size = 1.6, shape = 21, fill = "white", colour = muted, stroke = 0.4) +
  geom_errorbar(data = sum_comp, aes(celltype, ymin = low, ymax = high, colour = celltype),
                width = 0, linewidth = 0.8) +
  geom_point(data = sum_comp, aes(celltype, mean, colour = celltype), size = 3.0) +
  scale_colour_manual(values = cell_cols) +
  labs(x = NULL, y = "Change in composition\n(percentage points)") +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 6.5), axis.text.y = element_text(size = 6.5), legend.position = "none",
        axis.title.y = element_text(size = 8.0), plot.margin = margin(5, 6, 5, 11))

# All estimable prespecified programme-by-cell-class tests are retained.
exact <- read.csv(file.path(paired_dir, "prespecified_program_exact_signflip.csv"), check.names = FALSE)
program_labels <- c(
  ahr_aromatic_response = "AHR/aromatic response",
  cell_cycle = "Cell cycle",
  choline_phospholipid_metabolism = "Choline/phospholipid",
  hepatocyte_bile_acid_receptor_response = "Bile-acid receptor",
  hepatocyte_bile_acid_synthesis_transport = "Bile-acid synthesis/transport",
  type_i_interferon_response = "Type-I interferon"
)
program_short <- c(
  ahr_aromatic_response = "AHR",
  cell_cycle = "Cell cycle",
  choline_phospholipid_metabolism = "Choline",
  hepatocyte_bile_acid_receptor_response = "BA receptor",
  hepatocyte_bile_acid_synthesis_transport = "BA synthesis",
  type_i_interferon_response = "IFN-I"
)
exact$program_label <- unname(program_labels[exact$program])
exact$celltype <- factor(exact$celltype, levels = names(cell_cols))
exact$program_label <- factor(exact$program_label, levels = unname(program_labels))
p_e <- ggplot(exact, aes(program_label, celltype, fill = mean_tumor_minus_normal)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", mean_tumor_minus_normal)), family = font_family,
            size = 3.0, colour = ink) +
  scale_fill_gradient2(name = "Mean paired\ndifference", low = "#3B6FB6", mid = "white",
                       high = "#C94C4C", midpoint = 0) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1, size = 6.1), axis.text.y = element_text(size = 6.1), panel.grid = element_blank(),
        legend.position = "right")

exact$test_label <- paste(exact$celltype, exact$program_label, sep = " | ")
exact <- exact[order(exact$mean_tumor_minus_normal), ]
exact$test_label <- factor(exact$test_label, levels = exact$test_label)
p_f <- ggplot(exact, aes(mean_tumor_minus_normal, test_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = rule) +
  geom_errorbar(aes(xmin = mean_difference_ci95_low, xmax = mean_difference_ci95_high,
                    colour = celltype), orientation = "y", width = 0, linewidth = 0.65) +
  geom_point(aes(colour = celltype), size = 2.3) +
  scale_colour_manual(values = cell_cols) +
  labs(x = "Mean paired difference (95% CI)", y = NULL,
       subtitle = sprintf("Exact sign-flip: %d tests; min BH q = %.3f",
                          nrow(exact), min(exact$exact_signflip_q_global_bh))) +
  theme(axis.text.y = element_text(size = 6.0), legend.position = "none",
        plot.subtitle = element_text(size = 7.6, colour = muted, margin = margin(b = 2)),
        axis.title.x = element_text(size = 8.2))

main_fig <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_layout(widths = c(1, 1), heights = c(1.25, 1.22, 1.38), guides = "keep") +
  plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0.01, 0.99))
save_pub(main_fig, "Figure6_single_cell_atlas", width_mm = 210, height_mm = 255)

# Supplementary atlas diagnostics and robustness.
qc <- jsonlite::fromJSON(file.path(atlas_dir, "gse149614_atlas_qc.json"))
scree <- data.frame(PC = seq_along(qc$pca_variance_explained),
                    variance = 100 * unlist(qc$pca_variance_explained))
s_a <- ggplot(scree, aes(PC, variance)) + geom_col(fill = "#78909C", width = 0.78) +
  labs(x = "Principal component", y = "Variance explained (%)")

s_b <- ggplot(pca, aes(UMAP1, UMAP2, colour = patient)) +
  geom_point(size = 0.13, alpha = 0.35, stroke = 0) + coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Donor") +
  theme(legend.position = "right")

scores <- read.csv(file.path(paired_dir, "prespecified_program_scores.csv"))
best <- exact[which.min(exact$exact_signflip_q_global_bh), ]
sel <- scores[scores$celltype == as.character(best$celltype) & scores$program == best$program &
                scores$site %in% c("Normal", "Tumor"), ]
s_c <- ggplot(sel, aes(site, module_score, group = patient)) +
  geom_line(colour = "#87959B", linewidth = 0.55) +
  geom_point(aes(fill = site), shape = 21, size = 2.5, colour = "white", stroke = 0.35) +
  scale_fill_manual(values = site_cols) +
  labs(x = NULL, y = "Programme score",
       title = paste(as.character(best$celltype), program_short[best$program], sep = " | ")) +
  theme(legend.position = "none", plot.title = element_text(size = 7.8))

loo <- read.csv(file.path(paired_dir, "prespecified_program_leave_one_donor_out.csv"))
loo$label <- paste(loo$celltype, unname(program_short[loo$program]), sep = " | ")
s_d <- ggplot(loo, aes(omitted_patient, label, fill = mean_difference)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C94C4C", midpoint = 0) +
  labs(x = "Omitted donor", y = NULL, fill = "Mean paired\ndifference") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 6.8))

exact$test_label_short <- paste(exact$celltype, unname(program_short[exact$program]), sep = " | ")
exact$test_label_short <- factor(exact$test_label_short, levels = exact$test_label_short)
s_e <- ggplot(exact, aes(exact_signflip_q_global_bh, test_label_short, colour = celltype)) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = "#B44D4D") +
  geom_point(size = 2.2) + scale_colour_manual(values = cell_cols) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1)) +
  labs(x = "Global exact sign-flip BH q", y = NULL) +
  theme(axis.text.y = element_text(size = 6.8), legend.position = "none")

counts <- aggregate(Cell ~ patient + site, pca, length)
names(counts)[3] <- "cells"
s_f <- ggplot(counts, aes(patient, cells, fill = site)) +
  geom_col(position = "dodge", width = 0.75) + scale_fill_manual(values = site_cols) +
  scale_y_continuous(labels = comma) + labs(x = "Donor", y = "Retained cells", fill = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

supp_fig <- (s_a | s_b) / (s_c | s_d) / (s_e | s_f) +
  plot_layout(heights = c(0.85, 1.15, 1.10), widths = c(0.78, 1.22)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0.01, 0.99))
save_pub(supp_fig, "Supplementary_Figure3_single_cell_atlas_QC", height_mm = 225)

write.csv(exact, file.path(out_dir, "Figure6_source_data_exact_signflip_all_tests.csv"), row.names = FALSE)
write.csv(markers, file.path(out_dir, "Figure6_source_data_marker_dotplot.csv"), row.names = FALSE)
write.csv(wide, file.path(out_dir, "Figure6_source_data_donor_composition_differences.csv"), row.names = FALSE)
write.csv(data.frame(
  stage = c("before_complete_pair_requirement", "after_complete_pair_requirement"),
  rows = c(composition_rows_before_pairing, composition_rows_after_pairing),
  rule = c("donor-cell-class records after reshape", "both Normal and Tumor proportions observed")
), file.path(out_dir, "Figure6_composition_pairing_accounting.csv"), row.names = FALSE)
write.csv(centres, file.path(out_dir, "Figure6_source_data_umap_label_centres.csv"), row.names = FALSE)
