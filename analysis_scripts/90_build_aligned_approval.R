options(stringsAsFactors = FALSE)

root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out <- file.path(root, "07_figures", "layout_alignment_approval_2026-08-20")
scratch <- file.path(out, "_scratch")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)

PUBLICATION_FONT_FAMILY <- "Times New Roman"
FINAL_WIDTH_MM <- 183
MIN_FONT_SIZE_PT <- 8
font <- PUBLICATION_FONT_FAMILY
if (.Platform$OS.type == "windows") windowsFonts(`Times New Roman` = windowsFont("Times New Roman"))

save_candidate <- function(draw, stem, height_mm, pointsize = 13) {
  base <- file.path(out, stem); width_mm <- FINAL_WIDTH_MM
  png(paste0(base, ".png"), width = width_mm, height = height_mm, units = "mm",
      res = 300, pointsize = pointsize, family = "Times New Roman", bg = "white")
  draw(); dev.off()
  tiff(paste0(base, ".tiff"), width = width_mm, height = height_mm, units = "mm",
       res = 600, compression = "lzw", pointsize = pointsize, family = "Times New Roman", bg = "white")
  draw(); dev.off()
  cairo_pdf(paste0(base, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
            pointsize = pointsize, family = "Times New Roman")
  draw(); dev.off()
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(base, ".svg"), width = width_mm / 25.4,
                     height = height_mm / 25.4, pointsize = pointsize,
                     system_fonts = list(serif = PUBLICATION_FONT_FAMILY))
  } else {
    svg(paste0(base, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4,
        pointsize = pointsize, family = font)
  }
  draw(); dev.off()
}

copy_set <- function(src_base, stem) {
  for (ext in c("png", "tiff", "pdf", "svg")) {
    src <- paste0(src_base, ".", ext)
    if (!file.exists(src)) stop("Missing source export: ", src)
    file.copy(src, file.path(out, paste0(stem, ".", ext)), overwrite = TRUE)
  }
}

# Figure 1 — preserve the approved information architecture and enlarge the type.
old_out <- Sys.getenv("HCC_OUTPUT_DIR", unset = NA_character_)
Sys.setenv(HCC_OUTPUT_DIR = scratch)
source(file.path(root, "05_scripts", "65_redesign_figure1_information_architecture.R"), local = FALSE)
font <- "Times New Roman"
save_candidate(draw, "Figure_1_ALIGNED", 150, 13)

# Figures 2–5 — equal 2 x 2 panel cells, fixed tag anchor and generous plot margins.
source(file.path(root, "05_scripts", "68_redesign_figures2_to6_q4.R"), local = FALSE)
font <- "Times New Roman"
pt <- function(letter, title = NULL) {
  # Anchor tags to the panel cell, not to the inner plotting region.  This keeps
  # A/B and C/D on identical page-level x coordinates despite unequal y margins.
  x <- grconvertX(par("fig")[1] + .012, from = "ndc", to = "user")
  y <- grconvertY(par("fig")[4] - .014, from = "ndc", to = "user")
  text(x, y, letter, adj = c(0, 1), font = 2, cex = .96,
       col = C["ink"], family = font, xpd = NA)
}

normalise_panel_margins <- function(fun) {
  txt <- paste(deparse(body(fun), width.cutoff = 500), collapse = "\n")
  txt <- gsub(",1.4,", ",2.65,", txt, fixed = TRUE)
  txt <- gsub(", 1.4,", ", 2.65,", txt, fixed = TRUE)
  body(fun) <- parse(text = txt)[[1]]
  fun
}
draw2 <- normalise_panel_margins(draw2)
draw3 <- normalise_panel_margins(draw3)
draw4 <- normalise_panel_margins(draw4)
draw5 <- normalise_panel_margins(draw5)

# The decision cards in Figure 5 use short, two-line statements so the text
# remains inside the cards at final journal width.
txt5 <- paste(deparse(body(draw5), width.cutoff = 500), collapse = "\n")
txt5 <- gsub("3 programmes\\nthree-dataset concordant", "3 programmes\\nconcordant in 3 datasets", txt5, fixed = TRUE)
txt5 <- gsub("Respiratory chain\\ndirectional contradiction", "Respiratory chain\\ncontradictory direction", txt5, fixed = TRUE)
txt5 <- gsub("cex = 0.67", "cex = 0.48", txt5, fixed = TRUE)
txt5 <- gsub("cex = 0.64", "cex = 0.46", txt5, fixed = TRUE)
body(draw5) <- parse(text = txt5)[[1]]

# Short labels protect the plot region without changing analytical identities.
short_program <- c(cell_division = "Cell", respiratory_electron_transport_chain = "Resp.",
                   ribosome_biogenesis = "Ribo.", type_I_interferon_response = "IFN-I")
tc$program <- unname(short_program[as.character(tc$program)])
tc$program <- factor(tc$program, levels = unname(short_program))
rownames(m5) <- unname(short_program[rownames(m5)])
gse$program <- ifelse(as.character(gse$program) %in% names(short_program),
                      unname(short_program[as.character(gse$program)]), as.character(gse$program))

body(draw2)[[3]] <- quote(layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.08, .92), heights = c(1, 1)))
body(draw3)[[3]] <- quote(layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1, 1), heights = c(1, 1)))
body(draw4)[[3]] <- quote(layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.08, .92), heights = c(1, 1)))
body(draw5)[[3]] <- quote(layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.08, .92), heights = c(1, 1)))

save_candidate(draw2, "Figure_2_ALIGNED", 150, 13.5)
save_candidate(draw3, "Figure_3_ALIGNED", 150, 13.5)
save_candidate(draw4, "Figure_4_ALIGNED", 150, 13.5)
save_candidate(draw5, "Figure_5_ALIGNED", 150, 13.5)

# Figure 6 and Supplementary Figure 3 are the frozen atlas outputs. They already
# use equal two-column panel widths and a common Times New Roman theme.
copy_set(file.path(root, "07_figures", "single_cell_atlas_redesign", "Figure6_single_cell_atlas"), "Figure_6_ALIGNED")

# Figure 7 — re-execute the leakage-free benchmark on a taller canvas.
f7 <- readLines(file.path(root, "05_scripts", "87_draw_figure7_leakage_free_benchmark.R"), warn = FALSE)
f7 <- sub('outdir <- file.path\\(root, "07_figures", "second_scenario_crc_v2"\\)',
          'outdir <- out', f7)
f7 <- sub('width_mm <- 183; height_mm <- 170; pointsize <- 12',
          'width_mm <- 183; height_mm <- 185; pointsize <- 12.5', f7, fixed = TRUE)
# Detailed panel meanings remain in the legend. Removing long in-panel titles
# prevents the D/E headers from competing with their data at the lower edge.
f7 <- sub('  text\\(x\\+0.042, y-0.002, title,.*$', '  invisible(title)', f7)
f7 <- sub('text\\(\\.535,\\.190,', 'text(.515,.190,', f7)
f7 <- sub('text\\(\\.535,\\.158,', 'text(.515,.158,', f7)
f7 <- sub('text\\(\\.535,\\.118,', 'text(.515,.118,', f7)
f7 <- sub('text\\(\\.535,\\.087,', 'text(.515,.087,', f7)
f7 <- sub('cex=\\.60,col=pal\\["red"\\]', 'cex=.52,col=pal["red"]', f7)
eval(parse(text = f7), envir = .GlobalEnv)
copy_set(file.path(out, "Figure7_CRC_leakage_free_framework_benchmark"), "Figure_7_ALIGNED")

# Supplementary Figures 1–2 — retain large type but standardise the tag/title band.
sf <- readLines(file.path(root, "05_scripts", "63_redraw_supplementary_figures_q3_large_type.R"), warn = FALSE)
sf <- sub('pt<-function\\(l,t\\)\\{mtext\\(l,3,.18,adj=0,font=2,cex=1.45\\);mtext\\(t,3,.18,adj=.24,font=2,cex=1.10\\)\\}',
          'pt<-function(l,t){x<-grconvertX(par("fig")[1]+.012,"ndc","user");y<-grconvertY(par("fig")[4]-.010,"ndc","user");text(x,y,l,adj=c(0,1),font=2,cex=.92,xpd=NA)}', sf)
Sys.setenv(HCC_OUTPUT_DIR = out)
eval(parse(text = sf), envir = .GlobalEnv)
copy_set(file.path(out, "SupplementaryFigure1_microbiome_sensitivity_Q3"), "Supplementary_Figure_1_ALIGNED")
copy_set(file.path(out, "SupplementaryFigure2_screening_corrections_Q3"), "Supplementary_Figure_2_ALIGNED")
copy_set(file.path(root, "07_figures", "single_cell_atlas_redesign", "Supplementary_Figure3_single_cell_atlas_QC"), "Supplementary_Figure_3_ALIGNED")

unlink(file.path(out, "SupplementaryFigure1_microbiome_sensitivity_Q3.*"))
unlink(file.path(out, "SupplementaryFigure2_screening_corrections_Q3.*"))
unlink(file.path(out, "Figure7_CRC_leakage_free_framework_benchmark.*"))
if (!is.na(old_out)) Sys.setenv(HCC_OUTPUT_DIR = old_out) else Sys.unsetenv("HCC_OUTPUT_DIR")

cat(normalizePath(out, winslash = "/"), "\n")
