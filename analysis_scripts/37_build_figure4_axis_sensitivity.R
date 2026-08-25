options(stringsAsFactors = FALSE)

root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
}
infile <- file.path(root, "06_results", "axis_prioritisation", "axis_dimension_scores.csv")
outdir <- Sys.getenv("HCC_OUTPUT_DIR", unset = "")
if (!nzchar(outdir)) {
  outdir <- file.path(root, "06_results", "axis_prioritisation")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
x <- read.csv(infile, check.names = FALSE)

dims <- c(
  "microbial_reproducibility", "human_metabolite_observed",
  "microbe_metabolite_link", "metabolite_target_link",
  "host_reproduced_two_datasets", "single_cell_localisation",
  "major_contradiction_resolved", "robust_to_analytic_choices"
)

schemes <- list(
  equal = rep(1 / 8, 8),
  measured_layer_emphasis = c(0.18, 0.18, 0.12, 0.12, 0.14, 0.10, 0.08, 0.08),
  bridge_emphasis = c(0.10, 0.10, 0.22, 0.22, 0.12, 0.08, 0.08, 0.08),
  validation_emphasis = c(0.10, 0.10, 0.10, 0.10, 0.20, 0.16, 0.14, 0.10)
)

stopifnot(all(vapply(schemes, function(w) abs(sum(w) - 1) < 1e-9, logical(1))))
rows <- list()
for (nm in names(schemes)) {
  w <- schemes[[nm]]
  score <- as.numeric(as.matrix(x[, dims]) %*% w)
  ord <- order(-score, x$candidate_family)
  rank <- integer(length(score)); rank[ord] <- seq_along(ord)
  rows[[nm]] <- data.frame(
    weighting_scheme = nm,
    candidate_family = x$candidate_family,
    weighted_score = score,
    exploratory_rank = rank,
    essential_gate_pass = as.integer(x$essential_gate_pass),
    promotion_status = ifelse(x$essential_gate_pass == 1, "eligible", "NO-GO"),
    stringsAsFactors = FALSE
  )
}
res <- do.call(rbind, rows)
rownames(res) <- NULL
write.csv(res, file.path(outdir, "figure4_weight_sensitivity.csv"), row.names = FALSE, na = "")

rank_range <- aggregate(exploratory_rank ~ candidate_family, res,
                        function(z) paste0(min(z), "-", max(z)))
names(rank_range)[2] <- "rank_range_across_prespecified_schemes"
score_range <- aggregate(weighted_score ~ candidate_family, res,
                         function(z) paste0(format(min(z), digits = 4), "-", format(max(z), digits = 4)))
names(score_range)[2] <- "score_range_across_prespecified_schemes"
summary <- merge(rank_range, score_range, by = "candidate_family")
summary$essential_gate_pass_all_schemes <- 0L
summary$interpretation <- "Exploratory ranking only; failed essential gates prohibit axis promotion"
write.csv(summary, file.path(outdir, "figure4_rank_stability_summary.csv"), row.names = FALSE, na = "")

weights <- do.call(rbind, lapply(names(schemes), function(nm) {
  data.frame(weighting_scheme = nm, dimension = dims, weight = schemes[[nm]])
}))
write.csv(weights, file.path(outdir, "figure4_prespecified_weights.csv"), row.names = FALSE)

cat("Figure 4 sensitivity complete\n")
cat("Candidates:", nrow(x), "\n")
cat("Weighting schemes:", length(schemes), "\n")
cat("Eligible axes under essential gates:", sum(x$essential_gate_pass == 1), "\n")
