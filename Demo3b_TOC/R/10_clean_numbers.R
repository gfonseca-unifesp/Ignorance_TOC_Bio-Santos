# =============================================================================
# 10_clean_numbers.R — every number the report and manuscript quote, from the clean chain, in one file
# =============================================================================
# Writes outputs/comparison/clean_chain_numbers.txt (human-readable). Missing inputs are reported, not fatal.
Sys.setenv(DEMO3B_ITER = "iter2c_clean")
source("R/00_config.R")
I1 <- "iter1c_clean"; I2 <- "iter2c_clean"
out <- file.path(DIRS$compare, "clean_chain_numbers.txt")
sink(out, split = TRUE)
tb <- function(it, f) { p <- file.path(ROOT, "outputs", it, "tables", f); if (file.exists(p)) read.csv(p) else { cat("[missing]", it, f, "\n"); NULL } }
cm <- function(f) { p <- file.path(DIRS$compare, f); if (file.exists(p)) read.csv(p) else { cat("[missing]", f, "\n"); NULL } }
hdr <- function(s) cat("\n=====", s, "=====\n")
op <- options(width = 220, digits = 4)

hdr("data: North Sea rule and noise floor (01d)"); print(cm("northsea_rule_counts.csv"), row.names = FALSE)
for (it in c(I1, I2)) {
  hdr(paste("summary", it)); s <- tb(it, "toc_summary.csv"); if (!is.null(s)) print(s[, 1:2], row.names = FALSE)
  hdr(paste("calibration", it)); print(tb(it, "calibration_metrics_toc.csv"), row.names = FALSE)
  hdr(paste("calibration per region", it)); print(tb(it, "calibration_per_region_toc.csv"), row.names = FALSE)
  hdr(paste("variable importance", it)); v <- tb(it, "variable_importance_toc.csv"); if (!is.null(v)) print(v[order(-v$importance), ][1:8, ], row.names = FALSE)
  o <- tb(it, "calibration_outer_points_toc.csv"); if (!is.null(o)) cat("withheld RMSE", it, sqrt(mean((o$pred - o$obs)^2)), "\n")
}
hdr("diagnosis correlations iter1c"); print(tb(I1, "diagnosis_correlations_toc.csv"), row.names = FALSE)
hdr("diagnosis worst regions iter1c"); print(tb(I1, "diagnosis_worst_regions_toc.csv"), row.names = FALSE)
hdr("diagnosis residual R2 iter1c"); print(tb(I1, "diagnosis_residual_R2_toc.csv"), row.names = FALSE)
hdr("ignorance budget (09)"); print(cm("ignorance_budget.csv"), row.names = FALSE)
hdr("bootstrap budget (08k)"); b <- cm("iteration_bootstrap_budget.csv"); if (!is.null(b)) print(b[grepl("clean", b$first) | grepl("clean", b$second), ], row.names = FALSE)
hdr("bootstrap regions (08k)"); b <- cm("iteration_bootstrap_regions.csv"); if (!is.null(b)) print(b[grepl("clean", b$first) | grepl("clean", b$second), ], row.names = FALSE)
hdr("localization (08f)"); l <- cm("localization_by_scale.csv"); if (!is.null(l)) print(l[grepl("clean", l[[1]]), ], row.names = FALSE)
hdr("map change (09b)"); print(cm("toc_map_change_iter1_to_iter2.csv"), row.names = FALSE)
hdr("prediction artefact step test (08o)"); print(cm("prediction_artefact_step_test.csv"), row.names = FALSE)
hdr("prediction artefact by distance (08o)"); print(cm("prediction_artefact_by_distance.csv"), row.names = FALSE)
hdr("lithology edges check (08d)"); print(cm("lithology_edges_check.csv"), row.names = FALSE)
hdr("lithology edges realised error (08e)"); print(cm("lithology_edges_realised_error.csv"), row.names = FALSE)
for (it in c(I1, I2)) { hdr(paste("regionalization (08g)", it)); print(tb(it, "regionalization_by_region_toc.csv"), row.names = FALSE) }
for (it in c(I1, I2)) {   # F5.3: the collection test was repeated on the final map (iter2c, 10 random starts)
  hdr(paste("targeted sampling pooled (08h)", it)); print(tb(it, "targeted_sampling_pooled_toc.csv"), row.names = FALSE)
  hdr(paste("collection by region (08l)", it)); print(tb(it, "collection_by_region_toc.csv"), row.names = FALSE)
  hdr(paste("space-environment structure (08n)", it)); print(tb(it, "space_env_structure_toc.csv"), row.names = FALSE)
}
hdr("F5.4 representative novelty in the retrospective test (08ae)"); print(cm("f54_repnovelty_pooled_toc.csv"), row.names = FALSE)
hdr("F5.1 region-definition sensitivity (08ab)"); print(cm("f51_region_sensitivity_stability.csv"), row.names = FALSE)
hdr("F5.2 stratum interaction (08ac)"); print(cm("f52_stratum_interaction_gain.csv"), row.names = FALSE)
hdr("sampling priority: summary (08aa, F2)"); print(cm("sampling_priority_summary_toc.csv"), row.names = FALSE)
hdr("sampling priority: artefact check of the first sites (08aa, F2)"); print(cm("sampling_priority_artefact_check_toc.csv"), row.names = FALSE)
hdr("sampling priority: by region (08aa, F2)"); print(cm("sampling_priority_by_region_toc.csv"), row.names = FALSE)
hdr("sampling priority: first 10 sites (08aa, F2)")
sp <- cm("sampling_priority_sites_toc.csv")
if (!is.null(sp)) print(head(sp[, c("rank", "lon", "lat", "depth_m", "toc_predicted_pct", "DInorm_initial", "gain",
  "distance_nearest_observation_km", "region", "collection", "flag_outside_region_support")], 10), row.names = FALSE)
hdr("F5.3 pooled: iteration and number of starts (08ad)"); print(cm("f53_targeted_sampling_pooled_compare.csv"), row.names = FALSE)
hdr("F5.3 decision and strategy by region (08ad)"); print(cm("f53_collection_by_region_compare.csv"), row.names = FALSE)
hdr("environmental space PCA (08i) iter1c"); print(tb(I1, "environmental_space_pca_toc.csv"), row.names = FALSE)
options(op); sink()
msg("clean chain numbers written to %s", out)
