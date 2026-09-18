# =============================================================================
# 04_predict_aoa.R — global TOC prediction, DI, AOA and a calibrated ignorance map
# =============================================================================
# Ignorance map: CAST::errorProfiles fits a monotone (scam) relation between the
# dissimilarity index of the cross-validated predictions and their error. It is applied
# to the DI of every ocean cell inside the AOA, giving the expected RMSE (log10 TOC).
# Outside the AOA no error can be estimated: those cells are flagged, not coloured.

source("R/00_config.R")
suppressPackageStartupMessages({ library(caret); library(ranger); library(CAST) })
set.seed(CFG$seed)

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
stk   <- rast(path_toc_stack())
tmpl  <- stk[["sst_mean"]]
X     <- values(stk[[vars]], dataframe = TRUE)
cells <- which(complete.cases(X)); X <- X[cells, , drop = FALSE]
chunks <- split(seq_along(cells), cut(seq_along(cells), ceiling(length(cells) / 5e5), labels = FALSE))
msg("%d ocean cells; predictors: %s", length(cells), paste(vars, collapse = ", "))

pred <- numeric(length(cells))
for (ch in chunks) pred[ch] <- predict(model$finalModel, data = X[ch, , drop = FALSE], num.threads = CFG$cores)$predictions

tdi <- trainDI(model, verbose = FALSE, algorithm = "kd_tree")
saveRDS(tdi, file.path(DIRS$models, "trainDI_toc.rds"))
DI <- numeric(length(cells)); AOA <- integer(length(cells))
for (ch in chunks) { a <- aoa(X[ch, , drop = FALSE], model = model, trainDI = tdi, verbose = FALSE, algorithm = "kd_tree"); DI[ch] <- a$DI; AOA[ch] <- a$AOA }
msg("AOA threshold = %.4f", tdi$threshold)

ep <- errorProfiles(model, tdi, variable = "DI", calib = "scam")
saveRDS(ep, file.path(DIRS$models, "errorProfile_toc.rds"))
expRMSE <- rep(NA_real_, length(cells))
inA <- AOA == 1
expRMSE[inA] <- as.numeric(predict(ep, data.frame(DI = DI[inA])))

to_r <- function(v, name) { r <- rep(NA_real_, ncell(tmpl)); r[cells] <- v; r <- setValues(tmpl, r); names(r) <- name; r }
out <- function(r, stem, dt = "FLT4S") writeRaster(r, file.path(DIRS$rasters, sprintf("toc_%s.tif", stem)), overwrite = TRUE, datatype = dt, gdal = "COMPRESS=DEFLATE")
out(to_r(10^pred - CFG$toc_offset, "toc_pct"), "prediction_pct")
out(to_r(pred, "log10_toc"), "prediction_log10")
out(to_r(DI, "DI"), "DI"); out(to_r(AOA, "AOA"), "AOA", "INT1U")
out(to_r(expRMSE, "expected_RMSE_log10"), "ignorance_expectedRMSE")
obs_r <- rasterize(cbind(d$x_lon, d$y_lat), tmpl, field = 1, fun = "sum")
out(obs_r, "observation_cells", "INT2U")

# --- statistics --------------------------------------------------------------------------------
area  <- values(cellSize(tmpl, unit = "km"))[cells, 1]
depth <- values(stk[["depth"]])[cells, 1]
lat   <- yFromCell(tmpl, cells)
pct_out <- function(m) 100 * sum(area[m & AOA == 0]) / sum(area[m])
dcls <- cut(depth, c(-Inf, 200, 1000, 3000, Inf), labels = c("shelf <200 m", "slope 200-1000 m", "1000-3000 m", "abyss >3000 m"))
by_depth <- tapply(seq_along(cells), dcls, function(i) round(pct_out(seq_along(cells) %in% i), 1))

pr <- model$pred; pr <- pr[order(pr$rowIndex), ]
pr$DI <- tdi$trainDI[pr$rowIndex]; pr$abs_err <- abs(pr$pred - pr$obs)
rho <- cor.test(pr$DI, pr$abs_err, method = "spearman", exact = FALSE)
pr$bin <- cut(rank(pr$DI, ties.method = "first"), 10, labels = FALSE)
boot_ci <- function(e, fun, B = 1000) { v <- replicate(B, fun(sample(e, replace = TRUE))); quantile(v, c(0.025, 0.975)) }
bins <- do.call(rbind, lapply(split(pr, pr$bin), function(b) {
  ciR <- boot_ci(b$pred - b$obs, function(z) sqrt(mean(z^2)))
  data.frame(bin = b$bin[1], DI_median = median(b$DI), DInorm_median = median(b$DI) / tdi$threshold, n = nrow(b),
             RMSE = sqrt(mean((b$pred - b$obs)^2)), RMSE_lo = ciR[[1]], RMSE_hi = ciR[[2]],
             MAE = mean(b$abs_err), R2 = 1 - sum((b$pred - b$obs)^2) / sum((b$obs - mean(b$obs))^2)) }))
write.csv(bins, file.path(DIRS$tables, "DI_vs_CVerror_bins_toc.csv"), row.names = FALSE)
write.csv(pr[, c("rowIndex", "obs", "pred", "DI", "abs_err", "Resample")], file.path(DIRS$tables, "DI_vs_CVerror_points_toc.csv"), row.names = FALSE)
rho_bins <- cor.test(bins$DI_median, bins$RMSE, method = "spearman", exact = FALSE)
above <- pr$DI > tdi$threshold

cvperf <- read.csv(file.path(DIRS$tables, "cv_performance_toc.csv"))
gd <- read.csv(file.path(DIRS$tables, "geodist_toc.csv")); gmed <- tapply(gd$dist_km, gd$what, median)
S <- function(metric, value, description) data.frame(metric, value = as.character(value), description)
summ <- rbind(
  S("n_observation_cells", nrow(d), "0.1-deg cells with TOC observations"),
  S("predictors_selected", paste(vars, collapse = ";"), "CAST::ffs, spatial CV RMSE"),
  S("R2_log_spatialCV", round(cvperf$R2_log[cvperf$cv == "spatial_blocks"], 3), "pooled out-of-fold, 1000-km blocks"),
  S("R2_log_randomCV", round(cvperf$R2_log[cvperf$cv == "random"], 3), "random 5-fold CV (optimistic)"),
  S("RMSE_log_spatialCV", round(cvperf$RMSE_log[cvperf$cv == "spatial_blocks"], 3), "log10 TOC units"),
  S("noise_floor_withincell_sd_log", round(cvperf$noise_floor_withincell_sd_log[1], 3), "median SD of log10 TOC within 0.1-deg cells with >= 3 records"),
  S("AOA_DI_threshold", round(tdi$threshold, 4), ""),
  S("pct_ocean_outside_AOA", round(pct_out(rep(TRUE, length(cells))), 2), "area-weighted, whole ocean"),
  S("pct_outside_AOA_shelf", by_depth[["shelf <200 m"]], ""), S("pct_outside_AOA_slope", by_depth[["slope 200-1000 m"]], ""),
  S("pct_outside_AOA_1000_3000", by_depth[["1000-3000 m"]], ""), S("pct_outside_AOA_abyss", by_depth[["abyss >3000 m"]], ""),
  S("pct_outside_AOA_south", round(pct_out(lat < 0), 2), "Southern Hemisphere"), S("pct_outside_AOA_north", round(pct_out(lat >= 0), 2), "Northern Hemisphere"),
  S("spearman_DI_abs_error_CV_rho", round(unname(rho$estimate), 3), "CV points"), S("spearman_DI_abs_error_CV_p", signif(rho$p.value, 3), ""),
  S("spearman_binDI_RMSE_CV_rho", round(unname(rho_bins$estimate), 3), "10 DI deciles"),
  S("RMSE_CV_below_threshold", round(sqrt(mean((pr$pred - pr$obs)[!above]^2)), 3), sprintf("n = %d", sum(!above))),
  S("RMSE_CV_above_threshold", round(sqrt(mean((pr$pred - pr$obs)[above]^2)), 3), sprintf("n = %d", sum(above))),
  S("expected_RMSE_ocean_median", round(median(expRMSE, na.rm = TRUE), 3), "errorProfiles, inside AOA"),
  S("median_dist_CV_km", round(gmed[["CV test-to-train"]]), ""), S("median_dist_prediction_km", round(gmed[["prediction-to-sample (global ocean)"]]), ""),
  S("seed", CFG$seed, "")
)
write.csv(summ, file.path(DIRS$tables, "toc_summary.csv"), row.names = FALSE)
print(summ[, 1:2]); print(bins, digits = 3)
writeLines(capture.output(sessioninfo::session_info()), file.path(ROOT, "outputs", "sessionInfo.txt"))
