# =============================================================================
# 05_transfer.R — leave-region-out: does error grow with DI, and does the
#                 ignorance map predict the error it will meet?
# =============================================================================
# For each of CFG$n_regions regions (k-means on 3-D coordinates of the observation
# cells): refit the model without the region (same predictors and hyper-parameters,
# spatial-block CV on the remaining data), rebuild its AOA and error profile, then
# predict the withheld region. Two tests:
#  1. Realised absolute error vs DI / threshold in the withheld region (inside vs outside AOA).
#  2. Calibration of the ignorance map: expected RMSE (from the error profile fitted
#     without the region) vs realised RMSE, in bins of expected RMSE. This is the
#     falsifiable prediction of the manuscript: flagged-ignorant places should be where
#     new data change the inference most.

source("R/00_config.R")
suppressPackageStartupMessages({ library(caret); library(ranger); library(CAST) })
set.seed(CFG$seed)

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster

res <- list()
for (r in sort(unique(d$region))) {
  trn <- d[d$region != r, ]; tst <- d[d$region == r, ]
  f <- CreateSpacetimeFolds(trn, spacevar = "block", k = CFG$cv_k, seed = CFG$seed)
  set.seed(CFG$seed)
  m <- train(trn[, vars], trn$y, method = "ranger", metric = "RMSE",
             trControl = trainControl(method = "cv", index = f$index, indexOut = f$indexOut, savePredictions = "final"),
             tuneGrid = bt, num.trees = 300, importance = "permutation")
  tdi <- trainDI(m, verbose = FALSE, algorithm = "kd_tree")
  a   <- aoa(tst[, vars], model = m, trainDI = tdi, verbose = FALSE, algorithm = "kd_tree")
  ep  <- tryCatch(errorProfiles(m, tdi, variable = "DI", calib = "scam"), error = function(e) NULL)
  p   <- predict(m$finalModel, data = tst[, vars], num.threads = CFG$cores)$predictions
  exp_rmse <- if (is.null(ep)) NA_real_ else as.numeric(predict(ep, data.frame(DI = a$DI)))
  exp_rmse[a$AOA == 0] <- NA_real_      # the profile is only defined inside the AOA
  res[[length(res) + 1]] <- data.frame(region = r, cell = tst$cell, lon = tst$x_lon, lat = tst$y_lat, obs = tst$y, pred = p,
                                       DI = a$DI, DInorm = a$DI / tdi$threshold, AOA = a$AOA, expected_RMSE = exp_rmse)
  msg("region %d: n=%d, %.0f%% outside AOA, RMSE %.3f (inside %.3f, outside %.3f)", r, nrow(tst), 100 * mean(a$AOA == 0),
      sqrt(mean((p - tst$y)^2)), sqrt(mean((p - tst$y)[a$AOA == 1]^2)), sqrt(mean((p - tst$y)[a$AOA == 0]^2)))
}
tt <- bind_rows(res); tt$err <- tt$pred - tt$obs; tt$abs_err <- abs(tt$err)
write.csv(tt, file.path(DIRS$tables, "transfer_points_toc.csv"), row.names = FALSE)

rmse <- function(e) sqrt(mean(e^2))
boot <- function(e, B = 1000) quantile(replicate(B, rmse(sample(e, replace = TRUE))), c(0.025, 0.975))
regions <- tt %>% group_by(region) %>% summarise(lon = atan2(mean(sin(lon * pi / 180)), mean(cos(lon * pi / 180))) * 180 / pi,
  lat = mean(lat), n = n(), pct_outside_AOA = 100 * mean(AOA == 0), DInorm_median = median(DInorm),
  RMSE = rmse(err), RMSE_inside = rmse(err[AOA == 1]), RMSE_outside = ifelse(any(AOA == 0), rmse(err[AOA == 0]), NA_real_),
  R2 = 1 - sum(err^2) / sum((obs - mean(obs))^2), .groups = "drop")
write.csv(regions, file.path(DIRS$tables, "transfer_regions_toc.csv"), row.names = FALSE)

tt$bin <- cut(rank(tt$DInorm, ties.method = "first"), 10, labels = FALSE)
bins <- do.call(rbind, lapply(split(tt, tt$bin), function(b) { ci <- boot(b$err)
  data.frame(bin = b$bin[1], DInorm_median = median(b$DInorm), n = nrow(b), pct_outside_AOA = 100 * mean(b$AOA == 0),
             RMSE = rmse(b$err), RMSE_lo = ci[[1]], RMSE_hi = ci[[2]], MAE = mean(b$abs_err)) }))
write.csv(bins, file.path(DIRS$tables, "transfer_DI_bins_toc.csv"), row.names = FALSE)

cal <- tt[!is.na(tt$expected_RMSE), ]
cal$bin <- cut(rank(cal$expected_RMSE, ties.method = "first"), 10, labels = FALSE)
calib <- do.call(rbind, lapply(split(cal, cal$bin), function(b) { ci <- boot(b$err)
  data.frame(bin = b$bin[1], expected_RMSE = mean(b$expected_RMSE), n = nrow(b), realised_RMSE = rmse(b$err), lo = ci[[1]], hi = ci[[2]]) }))
write.csv(calib, file.path(DIRS$tables, "transfer_ignorance_calibration_toc.csv"), row.names = FALSE)

sp <- function(a, b) { ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)); c(unname(ct$estimate), ct$p.value) }
ci_in <- boot(tt$err[tt$AOA == 1]); ci_out <- boot(tt$err[tt$AOA == 0])
S <- function(metric, value, description) data.frame(metric, value = as.character(value), description)
add <- rbind(
  S("transfer_n_regions", nrow(regions), "k-means on 3-D coordinates"),
  S("transfer_pct_outside_AOA", round(100 * mean(tt$AOA == 0), 2), "withheld points outside the AOA of the model fitted without them"),
  S("transfer_RMSE_inside_AOA", sprintf("%.3f [%.3f-%.3f]", rmse(tt$err[tt$AOA == 1]), ci_in[1], ci_in[2]), "bootstrap 95% CI"),
  S("transfer_RMSE_outside_AOA", sprintf("%.3f [%.3f-%.3f]", rmse(tt$err[tt$AOA == 0]), ci_out[1], ci_out[2]), "bootstrap 95% CI"),
  S("transfer_spearman_DI_abs_error_rho", round(sp(tt$DInorm, tt$abs_err)[1], 3), "withheld points"),
  S("transfer_spearman_DI_abs_error_p", signif(sp(tt$DInorm, tt$abs_err)[2], 3), ""),
  S("transfer_spearman_binDI_RMSE_rho", round(sp(bins$DInorm_median, bins$RMSE)[1], 3), "10 DI/threshold deciles"),
  S("transfer_spearman_regions_pctOutside_RMSE_rho", round(sp(regions$pct_outside_AOA, regions$RMSE)[1], 3), "regions"),
  S("calibration_spearman_expected_realised_rho", round(sp(calib$expected_RMSE, calib$realised_RMSE)[1], 3), "10 bins of expected RMSE (inside AOA)"),
  S("calibration_slope_realised_on_expected", round(coef(lm(realised_RMSE ~ expected_RMSE, calib))[2], 3), "1 = perfectly calibrated")
)
f <- file.path(DIRS$tables, "toc_summary.csv"); summ <- read.csv(f, colClasses = "character")
write.csv(rbind(summ[!grepl("^(transfer_|calibration_)", summ$metric), ], add), f, row.names = FALSE)
print(as.data.frame(regions), digits = 3); print(bins, digits = 3); print(calib, digits = 3); print(add[, 1:2])
