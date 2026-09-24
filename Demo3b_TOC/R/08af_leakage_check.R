# =============================================================================
# 08af_leakage_check.R — MSv8 (author comments C67/C76): what does the nested calibration buy?
# =============================================================================
# The ignorance map is calibrated in a NESTED leave-one-region-out design: for each withheld region
# the forest, the dissimilarity normalisation and the error model are rebuilt from the other regions
# only (07_calibration.R). The cheaper alternative is to calibrate once on all data, with the errors
# of the spatial-block cross-validation. Then the data of a region inform the error model that will
# be read in that same region, and the dissimilarity of its cells is computed against themselves:
# two forms of spatial leakage. This script measures both, at the same points and against the same
# realised errors (the error each cell had when its region was withheld):
#   (1) nested          the calibration used for the map            (exp_C3 of the outer points)
#   (2) error model only the C3 model fitted on the block-CV errors of ALL data, evaluated at the
#                       nested dissimilarity and distance of those points
#   (3) as the map is read at those cells: the final map's expected error, where the dissimilarity
#                       is computed with the region's own data in the training set
# Functions fit_rf, cv_errors, fit_cal and exp_rmse are copied from 07_calibration.R so that the
# numbers come from the same code path.
#   DEMO3B_ITER=iter2c_clean Rscript R/08af_leakage_check.R
# Writes outputs/comparison/leakage_check_*.csv; nothing is overwritten.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(CAST); library(FNN); library(scam); library(terra) })
set.seed(CFG$seed)

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
d$depth <- pmax(d$depth, 0)
to_xyz   <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))
MAX_CAL <- 12000

fit_rf <- function(df, imp = "none") ranger(x = df[, vars, drop = FALSE], y = df$y, num.trees = 300,
  mtry = min(bt$mtry, length(vars)), min.node.size = bt$min.node.size, splitrule = "variance",
  importance = imp, num.threads = CFG$cores, seed = CFG$seed)

cv_errors <- function(df, fold_id, w, label) {          # 07_calibration.R
  ks <- sort(unique(fold_id)); pred <- rep(NA_real_, nrow(df))
  for (k in ks) { tr <- fold_id != k; m <- fit_rf(df[tr, ])
    pred[!tr] <- predict(m, df[!tr, vars, drop = FALSE], num.threads = CFG$cores)$predictions }
  CVtest <- lapply(ks, function(k) which(fold_id == k)); CVtrain <- lapply(ks, function(k) which(fold_id != k))
  tdi <- suppressMessages(trainDI(train = df[, vars, drop = FALSE], variables = vars, weight = w,
                                  CVtest = CVtest, CVtrain = CVtrain, verbose = FALSE))
  X <- to_xyz(df$x_lon, df$y_lat); dist <- numeric(nrow(df))
  for (i in seq_along(ks)) dist[CVtest[[i]]] <- chord_km(get.knnx(X[CVtrain[[i]], , drop = FALSE], X[CVtest[[i]], , drop = FALSE], k = 1)$nn.dist[, 1])
  data.frame(source = label, e2 = (pred - df$y)^2, DI = tdi$trainDI, dist_km = dist, pred = pred, depth = df$depth, alat = abs(df$y_lat))
}

fit_cal_C3 <- function(cal) {                            # 07_calibration.R, method C3
  cal <- cal[is.finite(cal$DI) & is.finite(cal$e2), ]
  if (nrow(cal) > MAX_CAL) cal <- cal[sample(nrow(cal), MAX_CAL), ]
  cal$e2 <- pmax(cal$e2, 1e-6); cal$ldist <- log10(cal$dist_km + 1)
  f <- e2 ~ s(DI, bs = "mpi") + s(ldist, bs = "mpi") + s(pred, bs = "tp", k = 6)
  fit <- tryCatch(scam(f, family = Gamma(link = "log"), data = cal), error = function(e) NULL)
  if (is.null(fit)) fit <- scam(f, family = Gamma(link = "log"), data = cal[cal$e2 < quantile(cal$e2, 0.999), ])
  fit
}
exp_rmse <- function(fit, nd) { nd$ldist <- log10(nd$dist_km + 1); sqrt(as.numeric(predict(fit, nd, type = "response"))) }

## ---- the three calibrations at the same points --------------------------------------------------
pts <- read.csv(file.path(DIRS$tables, "calibration_outer_points_toc.csv"))
m_all <- fit_rf(d, "permutation")
w_all <- as.data.frame(t(pmax(m_all$variable.importance[vars], 0)))
cal_block <- cv_errors(d, d$fold, w_all, "block")       # block CV on ALL data: the leaky calibration set
set.seed(CFG$seed)
fit_single <- fit_cal_C3(cal_block)
pts$exp_single <- exp_rmse(fit_single, data.frame(DI = pts$DI, dist_km = pts$dist_km, pred = pts$pred))
r_map <- rast(file.path(DIRS$rasters, "toc_ignorance_calibrated_expectedRMSE.tif"))
pts$exp_map <- extract(r_map, as.matrix(pts[, c("lon", "lat")]))[, 1]
pts <- pts[is.finite(pts$exp_C3) & is.finite(pts$exp_single) & is.finite(pts$exp_map), ]
msg("%d calibration points", nrow(pts))

realised <- sqrt(pts$e2)
row <- function(name, ex) data.frame(calibration = name, n = nrow(pts), expected_RMSE = sqrt(mean(ex^2)),
  realised_RMSE = sqrt(mean(pts$e2)), ratio_realised_expected = sqrt(mean(pts$e2)) / sqrt(mean(ex^2)),
  coverage90 = mean(realised <= 1.645 * ex),
  regions_covered_within_5_points = sum(tapply(seq_len(nrow(pts)), pts$region, function(i)
    abs(mean(realised[i] <= 1.645 * ex[i]) - 0.90) <= 0.05)))
out <- rbind(row("nested leave-one-region-out (the map)", pts$exp_C3),
             row("error model fitted on block-CV errors of all data", pts$exp_single),
             row("map read at these cells (dissimilarity also from all data)", pts$exp_map))
write.csv(out, file.path(DIRS$compare, "leakage_check_toc.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)

byreg <- do.call(rbind, lapply(split(seq_len(nrow(pts)), pts$region), function(i) data.frame(
  region = pts$region[i][1], n = length(i), realised_RMSE = sqrt(mean(pts$e2[i])),
  expected_nested = sqrt(mean(pts$exp_C3[i]^2)), expected_single = sqrt(mean(pts$exp_single[i]^2)),
  expected_map = sqrt(mean(pts$exp_map[i]^2)),
  coverage90_nested = mean(realised[i] <= 1.645 * pts$exp_C3[i]),
  coverage90_single = mean(realised[i] <= 1.645 * pts$exp_single[i]),
  coverage90_map = mean(realised[i] <= 1.645 * pts$exp_map[i]))))
write.csv(byreg, file.path(DIRS$compare, "leakage_check_by_region_toc.csv"), row.names = FALSE)
print(byreg, digits = 3, row.names = FALSE)
