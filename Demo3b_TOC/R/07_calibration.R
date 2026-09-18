# =============================================================================
# 07_calibration.R — calibrating the ignorance map (expected error) and testing it
# =============================================================================
# Problem (step 05): the CAST error profile, fitted on 1000-km block CV, ranked the
# error of withheld regions only weakly and underestimated it (slope 0.27).
# Error models compared, each producing expected RMSE (log10 TOC) = sqrt(E[e^2]):
#   C0  monotone profile of DI, from block-CV errors (analogue of CAST::errorProfiles)
#   C1  monotone profile of DI, from leave-region-out CV errors (mimics predicting unsampled regions)
#   C2  monotone in DI + monotone in distance to the nearest observation, block + region errors
#   C3  C2 + smooth of the predicted TOC level (heteroscedasticity)
#   C4  DI + smooths of log depth and |latitude| (heteroscedastic by stratum; block + region errors)
#   C5  C2 + smooths of log depth and |latitude|
# C4-C5 follow the error-structure test (08w, 08x): depth zone and latitude band structure the MAGNITUDE of
# out-of-sample error (confirmed in 87-100% of half-splits) but not its sign, so they enter the error model,
# not the mean model. Depth and |latitude| are clamped to the range of the calibration data when predicting.
# Fitted as scam GAMs on squared errors (Gamma, log link; 'mpi' = monotone increasing).
#
# Honest test: nested leave-one-region-out. For each withheld region r, the RF, the DI
# normalisation and every error model are built from the other regions only. Expected and
# realised errors in r are then compared. Metrics, pooled over withheld points in 10 bins of
# expected RMSE: calibration slope/intercept, mean absolute calibration error (MACE),
# realised/expected ratio, Spearman (bins, points), and coverage of a nominal 90% interval
# (|e| <= 1.645 * expected RMSE). Stratum calibration: realised vs expected RMSE and coverage within each depth
# zone (4) and latitude band (3); MACE_strata = mean |realised - expected| over these 7 levels.
# Selection (env CAL_SELECT, default "strata_rule"): among the methods whose coverage stays within 5 points of
# nominal in EVERY stratum (0.85-0.95), the one with the lowest MACE_strata builds the final map from all data;
# if none qualifies, the one with the most even coverage. The rule asks two things of a map, in order: its stated
# interval must mean what it says wherever it is read, and its expected error must be right in level within
# strata. Pooled MACE (the earlier criterion) is reported alongside. Setting CAL_SELECT to a metric name
# (e.g. "MACE", "MACE_strata") selects by that metric instead.

source("R/00_config.R")
suppressPackageStartupMessages({ library(caret); library(ranger); library(CAST); library(FNN); library(scam) })
set.seed(CFG$seed)

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz   <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster   # same as 05
d$depth <- pmax(d$depth, 0)
METHODS <- c("C0", "C1", "C2", "C3", "C4", "C5")
MAX_CAL <- 12000
CAL_SELECT <- Sys.getenv("CAL_SELECT", "strata_rule")
zone <- function(z) as.character(cut(z, c(-Inf, 200, 1000, 3000, Inf), labels = c("shelf <=200 m", "slope 200-1000 m", "1000-3000 m", "abyss >3000 m")))
band <- function(lat) as.character(cut(abs(lat), c(0, 30, 60, 90), include.lowest = TRUE, labels = c("tropics 0-30", "mid-latitudes 30-60", "high latitudes 60-90")))

fit_rf <- function(df, imp = "none") ranger(x = df[, vars, drop = FALSE], y = df$y, num.trees = 300,
  mtry = min(bt$mtry, length(vars)), min.node.size = bt$min.node.size, splitrule = "variance",
  importance = imp, num.threads = CFG$cores, seed = CFG$seed)

cv_errors <- function(df, fold_id, w, label) {
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

fit_cal <- function(cal, method) {
  cal <- cal[is.finite(cal$DI) & is.finite(cal$e2), ]
  if (nrow(cal) > MAX_CAL) cal <- cal[sample(nrow(cal), MAX_CAL), ]
  cal$e2 <- pmax(cal$e2, 1e-6); cal$ldist <- log10(cal$dist_km + 1)
  if (method %in% c("C4", "C5")) { cal <- cal[is.finite(cal$depth) & is.finite(cal$alat), ]; cal$ldepth <- log10(cal$depth + 1) }
  f <- switch(method,
    C0 = e2 ~ s(DI, bs = "mpi"),
    C1 = e2 ~ s(DI, bs = "mpi"),
    C2 = e2 ~ s(DI, bs = "mpi") + s(ldist, bs = "mpi"),
    C3 = e2 ~ s(DI, bs = "mpi") + s(ldist, bs = "mpi") + s(pred, bs = "tp", k = 6),
    C4 = e2 ~ s(DI, bs = "mpi") + s(ldepth, bs = "tp", k = 5) + s(alat, bs = "tp", k = 5),
    C5 = e2 ~ s(DI, bs = "mpi") + s(ldist, bs = "mpi") + s(ldepth, bs = "tp", k = 5) + s(alat, bs = "tp", k = 5))
  fit <- tryCatch(scam(f, family = Gamma(link = "log"), data = cal), error = function(e) NULL)
  if (is.null(fit)) {   # numerical failure (non-finite coefficients): retry without the most extreme squared errors
    msg("  scam %s failed; retrying without the top 0.1%% of squared errors", method)
    fit <- tryCatch(scam(f, family = Gamma(link = "log"), data = cal[cal$e2 < quantile(cal$e2, 0.999), ]), error = function(e) NULL)
  }
  if (is.null(fit)) {   # second retry: smaller bases
    msg("  scam %s failed again; retrying with smaller bases", method)
    f5 <- switch(method,
      C0 = e2 ~ s(DI, bs = "mpi", k = 5),
      C1 = e2 ~ s(DI, bs = "mpi", k = 5),
      C2 = e2 ~ s(DI, bs = "mpi", k = 5) + s(ldist, bs = "mpi", k = 5),
      C3 = e2 ~ s(DI, bs = "mpi", k = 5) + s(ldist, bs = "mpi", k = 5) + s(pred, bs = "tp", k = 4),
      C4 = e2 ~ s(DI, bs = "mpi", k = 5) + s(ldepth, bs = "tp", k = 4) + s(alat, bs = "tp", k = 4),
      C5 = e2 ~ s(DI, bs = "mpi", k = 5) + s(ldist, bs = "mpi", k = 5) + s(ldepth, bs = "tp", k = 4) + s(alat, bs = "tp", k = 4))
    fit <- tryCatch(scam(f5, family = Gamma(link = "log"), data = cal[cal$e2 < quantile(cal$e2, 0.999), ]), error = function(e) NULL)
  }
  if (is.null(fit)) stop("calibration fit failed for ", method)
  if (method %in% c("C4", "C5")) attr(fit, "rng") <- list(ldepth = range(cal$ldepth), alat = range(cal$alat))
  fit
}
cal_subset <- function(cal, method) switch(method, C0 = cal[cal$source == "block", ], C1 = cal[cal$source == "region", ], cal)
exp_rmse <- function(fit, nd) {
  nd$ldist <- log10(nd$dist_km + 1)
  rg <- attr(fit, "rng")
  if (!is.null(rg)) {   # no extrapolation of the depth and latitude smooths beyond the calibration data
    nd$ldepth <- pmin(pmax(log10(pmax(nd$depth, 0) + 1), rg$ldepth[1]), rg$ldepth[2])
    nd$alat   <- pmin(pmax(nd$alat, rg$alat[1]), rg$alat[2])
  }
  sqrt(as.numeric(predict(fit, nd, type = "response")))
}

# --- nested leave-one-region-out --------------------------------------------------------------
f_outer <- file.path(DIRS$models, "calibration_outer_toc_v2.rds")   # v2 = with C4-C5 (v1 files, C0-C3 only, are kept)
if (!file.exists(f_outer)) {
  outer <- list()
  for (r in sort(unique(d$region))) {
    f_r <- file.path(DIRS$models, sprintf("calibration_outer_region%d_toc_v2.rds", r))   # checkpoint per withheld region
    if (file.exists(f_r)) { outer[[length(outer) + 1]] <- readRDS(f_r); msg("outer region %d: loaded checkpoint", r); next }
    t0 <- Sys.time()
    trn <- d[d$region != r, ]; tst <- d[d$region == r, ]
    m_out <- fit_rf(trn, "permutation")
    w <- as.data.frame(t(pmax(m_out$variable.importance[vars], 0)))
    ks <- sort(unique(trn$fold))
    tdi_o <- suppressMessages(trainDI(train = trn[, vars], variables = vars, weight = w,
                                      CVtest = lapply(ks, function(k) which(trn$fold == k)),
                                      CVtrain = lapply(ks, function(k) which(trn$fold != k)), verbose = FALSE))
    a <- suppressMessages(aoa(tst[, vars], trainDI = tdi_o, verbose = FALSE))
    p <- predict(m_out, tst[, vars], num.threads = CFG$cores)$predictions
    dist_t <- chord_km(get.knnx(to_xyz(trn$x_lon, trn$y_lat), to_xyz(tst$x_lon, tst$y_lat), k = 1)$nn.dist[, 1])
    cal <- rbind(cv_errors(trn, trn$fold, w, "block"), cv_errors(trn, trn$region, w, "region"))
    res <- data.frame(region = r, lon = tst$x_lon, lat = tst$y_lat, obs = tst$y, pred = p, e2 = (p - tst$y)^2,
                      DI = a$DI, AOA = a$AOA, dist_km = dist_t, depth = tst$depth)
    nd <- data.frame(DI = a$DI, dist_km = dist_t, pred = p, depth = tst$depth, alat = abs(tst$y_lat))
    for (mth in METHODS) { set.seed(CFG$seed + 100 * r + match(mth, METHODS))   # subsample of cal independent of method order
      res[[paste0("exp_", mth)]] <- exp_rmse(fit_cal(cal_subset(cal, mth), mth), nd) }
    outer[[length(outer) + 1]] <- res; saveRDS(res, f_r)
    msg("outer region %d done (%.1f min): realised RMSE %.3f | expected %s", r,
        as.numeric(difftime(Sys.time(), t0, units = "mins")), sqrt(mean(res$e2)),
        paste(sprintf("%s %.3f", METHODS, sapply(METHODS, function(m) sqrt(mean(res[[paste0("exp_", m)]]^2)))), collapse = " "))
  }
  saveRDS(outer, f_outer)
}
ev <- bind_rows(readRDS(f_outer))
write.csv(ev, file.path(DIRS$tables, "calibration_outer_points_toc.csv"), row.names = FALSE)

# --- evaluation -------------------------------------------------------------------------------------
set.seed(CFG$seed)
strata_of <- function(ev) list(depth_zone = zone(ev$depth), latitude_band = band(ev$lat))
evaluate <- function(ex, e2, label, strata) {
  st <- do.call(rbind, lapply(names(strata), function(pn) { g <- strata[[pn]]
    do.call(rbind, lapply(sort(unique(na.omit(g))), function(lv) { j <- which(g == lv)
      data.frame(method = label, partition = pn, level = lv, n = length(j), realised = sqrt(mean(e2[j])), expected = sqrt(mean(ex[j]^2)),
                 ratio = sqrt(mean(e2[j])) / sqrt(mean(ex[j]^2)), coverage90 = mean(sqrt(e2[j]) <= 1.645 * ex[j])) })) }))
  b <- cut(rank(ex, ties.method = "first"), 10, labels = FALSE)
  bb <- do.call(rbind, lapply(split(seq_along(ex), b), function(i) {
    ci <- quantile(replicate(500, sqrt(mean(sample(e2[i], replace = TRUE)))), c(0.025, 0.975))
    data.frame(method = label, bin = b[i][1], n = length(i), expected = sqrt(mean(ex[i]^2)), realised = sqrt(mean(e2[i])),
               lo = ci[[1]], hi = ci[[2]]) }))
  lmb <- lm(realised ~ expected, bb)
  list(bins = bb, strata = st, metrics = data.frame(method = label, slope = unname(coef(lmb)[2]), intercept = unname(coef(lmb)[1]),
    MACE = mean(abs(bb$realised - bb$expected)), ratio_realised_expected = sqrt(mean(e2)) / sqrt(mean(ex^2)),
    spearman_bins = cor(bb$expected, bb$realised, method = "spearman"),
    spearman_points = cor(ex, sqrt(e2), method = "spearman"),
    coverage90 = mean(sqrt(e2) <= 1.645 * ex),
    MACE_strata = mean(abs(st$realised - st$expected)), coverage90_strata_dev = mean(abs(st$coverage90 - 0.90)),
    coverage90_strata_min = min(st$coverage90), coverage90_strata_max = max(st$coverage90)))
}
evs <- lapply(METHODS, function(m) evaluate(ev[[paste0("exp_", m)]], ev$e2, m, strata_of(ev)))
metrics <- do.call(rbind, lapply(evs, `[[`, "metrics")); bins <- do.call(rbind, lapply(evs, `[[`, "bins"))
strata <- do.call(rbind, lapply(evs, `[[`, "strata"))
# also inside-AOA only
metrics_in <- do.call(rbind, lapply(METHODS, function(m) { i <- ev$AOA == 1; evaluate(ev[[paste0("exp_", m)]][i], ev$e2[i], m, strata_of(ev[i, ]))$metrics }))
metrics_in$subset <- "inside AOA"; metrics$subset <- "all withheld points"
write.csv(rbind(metrics, metrics_in), file.path(DIRS$tables, "calibration_metrics_toc.csv"), row.names = FALSE)
write.csv(bins, file.path(DIRS$tables, "calibration_bins_toc.csv"), row.names = FALSE)
write.csv(strata, file.path(DIRS$tables, "calibration_strata_toc.csv"), row.names = FALSE)
if (CAL_SELECT == "strata_rule") {
  ok <- metrics$coverage90_strata_min >= 0.85 & metrics$coverage90_strata_max <= 0.95
  best <- if (any(ok)) metrics$method[ok][which.min(metrics$MACE_strata[ok])] else metrics$method[which.min(metrics$coverage90_strata_dev)]
  msg("methods with coverage within 5 points of nominal in every stratum: %s", paste(metrics$method[ok], collapse = " "))
} else best <- metrics$method[which.min(metrics[[CAL_SELECT]])]
print(rbind(metrics, metrics_in), digits = 3)
print(reshape(strata[, c("method", "partition", "level", "coverage90")], idvar = c("partition", "level"), timevar = "method", direction = "wide"), digits = 2, row.names = FALSE)
msg("best calibration by %s: %s (lowest pooled MACE: %s; lowest MACE_strata: %s)", CAL_SELECT, best,
    metrics$method[which.min(metrics$MACE)], metrics$method[which.min(metrics$MACE_strata)])
per_region <- ev %>% group_by(region) %>% summarise(n = n(), realised = sqrt(mean(e2)),
  expected_best = sqrt(mean(.data[[paste0("exp_", best)]]^2)), expected_C0 = sqrt(mean(exp_C0^2)), .groups = "drop") %>%
  mutate(ratio_best = realised / expected_best, ratio_C0 = realised / expected_C0)
write.csv(per_region, file.path(DIRS$tables, "calibration_per_region_toc.csv"), row.names = FALSE)
print(as.data.frame(per_region), digits = 3)

# --- final calibrated ignorance map (all data) -------------------------------------------------------------
tdi_full <- readRDS(file.path(DIRS$models, "trainDI_toc.rds"))
w_full <- tdi_full$weight
cal_full <- rbind(cv_errors(d, d$fold, w_full, "block"), cv_errors(d, d$region, w_full, "region"))
fit_best <- fit_cal(cal_subset(cal_full, best), best)
saveRDS(list(method = best, fit = fit_best), file.path(DIRS$models, "ignorance_calibration_toc.rds"))

stk_tmpl <- rast(file.path(DIRS$rasters, "toc_DI.tif"))
DIv   <- values(stk_tmpl)[, 1]
predv <- values(rast(file.path(DIRS$rasters, "toc_prediction_log10.tif")))[, 1]
cells <- which(!is.na(DIv) & !is.na(predv))
xy <- xyFromCell(stk_tmpl, cells)
Xtr <- to_xyz(d$x_lon, d$y_lat)
dist <- numeric(length(cells))
for (ch in split(seq_along(cells), ceiling(seq_along(cells) / 5e5)))
  dist[ch] <- chord_km(get.knnx(Xtr, to_xyz(xy[ch, 1], xy[ch, 2]), k = 1)$nn.dist[, 1])
dep_r <- rast(path_toc_stack())[["depth"]]
stopifnot(compareGeom(dep_r, stk_tmpl, stopOnError = FALSE))
depv <- values(dep_r)[cells, 1]; depv[!is.finite(depv)] <- median(d$depth)
expv <- exp_rmse(fit_best, data.frame(DI = DIv[cells], dist_km = dist, pred = predv[cells], depth = depv, alat = abs(xy[, 2])))
r_exp <- setValues(stk_tmpl, NA_real_); r_exp[cells] <- expv; names(r_exp) <- "expected_RMSE_calibrated"
r_dst <- setValues(stk_tmpl, NA_real_); r_dst[cells] <- dist; names(r_dst) <- "dist_nearest_obs_km"
writeRaster(r_exp, file.path(DIRS$rasters, "toc_ignorance_calibrated_expectedRMSE.tif"), overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
writeRaster(r_dst, file.path(DIRS$rasters, "toc_distance_nearest_observation_km.tif"), overwrite = TRUE, gdal = "COMPRESS=DEFLATE")

area <- values(cellSize(stk_tmpl, unit = "km"))[cells, 1]
S <- function(metric, value, description) data.frame(metric, value = as.character(value), description)
mb <- metrics[metrics$method == best, ]; m0 <- metrics[metrics$method == "C0", ]
add <- rbind(
  S("calibration_best_method", best, sprintf("lowest %s in nested leave-region-out", CAL_SELECT)),
  S("calibration_best_slope", round(mb$slope, 3), "realised ~ expected over 10 bins (1 = calibrated)"),
  S("calibration_best_ratio", round(mb$ratio_realised_expected, 3), "realised / expected RMSE (1 = unbiased)"),
  S("calibration_best_coverage90", round(mb$coverage90, 3), "nominal 0.90"),
  S("calibration_best_spearman_bins", round(mb$spearman_bins, 3), ""),
  S("calibration_best_MACE", round(mb$MACE, 4), "pooled, 10 bins"),
  S("calibration_best_MACE_strata", round(mb$MACE_strata, 4), "mean |realised - expected| over 4 depth zones + 3 latitude bands"),
  S("calibration_best_coverage90_strata_range", sprintf("%.2f-%.2f", mb$coverage90_strata_min, mb$coverage90_strata_max), "nominal 0.90 in every stratum"),
  S("calibration_C0_MACE_strata", round(m0$MACE_strata, 4), ""),
  S("calibration_C0_coverage90_strata_range", sprintf("%.2f-%.2f", m0$coverage90_strata_min, m0$coverage90_strata_max), ""),
  S("calibration_C0_slope", round(m0$slope, 3), "block-CV DI profile (CAST-like)"),
  S("calibration_C0_ratio", round(m0$ratio_realised_expected, 3), ""),
  S("calibration_C0_coverage90", round(m0$coverage90, 3), ""),
  S("calibrated_expected_RMSE_ocean_median", round(median(expv), 3), "area-unweighted"),
  S("calibrated_expected_RMSE_ocean_p90", round(quantile(expv, 0.9), 3), ""),
  S("pct_ocean_expected_RMSE_gt_0.5", round(100 * sum(area[expv > 0.5]) / sum(area), 2), "expected error > factor ~3.2")
)
f <- file.path(DIRS$tables, "toc_summary.csv"); summ <- read.csv(f, colClasses = "character")
write.csv(rbind(summ[!grepl("^calibrat", summ$metric), ], add), f, row.names = FALSE)
print(add[, 1:2])

# --- figure -------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(tidyterra); library(patchwork); library(rnaturalearth) })
lab <- setNames(sprintf("%s: slope %.2f | ratio %.2f | cover90 %.2f | strata %.2f-%.2f", metrics$method, metrics$slope,
                        metrics$ratio_realised_expected, metrics$coverage90, metrics$coverage90_strata_min, metrics$coverage90_strata_max), metrics$method)
desc <- c(C0 = "C0  DI profile, block CV", C1 = "C1  DI profile, leave-region-out CV",
          C2 = "C2  DI + distance", C3 = "C3  DI + distance + TOC level",
          C4 = "C4  DI + depth + |latitude|", C5 = "C5  DI + distance + depth + |latitude|")[METHODS]
bins$panel <- factor(paste0(desc[bins$method], "\n", sub("^C[0-9]: ", "", lab[bins$method])), levels = paste0(desc, "\n", sub("^C[0-9]: ", "", lab[METHODS])))
lim <- range(c(bins$expected, bins$lo, bins$hi))
pc <- ggplot(bins, aes(expected, realised)) + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), colour = "#D7301F", size = 0.2) + geom_line(colour = "#D7301F") +
  facet_wrap(~panel, nrow = 2) + coord_equal(xlim = lim, ylim = lim) +
  labs(x = "expected RMSE (log10 TOC), fitted without the withheld region", y = "realised RMSE in withheld region",
       title = "a  Calibration of candidate ignorance maps (nested leave-region-out; dashed = perfect calibration)") +
  theme_bw(base_size = 7.5) + theme(strip.text = element_text(size = 6.5), plot.title = element_text(face = "bold", size = 9))

PROJ <- "+proj=eqearth +datum=WGS84 +units=m"
tp <- project(r_exp, PROJ, res = 25000)
rp <- project(r_exp, tp, method = "bilinear")
lims <- as.numeric(global(rp, function(z) quantile(z, c(0.02, 0.98), na.rm = TRUE))[1, ])
land <- st_transform(ne_countries(scale = 50, returnclass = "sf"), PROJ)
obs_sf <- st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), PROJ)
pm <- ggplot() + geom_spatraster(data = rp) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = lims, oob = scales::squish, na.value = NA,
                       name = sprintf("expected RMSE, log10 TOC (method %s; 2nd-98th pct stretch)", best)) +
  geom_sf(data = land, fill = "grey80", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = obs_sf, size = 0.03, colour = "#00C8FF", alpha = 0.4, inherit.aes = FALSE) +
  coord_sf(crs = PROJ, expand = FALSE, datum = NA) + theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.key.width = unit(14, "mm"), legend.key.height = unit(2.5, "mm"),
        plot.title = element_text(face = "bold", size = 9)) +
  ggtitle(sprintf("b  Calibrated ignorance map of seafloor TOC (cyan: observation cells). Nominal 90%% interval covers %.0f%% of withheld observations",
                  100 * mb$coverage90))
fig <- (pc / pm) + plot_layout(heights = c(1.3, 1.4))
ggsave(file.path(DIRS$figs, "toc_calibration.png"), fig, width = 260, height = 290, units = "mm", dpi = 300, bg = "white")

# stratum calibration: coverage and realised/expected ratio per depth zone and latitude band
strata$level <- factor(strata$level, levels = c("shelf <=200 m", "slope 200-1000 m", "1000-3000 m", "abyss >3000 m",
                                               "tropics 0-30", "mid-latitudes 30-60", "high latitudes 60-90"))
strata$pn <- ifelse(strata$partition == "depth_zone", "depth zone", "latitude band")
strata$hl <- ifelse(strata$method == best, sprintf("%s (map)", best), ifelse(strata$method == "C0", "C0 (DI only)", "other"))
cols <- setNames(c("#D7301F", "#2C7FB8", "grey70"), c(sprintf("%s (map)", best), "C0 (DI only)", "other"))
if (best == "C0") cols <- cols[c(1, 3)]
ps <- lapply(c("coverage90", "ratio"), function(v) ggplot(strata, aes(level, .data[[v]], group = method, colour = hl)) +
  geom_hline(yintercept = if (v == "ratio") 1 else 0.9, linetype = 2, colour = "grey40") +
  geom_line(data = strata[strata$hl == "other", ], linewidth = 0.3) + geom_point(data = strata[strata$hl == "other", ], size = 0.8) +
  geom_line(data = strata[strata$hl != "other", ], linewidth = 0.7) + geom_point(data = strata[strata$hl != "other", ], size = 1.6) +
  geom_text(data = strata[strata$level %in% c("abyss >3000 m", "high latitudes 60-90"), ], aes(label = method), hjust = -0.3, size = 2, show.legend = FALSE) +
  facet_grid(~pn, scales = "free_x", space = "free_x") + scale_colour_manual(values = cols, name = NULL) +
  labs(x = NULL, y = if (v == "ratio") "realised / expected RMSE" else "coverage of nominal 90% interval") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top"))
fs <- (ps[[1]] | ps[[2]]) + plot_annotation(title = sprintf("Calibration within strata (nested leave-region-out; selection by %s)", CAL_SELECT),
                                            theme = theme(plot.title = element_text(face = "bold", size = 9)))
ggsave(file.path(DIRS$figs, "toc_calibration_strata.png"), fs, width = 230, height = 100, units = "mm", dpi = 300, bg = "white")
msg("calibration figures written")
