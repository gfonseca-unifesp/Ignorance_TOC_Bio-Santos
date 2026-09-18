# =============================================================================
# 08j_arctic_check.R — why did carbon-supply layers worsen transfer into the Siberian Arctic?
# =============================================================================
# Candidate cause: satellite surface POC and river-input fields are least constrained at high latitude
# (sea ice, river-borne coloured organic matter, very large Siberian rivers). Local proxies only (no new data):
#   sst_min (ice-season proxy), sss_mean (river freshening), river_poc, river_tss, poc_surf, depth.
# Tests on the withheld Siberian Arctic points (nested leave-region-out, 07_calibration.R):
#   (1) change in absolute error (iteration 1b - iteration 1) vs each proxy (Spearman; tercile means);
#   (2) share of Arctic points whose supply values fall outside the range of the other regions' data
#       (the model trained without the Arctic must extrapolate there).
Sys.setenv(DEMO3B_ITER = "iter1b_supply")
source("R/00_config.R")
key <- function(lon, lat) paste(round(lon, 3), round(lat, 3))
rd  <- function(it) read.csv(file.path(ROOT, "outputs", it, "tables", "calibration_outer_points_toc.csv"))
p1 <- rd("iter1_baseline"); p2 <- rd("iter1b_supply")
a <- data.frame(k = key(p1$lon, p1$lat), region = p1$region, lon = p1$lon, lat = p1$lat, ae1 = abs(p1$pred - p1$obs))
b <- data.frame(k = key(p2$lon, p2$lat), ae2 = abs(p2$pred - p2$obs), err2 = p2$pred - p2$obs)
p <- merge(a, b, by = "k")
reg <- read.csv(file.path(ROOT, "outputs", "iter1_baseline", "tables", "transfer_regions_toc.csv"))
arc <- reg$region[which.max(reg$lat)]
stk <- rast(path_toc_stack())
prox <- c("sst_min", "sss_mean", "river_poc", "river_tss", "poc_surf", "poc_flux", "depth", "dist_coast_km")
X <- stk[[prox]][cellFromXY(stk, cbind(p$lon, p$lat))]
p <- cbind(p, X)
A <- p[p$region == arc, ]; O <- p[p$region != arc, ]
msg("Arctic region %d: %d points; mean abs error 1 -> 1b: %.3f -> %.3f", arc, nrow(A), mean(A$ae1), mean(A$ae2))
A$dae <- A$ae2 - A$ae1
t1 <- do.call(rbind, lapply(prox, function(v) {
  q <- quantile(A[[v]], c(1/3, 2/3), na.rm = TRUE); g <- cut(A[[v]], c(-Inf, q, Inf), labels = c("low", "mid", "high"))
  m <- tapply(A$dae, g, mean)
  lo <- quantile(O[[v]], 0.01, na.rm = TRUE); hi <- quantile(O[[v]], 0.99, na.rm = TRUE)
  data.frame(proxy = v, rho_change_abs_error = cor(A[[v]], A$dae, method = "spearman", use = "complete.obs"),
             rho_signed_error_1b = cor(A[[v]], A$err2, method = "spearman", use = "complete.obs"),
             mean_change_low = m[["low"]], mean_change_mid = m[["mid"]], mean_change_high = m[["high"]],
             pct_arctic_outside_other_regions_1_99 = 100 * mean(A[[v]] < lo | A[[v]] > hi, na.rm = TRUE))
}))
write.csv(t1, file.path(DIRS$tables, "arctic_check_toc.csv"), row.names = FALSE)
print(t1, digits = 3, row.names = FALSE)
msg("mean signed error in the Arctic, iteration 1b: %.3f (positive = overprediction)", mean(A$err2))
msg("arctic check done")
