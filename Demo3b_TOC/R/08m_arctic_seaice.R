# =============================================================================
# 08m_arctic_seaice.R — does sea-ice cover explain why supply layers worsened the Siberian Arctic?
# =============================================================================
# Sea-ice concentration: Bio-ORACLE v3, siconc_baseline_2000_2020_depthsurf, variable siconc_mean
# (two decadal steps, averaged), ERDDAP griddap subset >= 60 N at stride 2 (0.1 deg).
# Tests on the withheld Siberian Arctic points (07_calibration.R, iteration 1 vs 1b):
#   (1) Spearman of sea-ice cover with the change in absolute error (1b - 1) and with the signed error of 1b;
#   (2) error change and signed error by sea-ice class (ice-free < 0.15, seasonal 0.15-0.5, persistent >= 0.5);
#   (3) surface-POC signal by ice class: is satellite POC lower under ice for the same observed TOC?
Sys.setenv(DEMO3B_ITER = "iter1b_supply")
source("R/00_config.R")
options(timeout = 1800)
f_nc <- file.path(DIRS$raw, "siconc_mean_baseline_2000_2020_arctic60N_0.1deg.nc")
if (!file.exists(f_nc)) {
  u <- paste0("https://erddap.bio-oracle.org/erddap/griddap/siconc_baseline_2000_2020_depthsurf.nc?",
              "siconc_mean%5B(2000-01-01T00:00:00Z):1:(2010-01-01T00:00:00Z)%5D%5B(60.025):2:(89.975)%5D%5B(-179.975):2:(179.975)%5D")
  dir.create(DIRS$raw, showWarnings = FALSE, recursive = TRUE)
  utils::download.file(u, f_nc, mode = "wb", quiet = TRUE, method = "libcurl")
}
msg("sea-ice file: %.1f MB", file.size(f_nc) / 1e6)
ice <- rast(f_nc)
if (nlyr(ice) > 1) ice <- mean(ice, na.rm = TRUE)
crs(ice) <- "EPSG:4326"
stk <- rast(path_toc_stack())

key <- function(lon, lat) paste(round(lon, 3), round(lat, 3))
rd  <- function(it) read.csv(file.path(ROOT, "outputs", it, "tables", "calibration_outer_points_toc.csv"))
p1 <- rd("iter1_baseline"); p2 <- rd("iter1b_supply")
p <- merge(data.frame(k = key(p1$lon, p1$lat), region = p1$region, lon = p1$lon, lat = p1$lat, obs = p1$obs, ae1 = abs(p1$pred - p1$obs)),
           data.frame(k = key(p2$lon, p2$lat), ae2 = abs(p2$pred - p2$obs), err2 = p2$pred - p2$obs), by = "k")
reg <- read.csv(file.path(ROOT, "outputs", "iter1_baseline", "tables", "transfer_regions_toc.csv"))
arc <- reg$region[which.max(reg$lat)]
A <- p[p$region == arc, ]
A$siconc <- extract(ice, cbind(A$lon, A$lat))[, 1]
A$poc_surf <- extract(stk[["poc_surf"]], cbind(A$lon, A$lat))[, 1]
A$dae <- A$ae2 - A$ae1
msg("Arctic points with sea-ice value: %d of %d; median siconc %.2f", sum(!is.na(A$siconc)), nrow(A), median(A$siconc, na.rm = TRUE))
A$ice_class <- cut(A$siconc, c(-Inf, 0.15, 0.5, Inf), labels = c("ice-free (<0.15)", "seasonal (0.15-0.5)", "persistent (>=0.5)"))
sp <- function(a, b) { ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)); c(rho = unname(ct$estimate), p = ct$p.value) }
tests <- rbind(
  data.frame(test = "rho(sea ice, change in abs error 1b-1)", t(sp(A$siconc, A$dae))),
  data.frame(test = "rho(sea ice, signed error 1b)", t(sp(A$siconc, A$err2))),
  data.frame(test = "rho(sea ice, surface POC)", t(sp(A$siconc, A$poc_surf))),
  data.frame(test = "rho(surface POC, observed TOC)", t(sp(A$poc_surf, A$obs))))
by_class <- do.call(rbind, lapply(split(A, A$ice_class), function(x) data.frame(ice_class = x$ice_class[1], n = nrow(x),
  mean_change_abs_error = mean(x$dae), mean_signed_error_1b = mean(x$err2), median_obs_log10TOC = median(x$obs),
  median_poc_surf = median(x$poc_surf, na.rm = TRUE))))
write.csv(tests, file.path(DIRS$tables, "arctic_seaice_tests_toc.csv"), row.names = FALSE)
write.csv(by_class, file.path(DIRS$tables, "arctic_seaice_by_class_toc.csv"), row.names = FALSE)
print(tests, digits = 3, row.names = FALSE); print(by_class, digits = 3, row.names = FALSE)
msg("arctic sea-ice check done")
