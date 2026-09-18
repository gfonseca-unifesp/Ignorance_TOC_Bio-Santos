# =============================================================================
# 08y_magnitude_by_stratum.R — depth zone and latitude structure the MAGNITUDE of error (08w, 08x): what to do?
# =============================================================================
# A magnitude-only structure means the model errs more in some strata without a systematic offset, so changing how
# depth or latitude enter the mean model is not expected to raise R2. Three checks decide the action:
#   (1) bias by stratum in withheld regions: is there really no systematic offset? (mean error per level)
#   (2) does the ignorance map already carry the structure? realised vs expected RMSE per level, for C0 and for the
#       lowest-MACE error model, with 90% interval coverage (|e| <= 1.645 expected)
#   (3) how much of it is irreducible? noise floor per level = median within-cell SD of log10 TOC (cells with >= 3 records)
# Withheld points = nested leave-region-out (07_calibration.R) of the clean chain.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
ITS <- c("iter1c_clean", "iter2c_clean")
cells <- read.csv(file.path(DIRS$data, "toc_cells.csv"))
key <- function(lon, lat) paste(round(lon, 3), round(lat, 3))
zone <- function(z) as.character(cut(z, c(-Inf, 200, 1000, 3000, Inf), labels = c("shelf <=200 m", "slope 200-1000 m", "1000-3000 m", "abyss >3000 m")))
band <- function(lat) as.character(cut(abs(lat), c(0, 30, 60, 90), include.lowest = TRUE, labels = c("tropics 0-30", "mid-latitudes 30-60", "high latitudes 60-90")))
out <- list()
for (it in ITS) {
  o <- read.csv(file.path(ROOT, "outputs", it, "tables", "calibration_outer_points_toc.csv"))
  tr <- read.csv(file.path(ROOT, "outputs", it, "data", "toc_training_table.csv"))
  o$depth <- tr$depth[match(key(o$lon, o$lat), key(tr$x_lon, tr$y_lat))]
  best <- read.csv(file.path(ROOT, "outputs", it, "tables", "toc_summary.csv")); best <- best$value[best$metric == "calibration_best_method"]
  o$e <- o$pred - o$obs; o$exp_best <- o[[paste0("exp_", best)]]
  for (pn in c("depth_zone", "latitude_band")) {
    g <- if (pn == "depth_zone") zone(o$depth) else band(o$lat)
    for (lv in sort(unique(na.omit(g)))) { j <- which(g == lv)
      out[[length(out) + 1]] <- data.frame(iteration = it, partition = pn, level = lv, n = length(j), mean_error = mean(o$e[j]),
        realised_rmse = sqrt(mean(o$e[j]^2)), expected_rmse_C0 = sqrt(mean(o$exp_C0[j]^2)), ratio_C0 = sqrt(mean(o$e[j]^2)) / sqrt(mean(o$exp_C0[j]^2)),
        coverage90_C0 = mean(abs(o$e[j]) <= 1.645 * o$exp_C0[j]), best_method = best,
        ratio_best = sqrt(mean(o$e[j]^2)) / sqrt(mean(o$exp_best[j]^2)), coverage90_best = mean(abs(o$e[j]) <= 1.645 * o$exp_best[j]))
    }
  }
}
R <- do.call(rbind, out)
cg <- cells[!is.na(cells$y_sd), ]
trc <- read.csv(file.path(ROOT, "outputs", "iter1c_clean", "data", "toc_training_table.csv"))
zc <- trc$depth[match(key(cg$x_lon, cg$y_lat), key(trc$x_lon, trc$y_lat))]
nf <- do.call(rbind, lapply(c("depth_zone", "latitude_band"), function(pn) {
  g <- if (pn == "depth_zone") zone(zc) else band(cg$y_lat); ok <- !is.na(g)
  a <- aggregate(cg$y_sd[ok], list(level = g[ok]), function(x) c(n = length(x), median = median(x)))
  data.frame(partition = pn, level = a$level, n_cells_with_sd = a$x[, 1], noise_floor = a$x[, 2]) }))
R <- merge(R, nf, by = c("partition", "level"), all.x = TRUE)
R$reducible_pct <- 100 * (1 - R$noise_floor^2 / R$realised_rmse^2)
R <- R[order(R$iteration, R$partition, R$level), ]
write.csv(R, file.path(DIRS$compare, "magnitude_by_stratum.csv"), row.names = FALSE)
print(R, digits = 3, row.names = FALSE)
msg("magnitude by stratum done")
