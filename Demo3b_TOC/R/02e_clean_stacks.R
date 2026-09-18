# =============================================================================
# 02e_clean_stacks.R — predictor stacks of the clean chain (iter1c_clean, iter2c_clean)
# =============================================================================
# Two corrections found by the provenance check (08q-08u; ANALYTIC_LOG), applied to the supply and supply+texture
# stacks while keeping the layer names, so that every later script runs unchanged:
#   dist_coast_km  -> distance to land masses >= 25,000 km2 (0.2 deg, geodesic, bilinear to 0.1 deg, log10 as in 02),
#                     because distance to any island drew straight-edged artefacts in the open ocean (08s, 08t)
#   sbt_mean, sbs_mean, o2b_mean, phyc_bot, sws_bot -> 0.5-deg moving mean where depth > 1500 m, original elsewhere,
#                     because the Bio-ORACLE bottom layers carry seams in deep water (08p, 08q)
Sys.setenv(DEMO3B_ITER = "iter2c_clean")
source("R/00_config.R")

base  <- rast(path_base_toc_stack()); ocean <- !is.na(base[["sst_mean"]])
BOT <- c("sbt_mean", "sbs_mean", "o2b_mean", "phyc_bot", "sws_bot")

f_dl <- file.path(CACHE, "predictor_dist_large_land_0.1deg.tif")          # also written by 08t (same code)
if (!file.exists(f_dl)) {
  land2 <- aggregate(ifel(ocean, NA, 1), fact = 2, fun = "max", na.rm = TRUE)
  pt <- patches(land2, directions = 8)
  za <- zonal(cellSize(land2, unit = "km"), pt, fun = "sum")
  big <- ifel(pt %in% za[[1]][za[[2]] >= 25000], 1, NA)
  dl <- mask(resample(distance(big) / 1000, base[["sst_mean"]], method = "bilinear"), ocean, maskvalues = FALSE)
  lo <- max(global(dl, function(z) quantile(z[z > 0], 0.001, na.rm = TRUE))[[1]], 1e-6)
  dl <- log10(dl + lo); names(dl) <- "dist_large_land"
  writeRaster(dl, f_dl, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
}
sup <- rast(file.path(CACHE, ITERS[["iter1b_supply"]]$stack))
f_sm <- file.path(CACHE, "predictors_toc_bottom_smoothed_deep_0.1deg.tif")   # also written by 08q (same code)
if (!file.exists(f_sm)) {
  dep <- sup[["depth"]]
  sm <- rast(lapply(BOT, function(k) { r <- sup[[k]]; ifel(dep > 1500, focal(r, w = 5, fun = "mean", na.rm = TRUE), r) }))
  names(sm) <- paste0(BOT, "_sm")
  writeRaster(sm, f_sm, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
}
dl <- rast(f_dl); sm <- rast(f_sm)

clean <- function(stk) {
  keep_na <- is.na(stk[["dist_coast_km"]])
  r <- cover(dl, stk[["dist_coast_km"]]); r <- mask(r, keep_na, maskvalues = TRUE); names(r) <- "dist_coast_km"; stk[["dist_coast_km"]] <- r
  for (k in BOT) { s <- mask(cover(sm[[paste0(k, "_sm")]], stk[[k]]), is.na(stk[[k]]), maskvalues = TRUE); names(s) <- k; stk[[k]] <- s }
  stk
}
write_stack <- function(stk, name) {
  path_out <- file.path(CACHE, name); tmp <- sub("[.]tif$", "_writing.tif", path_out)
  writeRaster(stk, tmp, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  unlink(path_out); stopifnot(file.rename(tmp, path_out))
  chk <- rast(path_out)
  msg("%s: %d layers; changed: dist_coast_km mean |diff| %.3f; sbt_mean mean |diff| %.4f", name, nlyr(chk),
      global(abs(chk[["dist_coast_km"]] - sup[["dist_coast_km"]]), "mean", na.rm = TRUE)[1, 1],
      global(abs(chk[["sbt_mean"]] - sup[["sbt_mean"]]), "mean", na.rm = TRUE)[1, 1])
}
sup_clean <- clean(sup)
stopifnot(identical(names(sup_clean), names(sup)))
write_stack(sup_clean, ITERS[["iter1c_clean"]]$stack)
tex <- rast(path_texture_stack())[[TEXTURE_VARS]]
write_stack(c(rast(file.path(CACHE, ITERS[["iter1c_clean"]]$stack)), tex), ITERS[["iter2c_clean"]]$stack)
msg("clean stacks written")
