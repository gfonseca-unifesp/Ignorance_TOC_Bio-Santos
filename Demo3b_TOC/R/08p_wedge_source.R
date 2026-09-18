# =============================================================================
# 08p_wedge_source.R — which predictor layer draws the SE Pacific wedge in the TOC predictions?
# =============================================================================
# 08o showed a sharp-edged wedge of low TOC in the SE Pacific that is already present in iteration 1 (carbon
# supply, no texture) and is deepened by texture in iteration 2. Map the preliminary model (no supply) and the
# candidate layers over the same box to find the layer with the same geometry.
Sys.setenv(DEMO3B_ITER = "iter2b_supply_texture")
source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(tidyterra); library(patchwork) })

e   <- ext(-140, -90, -50, 0)
stk <- rast(path_toc_stack())
pre <- rast(file.path(ROOT, "outputs", "iter1_baseline", "rasters", "toc_prediction_log10.tif"))
# rank stretch: a random forest splits on the ordering of values, so even differences in the 7th decimal can carry a
# geometry the model uses; colour = percentile rank of the value within the box
pm <- function(r, title) { rc <- crop(r, e); v <- values(rc)[, 1]; ok <- !is.na(v); v[ok] <- rank(v[ok], ties.method = "average") / sum(ok); values(rc) <- v
  ggplot() + geom_spatraster(data = rc) + scale_fill_viridis_c(limits = c(0, 1), na.value = "grey85") + coord_sf(expand = FALSE) +
    labs(title = title, fill = NULL, x = NULL, y = NULL) + theme_minimal(base_size = 7) +
    theme(legend.key.width = unit(2, "mm"), axis.text = element_blank(), plot.title = element_text(face = "bold", size = 7.5)) }
lay <- c(sed_thick = "sediment thickness", poc_surf = "surface POC", poc_flux = "POC flux at the seafloor", sws_bot = "bottom current speed",
         o2b_mean = "bottom oxygen", tri = "terrain ruggedness (TRI)", depth = "depth", sbt_mean = "bottom temperature", sbs_mean = "bottom salinity",
         chl_mean = "surface chlorophyll", sst_mean = "surface temperature", sss_mean = "surface salinity", phyc_bot = "phytoplankton near bottom",
         dist_coast_km = "distance to coast")
ee  <- ext(xmin(e) - 1, xmax(e) + 1, ymin(e) - 1, ymax(e) + 1)
triF <- focal(crop(stk[["tri"]], ee), w = 11, fun = "mean", na.rm = TRUE)            # ~1-degree mean roughness: survey-dependent texture of the bathymetry
depR <- focal(abs(crop(stk[["depth"]], ee) - focal(crop(stk[["depth"]], ee), w = 3, fun = "mean", na.rm = TRUE)), w = 11, fun = "mean", na.rm = TRUE)
fig <- wrap_plots(c(list(pm(pre, "log10 TOC, preliminary model (no supply)"), pm(triF, "TRI, 1-degree mean"), pm(depR, "depth micro-relief, 1-degree mean")),
                    lapply(names(lay), function(k) pm(stk[[k]], lay[[k]]))), ncol = 6)
ggsave(file.path(DIRS$compare, "fig_prediction_artefact_sources.png"), fig, width = 300, height = 190, units = "mm", dpi = 200, bg = "white")

# transects across the wedge: at the largest steps of the preliminary prediction, the step of each of its predictors,
# normalised by that predictor's median non-zero cell-to-cell change along the transect (values >> 1 = sharp edge)
stk0 <- rast(path_base_toc_stack())
m0 <- readRDS(file.path(ROOT, "outputs", "iter1_baseline", "models", "rf_spatialCV_toc.rds"))
v0 <- setdiff(names(m0$trainingData), c(".outcome", ".weights"))
tr <- do.call(rbind, lapply(c(-30, -35, -40), function(lat) {
  xy <- cbind(seq(-140, -95, by = 0.1), lat); p <- extract(pre, xy)[, 1]; X <- as.matrix(extract(stk0[[v0]], xy))
  dp <- abs(diff(p)); top <- order(dp, decreasing = TRUE)[1:4]
  do.call(rbind, lapply(top, function(t) data.frame(lat = lat, lon = xy[t, 1], prediction_step = dp[t],
    t(sapply(v0, function(v) { d <- abs(diff(X[, v])); md <- median(d[d > 0], na.rm = TRUE); d[t] / md })))))
}))
write.csv(tr, file.path(DIRS$compare, "prediction_artefact_wedge_transects.csv"), row.names = FALSE)
print(tr, digits = 3, row.names = FALSE)
msg("wedge source figure written")
