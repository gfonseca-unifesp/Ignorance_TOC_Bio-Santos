# =============================================================================
# 08_compare_iterations.R — knowledge gain between analytic iterations
# =============================================================================
# Reads each iteration's outputs (outputs/<iteration>/tables, rasters) and documents
# what changed: predictive skill, error vs dissimilarity, calibration of the ignorance
# map, per-region errors, and the map of change in expected error.

source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(tidyterra); library(rnaturalearth) })
dir.create(DIRS$compare, recursive = TRUE, showWarnings = FALSE)
iters <- names(ITERS)
tabp <- function(it, f) { p <- file.path(ROOT, "outputs", it, "tables", f); if (file.exists(p)) read.csv(p) else NULL }
rasp <- function(it, f) { p <- file.path(ROOT, "outputs", it, "rasters", f); if (file.exists(p)) rast(p) else NULL }
getv <- function(s, m) { v <- s$value[s$metric == m]; if (length(v)) v else NA }
num  <- function(x) suppressWarnings(as.numeric(sub(" .*", "", x)))

rows <- lapply(iters, function(it) {
  s <- tabp(it, "toc_summary.csv"); if (is.null(s)) return(NULL)
  cm <- tabp(it, "calibration_metrics_toc.csv"); best <- getv(s, "calibration_best_method")
  cb <- cm[cm$method == best & cm$subset == "all withheld points", ]
  data.frame(iteration = it, label = ITERS[[it]]$label,
             n_candidates = length(ITERS[[it]]$candidates), predictors_selected = getv(s, "predictors_selected"),
             R2_spatialCV = num(getv(s, "R2_log_spatialCV")), R2_randomCV = num(getv(s, "R2_log_randomCV")),
             RMSE_spatialCV = num(getv(s, "RMSE_log_spatialCV")), noise_floor = num(getv(s, "noise_floor_withincell_sd_log")),
             pct_ocean_outside_AOA = num(getv(s, "pct_ocean_outside_AOA")),
             rho_DI_error_CV = num(getv(s, "spearman_DI_abs_error_CV_rho")),
             transfer_RMSE_inside_AOA = num(getv(s, "transfer_RMSE_inside_AOA")), transfer_RMSE_outside_AOA = num(getv(s, "transfer_RMSE_outside_AOA")),
             calibration_method = best, calib_slope = cb$slope, calib_ratio = cb$ratio_realised_expected,
             calib_coverage90 = cb$coverage90, calib_rho_points = cb$spearman_points, calib_rho_bins = cb$spearman_bins, calib_MACE = cb$MACE)
})
cmp <- do.call(rbind, rows)
cmp$reducible_error_variance_pct <- round(100 * (1 - cmp$noise_floor^2 / cmp$RMSE_spatialCV^2), 1)
write.csv(cmp, file.path(DIRS$compare, "iteration_comparison.csv"), row.names = FALSE)
print(t(cmp))
# figures and the delta map compare one pair (default: the presentation sequence iter1b -> iter2b)
iters <- intersect(c(Sys.getenv("CMP_FROM", "iter1b_supply"), Sys.getenv("CMP_TO", "iter2b_supply_texture")), iters)

reg <- do.call(rbind, lapply(iters, function(it) { r <- tabp(it, "calibration_per_region_toc.csv"); tr <- tabp(it, "transfer_regions_toc.csv")
  if (is.null(r) || is.null(tr)) return(NULL); r$iteration <- it; merge(r, tr[, c("region", "lon", "lat", "R2")], by = "region") }))
write.csv(reg, file.path(DIRS$compare, "iteration_comparison_regions.csv"), row.names = FALSE)

if (length(unique(reg$iteration)) >= 2) {
  reg$region_lab <- sprintf("R%d (%.0f, %.0f)", reg$region, reg$lon, reg$lat)
  pa <- ggplot(reg, aes(realised, reorder(region_lab, realised))) +
    geom_line(aes(group = region_lab), colour = "grey70") + geom_point(aes(colour = iteration), size = 2) +
    geom_point(aes(x = expected_best, colour = iteration), shape = 4, size = 2) +
    labs(x = "RMSE in withheld region (dot = realised, x = expected by ignorance map)", y = NULL, colour = NULL,
         title = "a  Withheld regions: realised and expected error") +
    theme_classic(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  cal <- do.call(rbind, lapply(iters, function(it) { b <- tabp(it, "calibration_bins_toc.csv"); s <- tabp(it, "toc_summary.csv")
    if (is.null(b)) return(NULL); b <- b[b$method == getv(s, "calibration_best_method"), ]; b$iteration <- it; b }))
  lim <- range(c(cal$expected, cal$lo, cal$hi))
  pb <- ggplot(cal, aes(expected, realised, colour = iteration)) + geom_abline(linetype = 2, colour = "grey50") +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.2) + geom_line() + coord_equal(xlim = lim, ylim = lim) +
    labs(x = "expected RMSE", y = "realised RMSE", colour = NULL, title = "b  Calibration of the ignorance map") +
    theme_classic(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  eb <- do.call(rbind, lapply(iters, function(it) { a <- tabp(it, "DI_vs_CVerror_bins_toc.csv"); b <- tabp(it, "transfer_DI_bins_toc.csv")
    if (is.null(a)) return(NULL); rbind(data.frame(iteration = it, test = "spatial CV", x = a$DInorm_median, y = a$RMSE),
                                        data.frame(iteration = it, test = "withheld regions", x = b$DInorm_median, y = b$RMSE)) }))
  pc <- ggplot(eb, aes(x, y, colour = iteration, linetype = test)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
    geom_hline(yintercept = cmp$noise_floor[1], linetype = 3, colour = "grey40") + geom_line() + geom_point(size = 0.8) + scale_x_sqrt() +
    labs(x = "DI / AOA threshold", y = "RMSE (log10 TOC)", colour = NULL, linetype = NULL, title = "c  Error vs dissimilarity (dotted = noise floor)") +
    theme_classic(base_size = 8) + theme(legend.position = "bottom", legend.box = "vertical", plot.title = element_text(face = "bold"))
  r1 <- rasp(iters[1], "toc_ignorance_calibrated_expectedRMSE.tif"); r2 <- rasp(iters[length(iters)], "toc_ignorance_calibrated_expectedRMSE.tif")
  pd <- NULL
  if (!is.null(r1) && !is.null(r2)) {
    PROJ <- "+proj=eqearth +datum=WGS84 +units=m"
    dr <- r2 - r1; names(dr) <- "delta"
    writeRaster(dr, file.path(DIRS$compare, "delta_expected_RMSE_iter2_minus_iter1.tif"), overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
    tp <- project(dr, PROJ, res = 25000); dp <- project(dr, tp, method = "bilinear")
    q <- max(abs(as.numeric(global(dp, function(z) quantile(z, c(0.02, 0.98), na.rm = TRUE))[1, ])))
    land <- st_transform(ne_countries(scale = 50, returnclass = "sf"), PROJ)
    pd <- ggplot() + geom_spatraster(data = dp) +
      scale_fill_distiller(palette = "RdBu", limits = c(-q, q), oob = scales::squish, na.value = NA, name = "change in expected RMSE (iteration 2 - iteration 1)") +
      geom_sf(data = land, fill = "grey80", colour = NA, inherit.aes = FALSE) + coord_sf(crs = PROJ, expand = FALSE, datum = NA) +
      theme_void(base_size = 8) + theme(legend.position = "bottom", legend.key.width = unit(14, "mm"), plot.title = element_text(face = "bold")) +
      ggtitle("d  Where integrating sediment texture changed the ignorance map (blue = less expected error)")
  }
  sub <- sprintf("R2 spatial CV %.2f -> %.2f | calibration slope %.2f -> %.2f | ranking rho (points) %.2f -> %.2f | coverage90 %.2f -> %.2f",
                 cmp$R2_spatialCV[1], cmp$R2_spatialCV[nrow(cmp)], cmp$calib_slope[1], cmp$calib_slope[nrow(cmp)],
                 cmp$calib_rho_points[1], cmp$calib_rho_points[nrow(cmp)], cmp$calib_coverage90[1], cmp$calib_coverage90[nrow(cmp)])
  fig <- ((pa | pb | pc) / pd) + plot_layout(heights = c(1, 1.2)) +
    plot_annotation(title = "Knowledge gain from iteration 1 (baseline) to iteration 2 (+ sediment texture)", subtitle = sub)
  ggsave(file.path(DIRS$compare, "iteration_comparison.png"), fig, width = 270, height = 230, units = "mm", dpi = 300, bg = "white")
}
msg("comparison written to %s", DIRS$compare)
