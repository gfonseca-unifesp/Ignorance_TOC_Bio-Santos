# =============================================================================
# 06_figures.R — Demo 3b figures: TOC prediction, ignorance map, DI vs error, calibration
# =============================================================================
#  (a) predicted surface-sediment TOC (%), observation cells as dots
#  (b) ignorance map: expected RMSE from the DI-error profile inside the AOA;
#      grey = outside the AOA, where error cannot be estimated (extrapolation)
#  (c) error vs dissimilarity: spatial CV and withheld regions, with the noise floor
#  (d) calibration of the ignorance map in withheld regions (expected vs realised RMSE)

source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(tidyterra); library(patchwork); library(rnaturalearth) })

PROJ <- "+proj=eqearth +datum=WGS84 +units=m"
rr   <- function(stem) rast(file.path(DIRS$rasters, sprintf("toc_%s.tif", stem)))
tb   <- function(f) read.csv(file.path(DIRS$tables, f))
summ <- tb("toc_summary.csv"); val <- function(m) { v <- summ$value[summ$metric == m]; if (length(v)) v else "NA" }

tp   <- project(rr("prediction_pct"), PROJ, res = 25000)
toc  <- project(rr("prediction_pct"), tp, method = "bilinear")
ign  <- project(rr("ignorance_expectedRMSE"), tp, method = "bilinear")
aoa  <- project(rr("AOA"), tp, method = "near")
out_poly <- st_as_sf(as.polygons(ifel(aoa == 0, 1, NA), dissolve = TRUE))
obs  <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
obs_sf <- st_transform(st_as_sf(obs, coords = c("x_lon", "y_lat"), crs = 4326), PROJ)

edge  <- rbind(cbind(-180, seq(-90, 90, 1)), cbind(seq(-180, 180, 1), 90), cbind(180, seq(90, -90, -1)), cbind(seq(180, -180, -1), -90))
globe <- st_sfc(st_polygon(list(sf::sf_project(from = "EPSG:4326", to = PROJ, pts = edge))), crs = PROJ)
land  <- st_transform(ne_countries(scale = 50, returnclass = "sf"), PROJ)
map_layers <- function() list(
  geom_sf(data = land, fill = "grey80", colour = NA, inherit.aes = FALSE),
  geom_sf(data = globe, fill = NA, colour = "grey55", linewidth = 0.2, inherit.aes = FALSE),
  coord_sf(crs = PROJ, expand = FALSE, datum = NA), theme_void(base_size = 8),
  theme(legend.position = "bottom", legend.key.width = unit(14, "mm"), legend.key.height = unit(2.5, "mm"),
        plot.title = element_text(face = "bold", size = 9)))

pa <- ggplot() + geom_spatraster(data = toc) +
  scale_fill_viridis_c(option = "cividis", trans = "log10", limits = c(0.05, 5), oob = scales::squish, na.value = NA,
                       name = "predicted TOC (%)", breaks = c(0.05, 0.1, 0.5, 1, 5), labels = c("<=0.05", "0.1", "0.5", "1", ">=5")) +
  geom_sf(data = obs_sf, size = 0.05, colour = "#FF4D4D", alpha = 0.5, inherit.aes = FALSE) +
  map_layers() + ggtitle(sprintf("a  Global prediction of surface-sediment TOC (random forest; red dots: %s observation cells)",
                                 format(nrow(obs), big.mark = ",")))

ign_lim <- as.numeric(global(ign, function(z) quantile(z, c(0.02, 0.98), na.rm = TRUE))[1, ])
pb <- ggplot() + geom_spatraster(data = ign) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = ign_lim, oob = scales::squish, na.value = NA,
                       name = "expected RMSE (log10 TOC; 2nd-98th percentile stretch)") +
  geom_sf(data = out_poly, fill = "grey45", colour = NA, alpha = 0.85, inherit.aes = FALSE) +
  map_layers() + ggtitle(sprintf("b  Ignorance map: expected error inside the AOA; grey = outside the AOA (%s%% of the ocean), error not estimable",
                                 val("pct_ocean_outside_AOA")))

nf <- as.numeric(val("noise_floor_withincell_sd_log"))
cvb <- tb("DI_vs_CVerror_bins_toc.csv"); cvb$source <- "spatial CV"
trb <- tb("transfer_DI_bins_toc.csv");   trb$source <- "withheld regions"
eb  <- rbind(cvb[, c("DInorm_median", "RMSE", "RMSE_lo", "RMSE_hi", "source")], trb[, c("DInorm_median", "RMSE", "RMSE_lo", "RMSE_hi", "source")])
pc <- ggplot(eb, aes(DInorm_median, RMSE, colour = source)) +
  geom_hline(yintercept = nf, linetype = 3, colour = "grey40") +
  annotate("text", x = min(eb$DInorm_median), y = nf, label = "noise floor (within-cell SD)", vjust = -0.5, hjust = 0, size = 2.3, colour = "grey30") +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_line() + geom_pointrange(aes(ymin = RMSE_lo, ymax = RMSE_hi), size = 0.25) +
  scale_x_sqrt() + scale_colour_manual(values = c("spatial CV" = "#3E5C76", "withheld regions" = "#D7301F"), name = NULL) +
  labs(x = "dissimilarity (DI / AOA threshold; dashed = AOA edge)", y = "RMSE (log10 TOC), 95% bootstrap CI",
       title = "c  Error grows with dissimilarity",
       subtitle = sprintf("withheld regions, RMSE [95%% CI]\ninside AOA %s | outside %s", val("transfer_RMSE_inside_AOA"), val("transfer_RMSE_outside_AOA"))) +
  theme_classic(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 9))

cal <- tb("transfer_ignorance_calibration_toc.csv")
lim <- range(c(cal$expected_RMSE, cal$lo, cal$hi))
pd <- ggplot(cal, aes(expected_RMSE, realised_RMSE)) + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), colour = "#D7301F", size = 0.25) + geom_line(colour = "#D7301F") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "expected RMSE (ignorance map, fitted without the region)", y = "realised RMSE in the withheld region",
       title = "d  Does the ignorance map predict the error it meets?",
       subtitle = sprintf("10 bins of expected RMSE, withheld regions\nSpearman rho = %s | slope = %s (1 = calibrated)", val("calibration_spearman_expected_realised_rho"), val("calibration_slope_realised_on_expected"))) +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9))

stats_l <- paste0("Observation cells: ", val("n_observation_cells"), " (0.1 deg)\n",
                  "R2 (log10) spatial CV ", val("R2_log_spatialCV"), " | random CV ", val("R2_log_randomCV"), "\n",
                  "RMSE spatial CV ", val("RMSE_log_spatialCV"), " | noise floor ", val("noise_floor_withincell_sd_log"), " (log10 units)")
stats_r <- paste0("Ocean outside AOA: ", val("pct_ocean_outside_AOA"), "% (shelf ", val("pct_outside_AOA_shelf"), "%, slope ",
                  val("pct_outside_AOA_slope"), "%, abyss ", val("pct_outside_AOA_abyss"), "%)\n",
                  "Median distance to nearest observation: CV ", val("median_dist_CV_km"), " km | global map ", val("median_dist_prediction_km"), " km\n",
                  "Predictors (ffs): ", gsub(";", ", ", val("predictors_selected")))
pbox <- ggplot() + annotate("text", x = 0, y = 1, label = stats_l, hjust = 0, vjust = 1, size = 2.5, lineheight = 1.15) +
  annotate("text", x = 0.42, y = 1, label = stats_r, hjust = 0, vjust = 1, size = 2.5, lineheight = 1.15) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

main <- (pa / pb / (pc | pd) / pbox) + plot_layout(heights = c(1.15, 1.15, 0.95, 0.16))
ggsave(file.path(DIRS$figs, "toc_main.png"), main, width = 250, height = 300, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(DIRS$figs, "toc_main.pdf"), main, width = 250, height = 300, units = "mm", bg = "white")
msg("figures written")
