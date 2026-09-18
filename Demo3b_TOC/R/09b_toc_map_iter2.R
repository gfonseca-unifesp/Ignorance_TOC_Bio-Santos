# =============================================================================
# 09b_toc_map_iter2.R — the global TOC map after iteration 2, what changed, and its ignorance map
# =============================================================================
# (a) predicted TOC (%) of iteration 2 (base model with carbon supply + sediment texture)
# (b) change in predicted log10 TOC from iteration 1 (supply) to iteration 2 (+ texture)
# (c) calibrated expected error of iteration 2 (the ignorance map; C0, best nested calibration)
# Area-weighted summaries of the change by depth zone: outputs/comparison/toc_map_change_iter1_to_iter2.csv
Sys.setenv(DEMO3B_ITER = Sys.getenv("CMP_TO", "iter2b_supply_texture"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(tidyterra); library(rnaturalearth) })

IT1 <- Sys.getenv("CMP_FROM", "iter1b_supply"); IT2 <- Sys.getenv("CMP_TO", "iter2b_supply_texture")
rp  <- function(it, f) rast(file.path(ROOT, "outputs", it, "rasters", f))
y2  <- rp(IT2, "toc_prediction_log10.tif"); y1 <- rp(IT1, "toc_prediction_log10.tif")
pct <- clamp(rp(IT2, "toc_prediction_pct.tif"), lower = 0.01)
ign <- rp(IT2, "toc_ignorance_calibrated_expectedRMSE.tif")
dep <- rast(path_toc_stack())[["depth"]]
invisible(compareGeom(y1, y2)); invisible(compareGeom(y2, ign)); invisible(compareGeom(y2, dep))
dlt <- y2 - y1; names(dlt) <- "change"

# --- area-weighted summary of the change -----------------------------------------------------------------------
v <- data.frame(y1 = values(y1)[, 1], y2 = values(y2)[, 1], e = values(ign)[, 1], z = values(dep)[, 1],
                a = values(cellSize(y2, unit = "km"))[, 1])
v <- v[complete.cases(v[, c("y1", "y2")]), ]
v$d <- v$y2 - v$y1
v$zone <- cut(v$z, c(-Inf, 200, 1000, 3000, Inf), labels = c("shelf (<= 200 m)", "slope (200-1000 m)", "1000-3000 m", "abyss (> 3000 m)"))
wmed <- function(x, w) { o <- order(x); x[o][which(cumsum(w[o]) >= sum(w) / 2)[1]] }
summ_fun <- function(x, lab) { ke <- !is.na(x$e)
  data.frame(zone = lab, pct_ocean_area = 100 * sum(x$a) / sum(v$a),
             median_TOC_pct_iter1 = 10^wmed(x$y1, x$a) - CFG$toc_offset, median_TOC_pct_iter2 = 10^wmed(x$y2, x$a) - CFG$toc_offset,
             mean_change_log10 = weighted.mean(x$d, x$a),
             pct_area_higher_by_gt_0.1 = 100 * sum(x$a[x$d > 0.1]) / sum(x$a), pct_area_lower_by_gt_0.1 = 100 * sum(x$a[x$d < -0.1]) / sum(x$a),
             median_expected_RMSE_iter2 = wmed(x$e[ke], x$a[ke])) }
tab <- rbind(summ_fun(v, "whole ocean"), do.call(rbind, lapply(levels(v$zone), function(z) summ_fun(v[which(v$zone == z), ], z))))
write.csv(tab, file.path(DIRS$compare, "toc_map_change_iter1_to_iter2.csv"), row.names = FALSE)
print(tab, digits = 3, row.names = FALSE)

# --- maps --------------------------------------------------------------------------------------------------------
PROJ <- "+proj=eqearth +datum=WGS84 +units=m"
land <- st_transform(ne_countries(scale = 50, returnclass = "sf"), PROJ)
base_map <- function(r, sc, title, subtitle, legend = "right") ggplot() + geom_spatraster(data = project(r, PROJ, res = 25000, method = "bilinear")) + sc +
  geom_sf(data = land, fill = "grey82", colour = NA, inherit.aes = FALSE) +
  coord_sf(crs = PROJ, expand = FALSE, datum = NA) + theme_void(base_size = 8) + labs(title = title, subtitle = subtitle) +
  theme(plot.title = element_text(face = "bold", size = 9), legend.position = legend)
s2 <- read.csv(file.path(DIRS$tables, "toc_summary.csv")); gv <- function(m) s2$value[s2$metric == m]
op <- read.csv(file.path(DIRS$tables, "calibration_outer_points_toc.csv"))
eq <- as.numeric(quantile(v$e, c(0.01, 0.99), na.rm = TRUE))

pa <- base_map(pct, scale_fill_viridis_c(trans = "log10", limits = c(0.05, 5), oob = scales::squish, breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2.5, 5),
                                         labels = c("0.05", "0.1", "0.25", "0.5", "1", "2.5", "5"), na.value = NA, name = "TOC (%)"),
  "a  Seafloor organic carbon after iteration 2 (carbon supply + sediment texture)",
  sprintf("0.1-degree prediction of surface-sediment TOC; spatial-CV R2 %s; RMSE in withheld regions %.3f log10 units",
          gv("R2_log_spatialCV"), sqrt(mean((op$pred - op$obs)^2))))
pb <- base_map(dlt, scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-0.3, 0.3),
                                         oob = scales::squish, na.value = NA, name = "change\n(log10 TOC)"),
  "b  What integrating sediment texture changed", "iteration 2 minus iteration 1; red = more carbon predicted with texture", legend = "bottom")
pc <- base_map(ign, scale_fill_viridis_c(option = "magma", direction = -1, limits = eq, oob = scales::squish, na.value = NA, name = "expected RMSE\n(log10 TOC)"),
  "c  Ignorance map after iteration 2",
  sprintf("calibrated expected error; 90%% interval coverage %.2f in withheld regions", as.numeric(gv("calibration_best_coverage90"))), legend = "bottom")
fig <- pa / (pb | pc) + plot_layout(heights = c(1.15, 1))
ggsave(file.path(DIRS$compare, "fig_toc_map_iter2.png"), fig, width = 260, height = 230, units = "mm", dpi = 250, bg = "white")
msg("TOC map after iteration 2 written")
