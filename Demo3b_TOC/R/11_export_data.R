# =============================================================================
# 11_export_data.R — F4 of ROADMAP_MSv6: the Case 1 data package (answers comment M4)
# =============================================================================
# Writes to Ignorance_MS/data_release/ everything a reader needs to check the map and to reuse it:
# the training cells with their predictors and folds, the gridded outputs (raster and table), the
# environmental-space sample behind the ignorance-spaces figure, the ranked sampling sites, and the
# realised-versus-expected errors of the nested calibration.
#   DEMO3B_ITER=iter2c_clean Rscript R/11_export_data.R
# Raw NN-TOC labels are NOT redistributed here: only the aggregated cells used by the model (see the
# release README on what to cite and what not to redistribute).
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages(library(FNN))
REL <- file.path(dirname(ROOT), "data_release")
dir.create(REL, recursive = TRUE, showWarnings = FALSE)
msg("data release -> %s", REL)

## ---- 1. training cells ---------------------------------------------------------------
d <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)                                   # the 8 regions of 05 / 07 / 08g / 08h / 08aa
km <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic",
        "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
d$region <- RN[km$cluster]
front <- c("cell", "x_lon", "y_lat", "n_records", "n_sites", "y", "toc_median", "y_sd", "mask")
d <- d[, c(front, setdiff(names(d), c(front, "block", "fold", "region")), "block", "fold", "region")]
write.csv(d, file.path(REL, "case1_toc_training_cells.csv"), row.names = FALSE)
msg("training cells: %d rows, %d columns", nrow(d), ncol(d))

## ---- 2. gridded outputs ----------------------------------------------------------------
summ <- read.csv(file.path(DIRS$tables, "toc_summary.csv"))
thr  <- as.numeric(summ$value[summ$metric == "AOA_DI_threshold"])
r <- c(rast(file.path(DIRS$rasters, "toc_prediction_pct.tif")),
       rast(file.path(DIRS$rasters, "toc_prediction_log10.tif")),
       rast(file.path(DIRS$rasters, "toc_DI.tif")),
       rast(file.path(DIRS$rasters, "toc_DI.tif")) / thr,
       rast(file.path(DIRS$rasters, "toc_AOA.tif")),
       rast(file.path(DIRS$rasters, "toc_ignorance_calibrated_expectedRMSE.tif")),
       rast(file.path(DIRS$rasters, "toc_distance_nearest_observation_km.tif")))
names(r) <- c("toc_pct", "toc_log10", "DI", "DI_over_threshold", "inside_AOA",
              "expected_RMSE_calibrated_log10", "distance_nearest_observation_km")
writeRaster(r, file.path(REL, "case1_grid_iter2c.tif"), overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"))
msg("grid raster: %d layers, %d x %d cells", nlyr(r), nrow(r), ncol(r))

oc <- which(!is.na(values(r[["DI"]])[, 1]))
xy <- xyFromCell(r, oc)
g <- data.frame(cell = oc, lon = round(xy[, 1], 4), lat = round(xy[, 2], 4),
                round(as.data.frame(r[oc]), 4))
# an ocean cell is attached to the nearest region centroid; the regions are clusters of the OBSERVATION
# cells, so a cell far from every observation belongs to its region only nominally (see the README)
g$region <- RN[FNN::get.knnx(km$centers, to_xyz(g$lon, g$lat), k = 1)$nn.index[, 1]]
write.csv(g, gzfile(file.path(REL, "case1_grid_iter2c.csv.gz")), row.names = FALSE)
msg("grid table: %d ocean cells (non-NA DI)", nrow(g))

## ---- 3. the exclusion rule, in words and counts -----------------------------------------
writeLines(c(
  "# Case 1 — the North Sea exclusion rule",
  "",
  "The NN-TOC label compilation (Zenodo 10.5281/zenodo.11186224, CC-BY 4.0) holds 110,149 records of total",
  "organic carbon in surface sediments. One regional set inside the North Sea box is a **gridded product**,",
  "not a set of measurements, and is attributed to a personal communication:",
  "",
  "- coordinates on a regular 0.0298 deg grid; 98.5% of coordinates with 5 decimals;",
  "- 89% of TOC values with 5 or more decimals;",
  "- within-cell standard deviation 0.034 log10 units, against 0.12 elsewhere.",
  "",
  "**Rule applied here:** drop records inside the North Sea box whose TOC value carries 4 or more decimal",
  "places. It removes **84,265 records** and keeps **1,743** conventionally reported samples from that box.",
  "Keeping the product would train the model on another model's output and give the North Sea 40% of all cells.",
  "",
  "After the rule, 25,884 measurements aggregate into **12,944 cells** of 0.1 deg (mean of log10 TOC per cell).",
  "",
  "Source: Demo3b_TOC/R/01_data.R and Demo3b_TOC/README.md."),
  file.path(REL, "case1_toc_excluded_northsea_rule.md"))

## ---- 4. tables that travel with the map ---------------------------------------------------
file.copy(file.path(DIRS$compare, "environmental_space_sample_toc.csv.gz"),
          file.path(REL, "case1_fig_S10_environmental_space.csv.gz"), overwrite = TRUE)
file.copy(file.path(DIRS$compare, "environmental_space_loadings_toc.csv"),
          file.path(REL, "case1_fig_S10_environmental_space_loadings.csv"), overwrite = TRUE)
file.copy(file.path(DIRS$compare, "sampling_priority_sites_toc.csv"),
          file.path(REL, "case1_sampling_priority_sites.csv"), overwrite = TRUE)
file.copy(file.path(DIRS$compare, "sampling_priority_uncertainty_sites_toc.csv"),
          file.path(REL, "case1_sampling_priority_uncertainty_sites.csv"), overwrite = TRUE)
file.copy(file.path(DIRS$compare, "sampling_priority_curve_toc.csv"),
          file.path(REL, "case1_sampling_priority_curve.csv"), overwrite = TRUE)
cal <- read.csv(file.path(DIRS$tables, "calibration_outer_points_toc.csv"))
write.csv(cal, file.path(REL, "case1_calibration_outer_points.csv"), row.names = FALSE)
msg("calibration points: %d rows (withheld-region predictions with expected errors C0-C5)", nrow(cal))
msg("case 1 export done")
