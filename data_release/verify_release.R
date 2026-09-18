# =============================================================================
# verify_release.R — F4 acceptance: read every file of the release and check its size
# =============================================================================
#   Rscript verify_release.R      (run from this folder)
suppressPackageStartupMessages(library(terra))
ok <- TRUE
chk <- function(label, got, want) {
  pass <- isTRUE(all.equal(got, want))
  cat(sprintf("%-46s %12s  expected %-12s %s\n", label, format(got, big.mark = ","),
              format(want, big.mark = ","), if (pass) "OK" else "FAIL"))
  if (!pass) ok <<- FALSE
  invisible(pass)
}
nrows <- function(f) nrow(read.csv(f))

# the number of ocean cells is defined by the map itself: the non-NA cells of the DI layer
g <- rast("case1_grid_iter2c.tif")
n_ocean <- global(!is.na(g[["DI"]]), "sum", na.rm = TRUE)[1, 1]

chk("case1_toc_training_cells.csv (rows)",        nrows("case1_toc_training_cells.csv"), 12944)
chk("case1_calibration_outer_points.csv (rows)",  nrows("case1_calibration_outer_points.csv"), 12944)
chk("case1_grid_iter2c.tif (layers)",             nlyr(g), 7L)
chk("case1_grid_iter2c.csv.gz (rows)",            nrows(gzfile("case1_grid_iter2c.csv.gz")), n_ocean)
chk("case1_fig_S10_environmental_space (rows)",   nrows(gzfile("case1_fig_S10_environmental_space.csv.gz")), 60000)
chk("case1_sampling_priority_sites.csv (rows)",   nrows("case1_sampling_priority_sites.csv"), 50)
chk("case1_sampling_priority_uncertainty (rows)", nrows("case1_sampling_priority_uncertainty_sites.csv"), 50)
chk("case1_sampling_priority_curve.csv (rows)",   nrows("case1_sampling_priority_curve.csv"), 51)

s <- read.csv("case2_station_ignorance.csv")
chk("case2_station_ignorance.csv (rows)",         nrow(s), 198)
chk("case2_station_ignorance.csv (stations)",     length(unique(s$station)) / 2, 99)
chk("case2_station_ignorance.csv (training)",     sum(s$set == "train"), 160)
chk("case2_station_ignorance.csv (test)",         sum(s$set == "test"), 38)
chk("case2_grid_layers.tif (layers)",             nlyr(rast("case2_grid_layers.tif")), 8L)
chk("case2_sampling_priority_A_stations (rows)",  nrows("case2_sampling_priority_A_stations.csv"), 30)
chk("case2_sampling_priority_B_sites (rows)",     nrows("case2_sampling_priority_B_sites.csv"), 20)
chk("case2_sampling_priority_curve (rows)",       nrows("case2_sampling_priority_curve.csv"), 21)

cat("\n", if (ok) "release verified" else "RELEASE HAS PROBLEMS", "\n", sep = "")
if (!ok) quit(status = 1)
