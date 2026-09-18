# =============================================================================
# 01b_nntoc_composition.R — what did the published NN-TOC validation measure?
# =============================================================================
# NN-TOC (Parameswaran et al. 2025, GMD 18:2521) used all 110,149 label records, averaged them per
# 5-arcmin feature vector (21,125 entries) and validated with a random 85:15 split.
# Here, with the same label files:
#   (1) share of 5-arcmin cells that contain records of the North Sea gridded product (01_data.R rule);
#   (2) distance from a random 15% of cells to the nearest remaining cell, with and without the gridded
#       product, to compare with the distance from global map cells to the nearest observation (231 km,
#       03_model.R geodist diagnostic).
# Approximation: 5-arcmin binning stands in for NN-TOC's per-feature-vector averaging (their variance
# filter removed ~10% of cells).
source("R/00_config.R")
suppressPackageStartupMessages(library(FNN))
rd <- function(f) { d <- read.csv(file.path(DIRS$toc, f), check.names = FALSE); data.frame(lat = d$Latitude, lon = d$Longitude, toc = d[["TOC [%]"]]) }
x <- rbind(rd("toc_continentalshelves.csv"), rd("toc_deep.csv"))
x <- x[is.finite(x$lat) & is.finite(x$lon) & is.finite(x$toc), ]
x <- x[!duplicated(x[, c("lat", "lon", "toc")]), ]
ndec <- function(v, maxd = 6) vapply(v, function(z) { for (k in 0:maxd) if (abs(z * 10^k - round(z * 10^k)) < 1e-7) return(k); maxd + 1 }, numeric(1))
x$ns_grid <- (x$lat >= 50 & x$lat <= 62 & x$lon >= -5 & x$lon <= 13) & ndec(x$toc) >= 4
x$cell5 <- paste(floor(x$lon * 12), floor(x$lat * 12))
cells <- aggregate(cbind(n = 1, ns = ns_grid) ~ cell5, x, sum)
xy <- do.call(rbind, strsplit(cells$cell5, " ")); lon <- (as.numeric(xy[, 1]) + 0.5) / 12; lat <- (as.numeric(xy[, 2]) + 0.5) / 12
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))
P <- to_xyz(lon, lat)
split_stats <- function(keep, label) {
  set.seed(CFG$seed); idx <- which(keep); te <- sample(idx, round(0.15 * length(idx))); tr <- setdiff(idx, te)
  dd <- chord_km(get.knnx(P[tr, ], P[te, , drop = FALSE], k = 1)$nn.dist[, 1])
  data.frame(cell_set = label, n_cells = length(idx), median_km = median(dd), p90_km = unname(quantile(dd, 0.9)),
             pct_lt_10km = 100 * mean(dd < 10), pct_lt_50km = 100 * mean(dd < 50))
}
out <- list(
  composition = data.frame(records = nrow(x), records_gridded_product = sum(x$ns_grid), cells_5arcmin = nrow(cells),
                           pct_cells_any_gridded = 100 * mean(cells$ns > 0), pct_cells_only_gridded = 100 * mean(cells$ns == cells$n)),
  random_split_proximity = rbind(split_stats(rep(TRUE, nrow(cells)), "all 5-arcmin cells"),
                                 split_stats(cells$ns == 0, "without gridded product")))
dir.create(file.path(ROOT, "outputs", "comparison"), showWarnings = FALSE, recursive = TRUE)
write.csv(out$composition, file.path(ROOT, "outputs", "comparison", "nntoc_composition.csv"), row.names = FALSE)
write.csv(out$random_split_proximity, file.path(ROOT, "outputs", "comparison", "nntoc_random_split_proximity.csv"), row.names = FALSE)
print(out, digits = 3)
