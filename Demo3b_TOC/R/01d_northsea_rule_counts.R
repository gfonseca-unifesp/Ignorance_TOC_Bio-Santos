# =============================================================================
# 01d_northsea_rule_counts.R — F8.1: the counts of the North Sea exclusion rule, as a table
# =============================================================================
# 01_data.R applies the rule but logs only the total left after it; the numbers quoted in the text
# (records excluded, conventional samples kept in the box, noise-floor cells) were written in comments.
# This script reproduces the same filters on the raw NN-TOC labels and writes them to a table.
#   Rscript R/01d_northsea_rule_counts.R      -> outputs/comparison/northsea_rule_counts.csv
source("R/00_config.R")
rd <- function(f) { d <- read.csv(file.path(DIRS$toc, f), check.names = FALSE)
  data.frame(lat = d$Latitude, lon = d$Longitude, toc = d[["TOC [%]"]]) }
x <- rbind(rd("toc_continentalshelves.csv"), rd("toc_deep.csv"))
n_raw <- nrow(x)
x <- x[is.finite(x$lat) & is.finite(x$lon) & abs(x$lat) <= 90 & abs(x$lon) <= 180, ]
x <- x[is.finite(x$toc) & x$toc >= 0 & x$toc <= CFG$toc_max, ]
x <- x[!duplicated(x[, c("lat", "lon", "toc")]), ]
ndec <- function(v, maxd = 6) vapply(v, function(z) { for (k in 0:maxd) if (abs(z * 10^k - round(z * 10^k)) < 1e-7) return(k); maxd + 1 }, numeric(1))
ns <- x$lat >= 50 & x$lat <= 62 & x$lon >= -5 & x$lon <= 13          # the box of 01_data.R
excl <- ns & ndec(x$toc) >= 4
cells <- read.csv(file.path(DIRS$data, "toc_cells.csv"))
out <- data.frame(
  quantity = c("raw_records", "records_after_basic_qc", "records_in_north_sea_box", "excluded_gridded_product",
               "kept_in_box_conventional", "records_after_rule", "pct_records_excluded", "cells_0.1deg",
               "noise_floor_cells_3plus_sites"),
  value = c(n_raw, nrow(x), sum(ns), sum(excl), sum(ns & !excl), nrow(x) - sum(excl),
            round(100 * sum(excl) / n_raw, 1), nrow(cells), sum(cells$n_sites >= 3, na.rm = TRUE)))
write.csv(out, file.path(DIRS$compare, "northsea_rule_counts.csv"), row.names = FALSE)
print(out, row.names = FALSE)
