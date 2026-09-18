# =============================================================================
# 01_data.R — TOC observations: clean, snap to the ocean grid, aggregate per cell
# =============================================================================
# Source: NN-TOC label compilation (Seiter et al. 2004, Romankevich et al. 2009,
# van der Voort et al. 2021 / MOSAIC, Paradis et al. 2023, regional North Sea and
# Gulf of Mexico sets), extracted from Zenodo 10.5281/zenodo.11186224 (CC-BY 4.0):
#   data/raw/labels/toc_continentalshelves.csv, data/raw/labels/toc_deep.csv
# QC (2026-09-14): units consistent with Seiter (94.8% identical values at shared
# coordinates). The low shelf median reflects North Sea sands; the high "deep" median
# reflects upwelling margins on the slope. This is sampling composition, not units.
#
# One observation per 0.1-deg cell: mean of log10(TOC + offset). The within-cell SD
# (cells with >= 3 records) estimates small-scale variability, i.e. the noise floor
# no model at this resolution can remove.

source("R/00_config.R")
set.seed(CFG$seed)

rd <- function(f, mask) {
  d <- read.csv(file.path(DIRS$toc, f), check.names = FALSE)
  data.frame(mask = mask, lat = d$Latitude, lon = d$Longitude, toc = d[["TOC [%]"]])
}
x <- rbind(rd("toc_continentalshelves.csv", "shelf"), rd("toc_deep.csv", "deep"))
log_steps <- list(); step <- function(d, s) { log_steps[[length(log_steps) + 1]] <<- data.frame(step = s, n = nrow(d)); d }
x <- step(x, "raw records")
x <- x[is.finite(x$lat) & is.finite(x$lon) & abs(x$lat) <= 90 & abs(x$lon) <= 180, ] %>% step("valid coordinates")
x <- x[is.finite(x$toc) & x$toc >= 0 & x$toc <= CFG$toc_max, ] %>% step(sprintf("0 <= TOC <= %g%%", CFG$toc_max))
x <- x[!duplicated(x[, c("lat", "lon", "toc")]), ] %>% step("exact duplicates removed")

# Exclude the North Sea regional set (NN-TOC: "Wenyen Zhang, personal communication,
# 2023, HEREON"). It is not a citable source, and its records are a gridded product:
# a regular 0.0298-deg grid, 98.5% of coordinates with 5 decimals, 89% of TOC values
# with >= 5 decimals, and a within-cell SD of 0.034 log10 units vs 0.12 elsewhere.
# Using it would train on another model's output and give the North Sea 40% of all cells.
# Rule: inside the North Sea box, drop TOC values with >= 4 decimals. That keeps
# 1,743 conventionally reported samples (Seiter, MOSAIC). High-precision values
# elsewhere (e.g. Gulf of Mexico) have irregular coordinates and are kept.
ndec <- function(v, maxd = 6) vapply(v, function(z) { for (k in 0:maxd) if (abs(z * 10^k - round(z * 10^k)) < 1e-7) return(k); maxd + 1 }, numeric(1))
ns_box <- x$lat >= 50 & x$lat <= 62 & x$lon >= -5 & x$lon <= 13
x <- x[!(ns_box & ndec(x$toc) >= 4), ] %>% step("North Sea gridded product excluded")

stk   <- rast(path_base_stack())
ocean <- values(!is.na(stk[["sst_mean"]]))[, 1]
tmpl  <- stk[["sst_mean"]]
cells <- cellFromXY(tmpl, cbind(x$lon, x$lat))
bad   <- which(!is.na(cells) & !(ocean[cells] %in% TRUE))
if (length(bad)) {                                   # snap to nearest ocean cell among 8 neighbours
  adj <- adjacent(tmpl, cells[bad], directions = "queen")
  cells[bad] <- vapply(seq_along(bad), function(j) {
    nb <- adj[j, ]; nb <- nb[!is.na(nb) & ocean[nb] %in% TRUE]
    if (!length(nb)) return(NA_real_)
    xy <- xyFromCell(tmpl, nb); nb[which.min((xy[, 1] - x$lon[bad[j]])^2 + (xy[, 2] - x$lat[bad[j]])^2)]
  }, numeric(1))
}
x$cell <- cells
x <- x[!is.na(x$cell) & ocean[x$cell] %in% TRUE, ] %>% step("in ocean (snapped <= 1 cell)")

x$y <- log10(x$toc + CFG$toc_offset)
# Noise floor: the same sample is often reported by several compilations, so the SD is
# taken across distinct sites (coordinates rounded to 0.01 deg, site mean first),
# for cells with >= 3 distinct sites.
site_sd <- x %>% mutate(site = paste(round(lat, 2), round(lon, 2))) %>%
  group_by(cell, site) %>% summarise(ys = mean(y), .groups = "drop") %>%
  group_by(cell) %>% summarise(n_sites = n(), y_sd = if (n() >= 3) sd(ys) else NA_real_, .groups = "drop")
cellobs <- x %>% group_by(cell) %>% summarise(
  n_records = n(), y = mean(y),
  toc_median = median(toc), mask = names(sort(table(mask), decreasing = TRUE))[1], .groups = "drop") %>%
  left_join(site_sd, by = "cell")
xy <- xyFromCell(tmpl, cellobs$cell); cellobs$x_lon <- xy[, 1]; cellobs$y_lat <- xy[, 2]
log_steps[[length(log_steps) + 1]] <- data.frame(step = "aggregated to 0.1-deg cells", n = nrow(cellobs))

write.csv(cellobs, file.path(DIRS$data, "toc_cells.csv"), row.names = FALSE)
clog <- bind_rows(log_steps)
write.csv(clog, file.path(DIRS$tables, "toc_cleaning_log.csv"), row.names = FALSE)
print(clog)
nf <- cellobs$y_sd[!is.na(cellobs$y_sd)]
msg("cells: %d | records per cell median %g, max %d | within-cell SD of log10 TOC across distinct sites (cells with >= 3 sites, n = %d): median %.3f",
    nrow(cellobs), median(cellobs$n_records), max(cellobs$n_records), length(nf), median(nf))
