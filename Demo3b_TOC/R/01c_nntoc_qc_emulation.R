# =============================================================================
# 01c_nntoc_qc_emulation.R — emulate the NN-TOC label QC to estimate the share of gridded-product entries
# =============================================================================
# Parameswaran et al. (2025) describe their label QC as: remove exact duplicates (lon, lat, TOC); group
# labels that share a feature vector (5-arcmin grid); exclude groups whose SD exceeds 20% of the maximum;
# average the remaining groups. Separate shelf and deep sets; 21,125 entries in total.
# Their input files hold every raw record, and the reduction happens in their training code, so here the
# QC is emulated on the same label files. Two readings of "20% of the maximum" are tested (within-group
# maximum; overall maximum). The reading whose entry count matches 21,125 is the plausible one.
source("R/00_config.R")
rd <- function(f, dom) { d <- read.csv(file.path(DIRS$toc, f), check.names = FALSE); data.frame(dom = dom, lat = d$Latitude, lon = d$Longitude, toc = d[["TOC [%]"]]) }
x <- rbind(rd("toc_continentalshelves.csv", "shelf"), rd("toc_deep.csv", "deep"))
x <- x[is.finite(x$lat) & is.finite(x$lon) & is.finite(x$toc), ]
x <- x[!duplicated(x[, c("dom", "lon", "lat", "toc")]), ]
ndec <- function(v, maxd = 6) vapply(v, function(z) { for (k in 0:maxd) if (abs(z * 10^k - round(z * 10^k)) < 1e-7) return(k); maxd + 1 }, numeric(1))
x$grid <- (x$lat >= 50 & x$lat <= 62 & x$lon >= -5 & x$lon <= 13) & ndec(x$toc) >= 4
x$g <- paste(x$dom, floor(x$lon * 12), floor(x$lat * 12))
G <- x %>% group_by(dom, g) %>%
  summarise(n = n(), sd = if (n() > 1) sd(toc) else 0, mx = max(toc), grid_all = all(grid), grid_any = any(grid), .groups = "drop")
gmax <- tapply(x$toc, x$dom, max)
variants <- list(`SD <= 20% of within-group max` = G$sd <= 0.2 * G$mx,
                 `SD <= 20% of overall max (per set)` = G$sd <= 0.2 * gmax[G$dom],
                 `no variance filter` = rep(TRUE, nrow(G)))
out <- do.call(rbind, lapply(names(variants), function(v) { k <- variants[[v]]
  do.call(rbind, lapply(c("shelf", "deep", "all"), function(s) { i <- k & (s == "all" | G$dom == s)
    data.frame(qc_variant = v, set = s, entries = sum(i), pct_entries_only_gridded = 100 * mean(G$grid_all[i]),
               pct_entries_any_gridded = 100 * mean(G$grid_any[i]), difference_to_published_21125 = if (s == "all") sum(i) - 21125 else NA)
  }))
}))
dir.create(DIRS$compare, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(DIRS$compare, "nntoc_qc_emulation.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
