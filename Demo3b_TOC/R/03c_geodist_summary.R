# =============================================================================
# 03c_geodist_summary.R — F8.1: medians of the geographic-distance diagnostic, as a table
# =============================================================================
# 03_model.R writes the point-level CAST::geodist output; the text quotes its medians (a randomly
# withheld cell lies ~7 km from training data; the map predicts cells ~231 km from any observation).
#   Rscript R/03c_geodist_summary.R   -> outputs/comparison/geodist_summary_toc.csv
source("R/00_config.R")
g <- read.csv(file.path(ROOT, "outputs", "iter1_baseline", "tables", "geodist_toc.csv"))
out <- aggregate(dist_km ~ what, g, function(v) c(n = length(v), median = median(v), p90 = unname(quantile(v, 0.9))))
out <- data.frame(what = out$what, out$dist_km)
write.csv(out, file.path(DIRS$compare, "geodist_summary_toc.csv"), row.names = FALSE)
print(out, row.names = FALSE)
