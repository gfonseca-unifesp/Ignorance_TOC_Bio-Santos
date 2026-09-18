# =============================================================================
# 08f_localization_by_scale.R — at what spatial scale does the ignorance map locate error?
# =============================================================================
# Nested leave-region-out points (07_calibration.R). For each iteration and error model (C0-C3):
# Spearman between expected and realised error at point level, and between expected and realised
# RMSE aggregated in grid cells of 1, 2, 5 and 10 deg (within withheld regions; cells with >= 10
# points) and per region (8). Bootstrap 95% CI over units (500 resamples).
Sys.setenv(DEMO3B_ITER = "iter2_texture")
source("R/00_config.R")
set.seed(CFG$seed)
rd <- function(it) read.csv(file.path(ROOT, "outputs", it, "tables", "calibration_outer_points_toc.csv"))
boot_rho <- function(a, b, B = 500) { n <- length(a); if (n < 5) return(c(NA, NA))
  r <- replicate(B, { i <- sample(n, n, replace = TRUE); suppressWarnings(cor(a[i], b[i], method = "spearman")) })
  quantile(r, c(0.025, 0.975), na.rm = TRUE) }
res <- list()
its_done <- names(ITERS)[file.exists(file.path(ROOT, "outputs", names(ITERS), "tables", "calibration_outer_points_toc.csv"))]
for (it in its_done) {
  p <- rd(it); p$e2 <- (p$pred - p$obs)^2
  for (m in sub("^exp_", "", grep("^exp_C[0-9]$", names(p), value = TRUE))) {
    ex <- p[[paste0("exp_", m)]]
    add <- function(scale, real, expd) { ci <- boot_rho(real, expd)
      res[[length(res) + 1]] <<- data.frame(iteration = it, method = m, scale = scale, n_units = length(real),
        spearman = cor(real, expd, method = "spearman"), ci_lo = ci[1], ci_hi = ci[2]) }
    add("point", sqrt(p$e2), ex)
    for (s in c(1, 2, 5, 10)) {
      g <- paste(p$region, floor(p$lon / s), floor(p$lat / s)); k <- table(g); keep <- names(k)[k >= 10]
      add(sprintf("%g deg", s), sqrt(tapply(p$e2, g, mean))[keep], sqrt(tapply(ex^2, g, mean))[keep])
    }
    add("region (8)", sqrt(tapply(p$e2, p$region, mean)), sqrt(tapply(ex^2, p$region, mean)))
  }
}
res <- do.call(rbind, res); rownames(res) <- NULL
write.csv(res, file.path(DIRS$compare, "localization_by_scale.csv"), row.names = FALSE)
res$txt <- sprintf("%5.2f [%5.2f,%5.2f] n=%d", res$spearman, res$ci_lo, res$ci_hi, res$n_units)
for (it in unique(res$iteration)) { cat("\n==", it, "==\n")
  w <- reshape(res[res$iteration == it, c("method", "scale", "txt")], idvar = "scale", timevar = "method", direction = "wide")
  print(w, row.names = FALSE) }
