# =============================================================================
# 08z_calibration_selection.R — which error model should build the ignorance map?
# =============================================================================
# 07_calibration.R fits six error models (C0-C5; C4-C5 add depth and |latitude|, following the error-structure
# test in 08w/08x) and evaluates them in nested leave-region-out. Two things can be asked of an ignorance map:
#   (i)  its expected error is right on average within strata  -> MACE_strata (mean |realised - expected| RMSE
#        over 4 depth zones + 3 latitude bands);
#   (ii) its stated interval means what it says within strata  -> coverage_dev (mean |coverage90 - 0.90|).
# A heteroscedastic mean can improve (i) and still fail (ii) if the error distribution is heavy-tailed where the
# expected error is small. This script bootstraps both criteria over the withheld regions (the unit of the test),
# so the choice is made with uncertainty rather than on point estimates.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
set.seed(CFG$seed)
ITS <- strsplit(Sys.getenv("CAL_SEL_ITERS", "iter1c_clean,iter2c_clean"), ",")[[1]]
B <- 1000
zone <- function(z) as.character(cut(z, c(-Inf, 200, 1000, 3000, Inf), labels = c("shelf <=200 m", "slope 200-1000 m", "1000-3000 m", "abyss >3000 m")))
band <- function(lat) as.character(cut(abs(lat), c(0, 30, 60, 90), include.lowest = TRUE, labels = c("tropics 0-30", "mid-latitudes 30-60", "high latitudes 60-90")))

crit <- function(e2, ex, g) {           # g = stratum label per point
  lv <- sort(unique(na.omit(g)))
  cov <- vapply(lv, function(l) { j <- which(g == l); mean(sqrt(e2[j]) <= 1.645 * ex[j]) }, 0)
  mac <- vapply(lv, function(l) { j <- which(g == l); abs(sqrt(mean(e2[j])) - sqrt(mean(ex[j]^2))) }, 0)
  # rescaled interval: k makes the POOLED coverage exactly 0.90, so what is left is how evenly the interval holds
  # across strata (a heteroscedastic mean should flatten this even when the global scale is off)
  k <- unname(quantile(sqrt(e2) / pmax(ex, 1e-6), 0.90))
  cvk <- vapply(lv, function(l) { j <- which(g == l); mean(sqrt(e2[j]) <= k * ex[j]) }, 0)
  c(MACE_strata = mean(mac), coverage_dev = mean(abs(cov - 0.90)), coverage_min = min(cov),
    coverage_pooled = mean(sqrt(e2) <= 1.645 * ex), ratio_pooled = sqrt(mean(e2)) / sqrt(mean(ex^2)),
    k = k, coverage_dev_k = mean(abs(cvk - 0.90)), coverage_min_k = min(cvk))
}
out <- list(); wins <- list()
for (it in ITS) {
  f <- file.path(ROOT, "outputs", it, "tables", "calibration_outer_points_toc.csv")
  if (!file.exists(f)) { msg("skip %s (no outer points)", it); next }
  o <- read.csv(f)
  o$g <- paste(zone(o$depth), band(o$lat))          # 12 strata (depth zone x latitude band), non-empty ones only
  o$g1 <- zone(o$depth); o$g2 <- band(o$lat)
  MS <- sub("^exp_", "", grep("^exp_C[0-9]$", names(o), value = TRUE))
  regs <- sort(unique(o$region))
  idx <- lapply(regs, function(r) which(o$region == r))
  bs <- replicate(B, unlist(idx[sample(seq_along(regs), length(regs), replace = TRUE)]), simplify = FALSE)
  for (m in MS) {
    ex <- o[[paste0("exp_", m)]]
    g7 <- c(o$g1, o$g2); e7 <- c(o$e2, o$e2); x7 <- c(ex, ex)    # 4 depth zones + 3 latitude bands, as in 07
    pt <- crit(e7, x7, g7)
    bo <- vapply(bs, function(i) crit(c(o$e2[i], o$e2[i]), c(ex[i], ex[i]), c(o$g1[i], o$g2[i])), numeric(8))
    ci <- apply(bo, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
    out[[length(out) + 1]] <- data.frame(iteration = it, method = m, t(pt),
      MACE_strata_lo = ci[1, 1], MACE_strata_hi = ci[2, 1], coverage_dev_lo = ci[1, 2], coverage_dev_hi = ci[2, 2],
      coverage_pooled_lo = ci[1, 4], coverage_pooled_hi = ci[2, 4],
      coverage_dev_k_lo = ci[1, 7], coverage_dev_k_hi = ci[2, 7])
    wins[[paste(it, m)]] <- bo
  }
  # how often is each method the best on each criterion, across bootstrap replicates?
  for (cr in c("MACE_strata", "coverage_dev", "coverage_dev_k")) {
    M <- sapply(MS, function(m) wins[[paste(it, m)]][cr, ])
    w <- table(factor(MS[apply(M, 1, which.min)], levels = MS)) / B
    msg("%s | %s: %s", it, cr, paste(sprintf("%s %.2f", MS, as.numeric(w)), collapse = "  "))
  }
}
R <- do.call(rbind, out)
write.csv(R, file.path(DIRS$compare, "calibration_selection.csv"), row.names = FALSE)
print(R, digits = 3, row.names = FALSE)
msg("calibration selection written")
