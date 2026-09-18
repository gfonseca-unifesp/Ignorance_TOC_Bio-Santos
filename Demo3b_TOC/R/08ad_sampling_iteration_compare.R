# =============================================================================
# 08ad_sampling_iteration_compare.R — F5.3: does the collection test still say the same
# thing on the final map (iter2c_clean) as on the base model (iter1c_clean)?
# =============================================================================
# Two things change at once between the published test and the new one:
#   (a) the iteration: iter1c_clean (base + carbon supply) -> iter2c_clean (+ sample-based lithology);
#   (b) the number of random starts per region: 3 -> 10.
# Starts are seeded by region and start index, so starts 1-3 of the 10-start run are the same runs a
# 3-start run would give. Restricting the new runs to reps 1-3 therefore isolates (a) from (b).
# Reads only existing run tables; writes to outputs/comparison/ (nothing is overwritten).
Sys.setenv(DEMO3B_ITER = "iter2c_clean")
source("R/00_config.R")
I1 <- "iter1c_clean"; I2 <- "iter2c_clean"
tb <- function(it, f) read.csv(file.path(ROOT, "outputs", it, "tables", f))

# --- recompute the pooled and per-region summaries from the raw runs, for any subset of starts ------
summarise_runs <- function(R, reps) {
  R <- R[R$rep %in% reps, ]
  st <- R[R$pct_pool == 20, c("region", "rep", "strategy", "rmse")]; names(st)[4] <- "rmse_start"
  G <- merge(R[R$pct_pool > 20, ], st, by = c("region", "rep", "strategy"))
  G$gain <- G$rmse_start - G$rmse
  rnd <- G[G$strategy == "random", c("region", "rep", "pct_pool", "gain")]; names(rnd)[4] <- "gain_random"
  G <- merge(G, rnd, by = c("region", "rep", "pct_pool"))
  G$gain_minus_random <- G$gain - G$gain_random
  set.seed(CFG$seed)
  pooled <- do.call(rbind, lapply(split(G, list(G$strategy, G$pct_pool), drop = TRUE), function(x) {
    b <- replicate(2000, mean(sample(x$gain_minus_random, replace = TRUE)))
    data.frame(strategy = x$strategy[1], pct_pool = x$pct_pool[1], n_runs = nrow(x),
               mean_gain = mean(x$gain), mean_gain_minus_random = mean(x$gain_minus_random),
               ci_lo = quantile(b, 0.025), ci_hi = quantile(b, 0.975),
               ratio_vs_random = mean(x$gain) / mean(x$gain_random[x$strategy == x$strategy[1]]),
               regions_better_than_random = sum(tapply(x$gain_minus_random, x$region, mean) > 0))
  }))
  list(pooled = pooled, by_region = aggregate(cbind(rmse, gain, gain_minus_random) ~ region + strategy + pct_pool, G, mean))
}

# --- the decision rule of 08l, applied to any per-region table -------------------------------------
LAB <- c(novelty = "directed: environmental space", uncertainty = "directed: model uncertainty", space = "directed: geographic space")
recommend <- function(h7) do.call(rbind, lapply(split(h7, h7$region), function(x) {
  s <- do.call(rbind, lapply(names(LAB), function(k) { y <- x[x$strategy == k, ]
    data.frame(strategy = k, adv40 = y$gain_minus_random[y$pct_pool == 40], adv60 = y$gain_minus_random[y$pct_pool == 60]) }))
  s$mean_adv <- (s$adv40 + s$adv60) / 2
  ok <- s[s$adv40 > 0 & s$adv60 > 0, ]
  best <- if (nrow(ok)) ok[which.max(ok$mean_adv), ] else NULL
  data.frame(region = x$region[1], collection = if (is.null(best)) "random suffices" else unname(LAB[[best$strategy]]),
             advantage = if (is.null(best)) NA else best$mean_adv)
}))

R1 <- readRDS(file.path(ROOT, "outputs", I1, "models", "targeted_sampling_runs_toc.rds"))
R2 <- readRDS(file.path(ROOT, "outputs", I2, "models", "targeted_sampling_runs_toc.rds"))
s1  <- summarise_runs(R1, 1:3)    # published test
s2a <- summarise_runs(R2, 1:3)    # new iteration, same starts  -> isolates the iteration
s2b <- summarise_runs(R2, 1:10)   # new iteration, 10 starts    -> isolates the extra starts

pooled <- rbind(cbind(run = paste0(I1, " (3 starts)"),  s1$pooled),
                cbind(run = paste0(I2, " (3 starts)"),  s2a$pooled),
                cbind(run = paste0(I2, " (10 starts)"), s2b$pooled))
write.csv(pooled, file.path(DIRS$compare, "f53_targeted_sampling_pooled_compare.csv"), row.names = FALSE)
cat("\n===== pooled gain over random (log10 TOC RMSE) =====\n")
print(pooled[order(pooled$pct_pool, pooled$strategy), c("run", "strategy", "pct_pool", "n_runs", "mean_gain",
      "ratio_vs_random", "mean_gain_minus_random", "ci_lo", "ci_hi", "regions_better_than_random")], digits = 3, row.names = FALSE)

# --- per region: decision path (08g) and recommended strategy --------------------------------------
rule1 <- tb(I1, "regionalization_by_region_toc.csv"); rule2 <- tb(I2, "regionalization_by_region_toc.csv")
reg <- Reduce(function(a, b) merge(a, b, by = "region"), list(
  data.frame(region = rule1$region, decision_iter1c = rule1$decision, T_minus_G_iter1c = rule1$T_minus_G),
  data.frame(region = rule2$region, decision_iter2c = rule2$decision, T_minus_G_iter2c = rule2$T_minus_G),
  setNames(recommend(s1$by_region),  c("region", "strategy_iter1c_3", "adv_iter1c_3")),
  setNames(recommend(s2a$by_region), c("region", "strategy_iter2c_3", "adv_iter2c_3")),
  setNames(recommend(s2b$by_region), c("region", "strategy_iter2c_10", "adv_iter2c_10"))))
reg$decision_changed <- reg$decision_iter1c != reg$decision_iter2c
reg$strategy_changed_by_iteration <- reg$strategy_iter1c_3 != reg$strategy_iter2c_3
reg$strategy_changed_by_starts    <- reg$strategy_iter2c_3 != reg$strategy_iter2c_10
write.csv(reg, file.path(DIRS$compare, "f53_collection_by_region_compare.csv"), row.names = FALSE)
cat("\n===== decision path and recommended strategy, by region =====\n")
print(reg[, c("region", "decision_iter1c", "decision_iter2c", "strategy_iter1c_3", "strategy_iter2c_3", "strategy_iter2c_10")],
      row.names = FALSE, right = FALSE)
cat("\nchanged: decision", sum(reg$decision_changed), "/ 8 | strategy by iteration",
    sum(reg$strategy_changed_by_iteration), "/ 8 | strategy by extra starts", sum(reg$strategy_changed_by_starts), "/ 8\n")
msg("F5.3 comparison written to %s", DIRS$compare)
