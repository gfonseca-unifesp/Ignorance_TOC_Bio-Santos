# =============================================================================
# 08v_decision_margin.R — a minimum margin for the decision rule, from boxplot statistics of block-level differences
# =============================================================================
# 08g decides from region-level differences with block-bootstrap CIs. With ~130-300 test blocks per region, tiny
# differences (G - R = 0.004) become "significant" without being meaningful. Here each region's test blocks (300 km)
# give a distribution of paired differences in RMSE (blocks with >= 5 test cells):
#   G - R    > 0: the regional model beats the global model with regional data (regionalize?)
#   T - G    > 0: the region's data reduce transfer error (collect?)
#   G50 - G  > 0: the second half of the region's data still helps (learning curve not flat?)
# Boxplot statistics per region: median, hinges (Q1, Q3), IQR, notch = median +/- 1.58 IQR / sqrt(n)
# (McGill et al. 1978), share of blocks with a positive difference. Candidate rules for "regionalize":
#   R0 current   bootstrap 95% CI of the region-level G - R excludes 0 (08g)
#   R1 notch     the notch of the block differences lies above 0
#   R2 hinge     Q1 > 0: the regional model wins in at least 75% of blocks
#   R3 notch + practical margin: notch above 0 and median >= 5% of the global model's RMSE in the region
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2) })
ITS <- c("iter1c_clean", "iter2c_clean")
bstat <- function(x) { q <- quantile(x, c(0.25, 0.5, 0.75)); iqr <- q[[3]] - q[[1]]; n <- length(x)
  c(n_blocks = n, median = q[[2]], q1 = q[[1]], q3 = q[[3]], iqr = iqr, notch_lo = q[[2]] - 1.58 * iqr / sqrt(n),
    notch_hi = q[[2]] + 1.58 * iqr / sqrt(n), pct_positive = 100 * mean(x > 0)) }
res <- list(); blk_all <- list()
for (it in ITS) {
  p <- read.csv(file.path(ROOT, "outputs", it, "tables", "regionalization_points_toc.csv"))
  rg <- read.csv(file.path(ROOT, "outputs", it, "tables", "regionalization_by_region_toc.csv"))
  rm_ <- function(e) sqrt(mean(e^2))
  blk <- do.call(rbind, lapply(split(p, list(p$region, p$block), drop = TRUE), function(x)
    data.frame(region = x$region[1], block = x$block[1], n = nrow(x), T = rm_(x$T - x$obs), G = rm_(x$G - x$obs), R = rm_(x$R - x$obs), G50 = rm_(x$G50 - x$obs))))
  blk <- blk[blk$n >= 5, ]
  blk$GR <- blk$G - blk$R; blk$TG <- blk$T - blk$G; blk$G50G <- blk$G50 - blk$G
  blk$iteration <- it; blk_all[[it]] <- blk
  for (r in sort(unique(blk$region))) {
    b <- blk[blk$region == r, ]; lab <- rg$region[r]
    gG <- rg$RMSE_global_G[r]
    for (m in c("GR", "TG", "G50G")) {
      s <- bstat(b[[m]])
      res[[length(res) + 1]] <- data.frame(iteration = it, region = lab, difference = m, t(s), region_RMSE_G = gG,
        practical_margin_5pct = 0.05 * gG, decision_08g = rg$decision[r],
        bootstrap_ci_lo = switch(m, GR = rg$G_minus_R_lo[r], TG = rg$T_minus_G_lo[r], G50G = rg$G50_minus_G_lo[r]),
        bootstrap_ci_hi = switch(m, GR = rg$G_minus_R_hi[r], TG = rg$T_minus_G_hi[r], G50G = rg$G50_minus_G_hi[r]))
    }
  }
}
res <- do.call(rbind, res)
gr <- res[res$difference == "GR", ]
gr$R0_bootstrap <- gr$bootstrap_ci_lo > 0
gr$R1_notch <- gr$notch_lo > 0
gr$R2_hinge <- gr$q1 > 0
gr$R3_notch_practical <- gr$notch_lo > 0 & gr$median >= gr$practical_margin_5pct
write.csv(res, file.path(DIRS$compare, "decision_margin_boxplot_stats.csv"), row.names = FALSE)
write.csv(gr, file.path(DIRS$compare, "decision_margin_regionalize_rules.csv"), row.names = FALSE)
print(gr[, c("iteration", "region", "n_blocks", "median", "q1", "q3", "notch_lo", "notch_hi", "pct_positive", "practical_margin_5pct",
             "R0_bootstrap", "R1_notch", "R2_hinge", "R3_notch_practical")], digits = 3, row.names = FALSE)
cat("\n--- T - G and G50 - G: notch above 0? ---\n")
tg <- res[res$difference != "GR", ]; tg$notch_above_0 <- tg$notch_lo > 0; tg$hinge_above_0 <- tg$q1 > 0
print(tg[, c("iteration", "region", "difference", "median", "q1", "notch_lo", "pct_positive", "notch_above_0", "hinge_above_0", "bootstrap_ci_lo")], digits = 3, row.names = FALSE)

B <- do.call(rbind, blk_all)
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
B$region_lab <- factor(RN[B$region], levels = rev(RN))
B$iteration <- factor(ifelse(B$iteration == "iter1c_clean", "iteration 1", "iteration 2"))
pm <- data.frame(region_lab = factor(RN, levels = rev(RN)), m = 0.05 * read.csv(file.path(ROOT, "outputs", "iter1c_clean", "tables", "regionalization_by_region_toc.csv"))$RMSE_global_G)
fig <- ggplot(B, aes(GR, region_lab, fill = iteration)) + geom_vline(xintercept = 0, colour = "grey40") +
  geom_boxplot(notch = TRUE, outlier.size = 0.4, outlier.alpha = 0.4, position = position_dodge(width = 0.75), width = 0.6) +
  geom_point(data = pm, aes(x = m, y = region_lab), inherit.aes = FALSE, shape = 124, size = 4, colour = "#C2410C") +
  coord_cartesian(xlim = c(-0.15, 0.15)) + scale_fill_manual(values = c("#9FC9CE", "#E8C9A6"), name = NULL) +
  labs(x = "block-level RMSE difference G - R (log10 TOC); > 0 = regional model better; orange tick = 5% of the global model's RMSE",
       y = NULL, title = "Does a regional model beat the global model with regional data?",
       subtitle = "notched boxplots of paired differences over 300-km test blocks (>= 5 cells); notches ~95% interval of the median") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(DIRS$compare, "fig_decision_margin.png"), fig, width = 200, height = 150, units = "mm", dpi = 250, bg = "white")
msg("decision margin done")
