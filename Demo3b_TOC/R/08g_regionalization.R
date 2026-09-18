# =============================================================================
# 08g_regionalization.R — interpolation vs transfer error by region:
# collect data and load it into the global model, regionalize the model, or neither?
# =============================================================================
# For each of the 8 regions, the SAME within-region test blocks (300-km Equal Earth blocks,
# up to 5 folds) are predicted by:
#   T   transfer : global model trained without any data from the region (trust in a data-less region)
#   G   global   : global model trained on all other regions + the region's other folds (interpolation)
#   R   regional : model trained only on the region's other folds (a regional model)
#   G25, G50     : global model with 25% / 50% of the region's training blocks (learning curve:
#                  does more data from the region keep reducing error?)
# Noise floor per region: median within-cell SD of log10 TOC (cells with >= 3 sites).
# Uncertainty: bootstrap over test blocks (1000 resamples) of RMSE differences.
# Decision per region, evaluated in this order:
#   stop (L5)                  reducible share of G's squared error < 25%
#   regionalize (L2)           R better than G (CI of G - R above 0)
#   collect and load data (L3) T worse than G (CI of T - G above 0) and G50 worse than G (more data still helps)
#   new drivers / resolution   none of the above, with error well above the noise floor
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1_baseline"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(CAST); library(ggplot2); library(patchwork) })
set.seed(CFG$seed)

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster   # same as 05 / 07
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
d$rblock <- paste(d$region, floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")

fit <- function(idx) ranger(x = d[idx, vars, drop = FALSE], y = d$y[idx], num.trees = 300,
  mtry = min(bt$mtry, length(vars)), min.node.size = bt$min.node.size, splitrule = "variance",
  num.threads = CFG$cores, seed = CFG$seed)
prd <- function(m, idx) predict(m, d[idx, vars, drop = FALSE], num.threads = CFG$cores)$predictions

f_pts <- file.path(DIRS$models, "regionalization_points_toc.rds")
if (!file.exists(f_pts)) {
  pts <- list()
  for (r in sort(unique(d$region))) {
    t0 <- Sys.time()
    inr <- which(d$region == r); outr <- which(d$region != r)
    blk <- unique(d$rblock[inr]); k <- min(5, length(blk))
    set.seed(CFG$seed + r); bfold <- setNames(sample(rep_len(seq_len(k), length(blk))), blk)
    fold <- bfold[d$rblock[inr]]
    mT <- fit(outr)
    for (kk in seq_len(k)) {
      te <- inr[fold == kk]; trr <- inr[fold != kk]
      if (length(te) < 5 || length(trr) < 20) next
      res <- data.frame(region = r, fold = kk, block = d$rblock[te], obs = d$y[te],
                        T = prd(mT, te), G = prd(fit(c(outr, trr)), te), R = prd(fit(trr), te))
      tb <- unique(d$rblock[trr])
      for (fr in c(25, 50)) {
        p <- sapply(1:2, function(draw) { set.seed(CFG$seed + 100 * r + 10 * kk + draw)
          keep <- sample(tb, max(1, round(length(tb) * fr / 100)))
          prd(fit(c(outr, trr[d$rblock[trr] %in% keep])), te) })
        res[[paste0("G", fr)]] <- rowMeans(p)
      }
      pts[[length(pts) + 1]] <- res
    }
    msg("region %d done (%.1f min, %d blocks, %d folds)", r, as.numeric(difftime(Sys.time(), t0, units = "mins")), length(blk), k)
  }
  saveRDS(do.call(rbind, pts), f_pts)
}
P <- readRDS(f_pts)
write.csv(P, file.path(DIRS$tables, "regionalization_points_toc.csv"), row.names = FALSE)

rmse <- function(p, o) sqrt(mean((p - o)^2))
MOD <- c("T", "G25", "G50", "G", "R")
tab <- do.call(rbind, lapply(split(P, P$region), function(x) {
  r <- x$region[1]; floor_r <- median(d$y_sd[d$region == r], na.rm = TRUE)
  est <- sapply(MOD, function(m) rmse(x[[m]], x$obs))
  bl <- unique(x$block); bi <- split(seq_len(nrow(x)), x$block)
  set.seed(CFG$seed + r)
  bs <- t(replicate(1000, { i <- unlist(bi[sample(bl, length(bl), replace = TRUE)], use.names = FALSE)
    s <- sapply(MOD, function(m) rmse(x[[m]][i], x$obs[i])); c(TG = s[["T"]] - s[["G"]], GR = s[["G"]] - s[["R"]], G50G = s[["G50"]] - s[["G"]]) }))
  ci <- apply(bs, 2, quantile, c(0.025, 0.975))
  reducible <- 100 * (1 - floor_r^2 / est[["G"]]^2)
  decision <- if (reducible < 25) "stop: near noise floor (L5)" else
              if (ci[1, "GR"] > 0) "regionalize the model (L2)" else
              if (ci[1, "TG"] > 0 && ci[1, "G50G"] > 0) "collect data and load into the global model (L3)" else
              "new drivers or finer resolution (L2/L4)"
  data.frame(region = RN[r], n_test = nrow(x), n_blocks = length(bl), noise_floor = floor_r,
             RMSE_transfer_T = est[["T"]], RMSE_G25 = est[["G25"]], RMSE_G50 = est[["G50"]], RMSE_global_G = est[["G"]], RMSE_regional_R = est[["R"]],
             T_minus_G = est[["T"]] - est[["G"]], T_minus_G_lo = ci[1, "TG"], T_minus_G_hi = ci[2, "TG"],
             G_minus_R = est[["G"]] - est[["R"]], G_minus_R_lo = ci[1, "GR"], G_minus_R_hi = ci[2, "GR"],
             G50_minus_G = est[["G50"]] - est[["G"]], G50_minus_G_lo = ci[1, "G50G"], G50_minus_G_hi = ci[2, "G50G"],
             reducible_pct_G = reducible, decision = decision)
}))
write.csv(tab, file.path(DIRS$tables, "regionalization_by_region_toc.csv"), row.names = FALSE)
print(tab[, c("region", "n_test", "noise_floor", "RMSE_transfer_T", "RMSE_G50", "RMSE_global_G", "RMSE_regional_R",
              "T_minus_G_lo", "G_minus_R_lo", "G50_minus_G_lo", "reducible_pct_G", "decision")], digits = 3, row.names = FALSE)

# --- figure ---------------------------------------------------------------------------------------------------
long <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) data.frame(region = tab$region[i],
  model = factor(c("transfer (no regional data)", "global + regional data", "regional model"), levels = c("transfer (no regional data)", "global + regional data", "regional model")),
  RMSE = c(tab$RMSE_transfer_T[i], tab$RMSE_global_G[i], tab$RMSE_regional_R[i]), floor = tab$noise_floor[i], decision = tab$decision[i])))
long$region_lab <- paste0(long$region, "\n", long$decision)
pa <- ggplot(long, aes(RMSE, region_lab, colour = model)) +
  geom_point(aes(x = floor), shape = "|", size = 5, colour = "grey45") +
  geom_point(size = 2.4, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("#C2410C", "#1D6A73", "#7F77DD"), name = NULL) +
  labs(x = "RMSE on the same within-region test blocks (log10 TOC); | = noise floor", y = NULL,
       title = "a  Transfer vs interpolation vs regional model, and the resulting decision") +
  theme_classic(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 9))
lc <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) data.frame(region = tab$region[i], fraction = c(0, 25, 50, 100),
  RMSE = c(tab$RMSE_transfer_T[i], tab$RMSE_G25[i], tab$RMSE_G50[i], tab$RMSE_global_G[i]))))
pb <- ggplot(lc, aes(fraction, RMSE, colour = region)) + geom_line() + geom_point(size = 1.5) +
  scale_x_continuous(breaks = c(0, 25, 50, 100)) +
  labs(x = "% of the region's training blocks added to the global model", y = "RMSE (log10 TOC)",
       title = "b  Learning curve: value of more data from each region") +
  theme_classic(base_size = 8) + theme(legend.position = "right", plot.title = element_text(face = "bold", size = 9))
ggsave(file.path(DIRS$figs, "regionalization_toc.png"), pa / pb + plot_layout(heights = c(1.4, 1)),
       width = 200, height = 210, units = "mm", dpi = 250, bg = "white")
msg("regionalization done")
