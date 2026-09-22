# =============================================================================
# 08ae_repnovelty_sampling.R — F5.4 of ROADMAP_MSv6: does the criterion that ranks the map's sites
# also reduce ERROR, and not only dissimilarity?
# =============================================================================
# 08aa ranks sampling sites by REPRESENTATIVE NOVELTY: a site is worth more when it represents a large
# area of unsampled conditions. The gain curve of that figure measures environmental ignorance (DI) and
# is, by construction, favourable to novelty. The retrospective test of 08h measures ERROR, but it never
# included this criterion — so the map could only claim to reduce dissimilarity.
#
# This script adds the missing strategy to the SAME experiment. It reproduces the design of 08h exactly
# (same regions, same 300-km blocks, same test blocks, same random starts, same model settings, same
# seeds), runs only the new strategy, and merges the result with the runs 08h already produced:
#   repnovelty  the blocks whose addition most reduces the dissimilarity of the whole region, i.e.
#               sum over the region's cells of max(0, d_c - d(c, block)), where d_c is the distance
#               from cell c to the nearest cell already sampled, in the importance-weighted
#               standardised predictor space (the space the DI uses).
# As for novelty and uncertainty in 08h, the blocks of one step are the top n by that score (the score
# is recomputed from scratch at each step, after the previous blocks were added).
#   DEMO3B_ITER=iter2c_clean Rscript R/08ae_repnovelty_sampling.R
# Writes outputs/comparison/f54_*; nothing is overwritten.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(FNN) })
set.seed(CFG$seed)
t_start <- Sys.time()
N_STARTS <- as.integer(Sys.getenv("DEMO3B_SAMPLING_STARTS", "10"))
# MSv7: DEMO3B_FIG_ONLY=1 redraws the figure from the saved runs (no refit); DEMO3B_FIG_TAG suffixes its file name
FIG_ONLY <- Sys.getenv("DEMO3B_FIG_ONLY", "0") == "1"
FIG_TAG  <- Sys.getenv("DEMO3B_FIG_TAG", "")

## ---- the setup of 08h, reproduced line for line -------------------------------------
model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz   <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
d$rblock <- paste(d$region, floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic",
        "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
imp <- pmax(model$finalModel$variable.importance[vars], 0); imp <- imp / sum(imp)
Zw  <- sweep(scale(as.matrix(d[, vars])), 2, imp, "*")

fit <- function(idx) ranger(x = d[idx, vars, drop = FALSE], y = d$y[idx], num.trees = 300,
  mtry = min(bt$mtry, length(vars)), min.node.size = bt$min.node.size, splitrule = "variance",
  quantreg = TRUE, num.threads = CFG$cores, seed = CFG$seed)
test_rmse <- function(m, te) sqrt(mean((predict(m, d[te, vars, drop = FALSE], type = "response", num.threads = CFG$cores)$predictions - d$y[te])^2))

# representative novelty at block level
select_repnovelty <- function(cand, sampled_idx, inr, n_add) {
  ci <- inr[d$rblock[inr] %in% cand]; cb <- d$rblock[ci]
  dmin <- get.knnx(Zw[sampled_idx, , drop = FALSE], Zw[inr, , drop = FALSE], k = 1)$nn.dist[, 1]
  blocks <- unique(cb)
  sc <- vapply(blocks, function(b) {
    dnew <- get.knnx(Zw[ci[cb == b], , drop = FALSE], Zw[inr, , drop = FALSE], k = 1)$nn.dist[, 1]
    sum(pmax(dmin - dnew, 0))
  }, numeric(1))
  names(sort(sc, decreasing = TRUE))[seq_len(n_add)]
}

## ---- run the new strategy on the same starts ----------------------------------------
old <- read.csv(file.path(DIRS$tables, "targeted_sampling_runs_toc.csv"))     # 08h, same iteration
if (FIG_ONLY) {
  NEW <- read.csv(file.path(DIRS$compare, "f54_repnovelty_runs_toc.csv"))
  R <- rbind(old[, names(NEW)], NEW)
} else {
runs <- list()
for (r in sort(unique(d$region))) {
  t0 <- Sys.time()
  inr <- which(d$region == r); outr <- which(d$region != r)
  blk <- unique(d$rblock[inr]); set.seed(CFG$seed + r)
  test_b <- sample(blk, round(0.2 * length(blk))); pool <- setdiff(blk, test_b)
  te <- inr[d$rblock[inr] %in% test_b]
  for (rep in seq_len(N_STARTS)) {
    set.seed(CFG$seed + 10 * r + rep)
    start_b <- sample(pool, max(1, round(0.2 * length(pool))))
    # the start is the same run 08h already scored; take its RMSE from there instead of refitting
    r0 <- old$rmse[old$region == RN[r] & old$rep == rep & old$strategy == "random" & old$pct_pool == 20]
    stopifnot(length(r0) == 1)
    set.seed(CFG$seed + 100 * r + 10 * rep + 6)          # strategy index 6 (08h uses 1-5)
    samp <- start_b
    runs[[length(runs) + 1]] <- data.frame(region = RN[r], rep = rep, strategy = "repnovelty",
                                           pct_pool = 20, n_blocks = length(samp), rmse = r0)
    for (target in c(0.4, 0.6)) {
      n_add <- round(target * length(pool)) - length(samp)
      if (n_add < 1) next
      sampled_idx <- c(outr, inr[d$rblock[inr] %in% samp])
      samp <- c(samp, select_repnovelty(setdiff(pool, samp), sampled_idx, inr, n_add))
      m <- fit(c(outr, inr[d$rblock[inr] %in% samp]))
      runs[[length(runs) + 1]] <- data.frame(region = RN[r], rep = rep, strategy = "repnovelty",
                                             pct_pool = 100 * target, n_blocks = length(samp), rmse = test_rmse(m, te))
    }
  }
  msg("region %d done (%.1f min)", r, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
NEW <- do.call(rbind, runs)
write.csv(NEW, file.path(DIRS$compare, "f54_repnovelty_runs_toc.csv"), row.names = FALSE)

## ---- pooled comparison with the strategies of 08h -----------------------------------
R <- rbind(old[, names(NEW)], NEW)
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
             mean_gain = mean(x$gain), ratio_vs_random = mean(x$gain) / mean(x$gain_random),
             mean_gain_minus_random = mean(x$gain_minus_random),
             ci_lo = quantile(b, 0.025), ci_hi = quantile(b, 0.975),
             regions_better_than_random = sum(tapply(x$gain_minus_random, x$region, mean) > 0))
}))
write.csv(pooled, file.path(DIRS$compare, "f54_repnovelty_pooled_toc.csv"), row.names = FALSE)
print(pooled[order(pooled$pct_pool, -pooled$mean_gain), ], digits = 3, row.names = FALSE)

by_reg <- aggregate(gain_minus_random ~ region + strategy + pct_pool, G, mean)
write.csv(by_reg, file.path(DIRS$compare, "f54_repnovelty_by_region_toc.csv"), row.names = FALSE)
msg("F5.4 done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
}

## ---- figure: the learning curves of 08h with the new strategy added (Supplementary Figure) --------
suppressPackageStartupMessages(library(ggplot2))
STR6 <- c("random", "space", "novelty", "repnovelty", "uncertainty", "oracle")
lc <- aggregate(rmse ~ region + strategy + pct_pool, R, mean)
lc$strategy <- factor(lc$strategy, levels = STR6)
p <- ggplot(lc, aes(pct_pool, rmse, colour = strategy)) + geom_line() + geom_point(size = 1.2) +
  facet_wrap(~region, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = c(random = "grey45", space = "#1D6A73", novelty = "#7F77DD",
                                 repnovelty = "#08306B", uncertainty = "#C2410C", oracle = "#111111"),
                      labels = c(random = "random", space = "geographic gaps", novelty = "environmental novelty",
                                 repnovelty = "representative novelty", uncertainty = "model uncertainty",
                                 oracle = "largest-error benchmark (uses TOC values)")) +
  scale_x_continuous(breaks = c(20, 40, 60)) +
  labs(x = "% of the region's candidate blocks sampled", y = "RMSE on fixed test blocks (log10 TOC)", colour = NULL,
       title = "Directed vs random addition of data, 8 regions x 10 random starts (iteration 2)") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(DIRS$compare, paste0("f54_targeted_sampling_curves_toc", FIG_TAG, ".png")), p, width = 220, height = 130, units = "mm", dpi = 250, bg = "white")
msg("learning-curve figure with six strategies written")
