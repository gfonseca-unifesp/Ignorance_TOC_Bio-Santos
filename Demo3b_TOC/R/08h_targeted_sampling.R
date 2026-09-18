# =============================================================================
# 08h_targeted_sampling.R — H7: can new data be directed to where they add most to the model?
# =============================================================================
# Retrospective sampling experiment inside each region (the manuscript's falsifiable prediction:
# data placed where ignorance is flagged should reduce error more than data placed at random).
#   - The region's 300-km blocks are split into fixed test blocks (20%) and a candidate pool (80%).
#   - A random 20% of the pool starts as "already sampled"; the global model is trained on all other
#     regions + the sampled blocks.
#   - Blocks are added in two steps (to 40% and 60% of the pool) by one of five strategies. All but the
#     oracle use only information available before sampling:
#       random      random blocks (baseline)
#       space       geographic gap filling: sequentially, the block farthest from all data held so far
#       novelty     environmental gap filling: blocks most dissimilar to the training data in predictor space
#                   (importance-weighted standardised nearest-neighbour distance, the logic of the DI)
#       uncertainty blocks with the widest 80% prediction interval of a quantile regression forest
#       oracle      blocks with the largest current absolute error (uses the TOC values; an upper bound only)
#   - Skill = RMSE on the fixed test blocks; N random starts per region (DEMO3B_SAMPLING_STARTS, default 3;
#     F5.3 of ROADMAP_MSv6 runs iter2c_clean with 10). Seeds depend on region and start, so the first three
#     starts of a 10-start run are the same runs as a 3-start run.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1_baseline"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(FNN); library(ggplot2) })
set.seed(CFG$seed)

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz   <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster   # same as 05 / 07 / 08g
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
d$rblock <- paste(d$region, floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")

imp <- pmax(model$finalModel$variable.importance[vars], 0); imp <- imp / sum(imp)
Zw  <- sweep(scale(as.matrix(d[, vars])), 2, imp, "*")
X3  <- to_xyz(d$x_lon, d$y_lat)
STRATS <- c("random", "space", "novelty", "uncertainty", "oracle")
N_STARTS <- as.integer(Sys.getenv("DEMO3B_SAMPLING_STARTS", "3"))   # F5.3: random starts per region

fit <- function(idx) ranger(x = d[idx, vars, drop = FALSE], y = d$y[idx], num.trees = 300,
  mtry = min(bt$mtry, length(vars)), min.node.size = bt$min.node.size, splitrule = "variance",
  quantreg = TRUE, num.threads = CFG$cores, seed = CFG$seed)
test_rmse <- function(m, te) sqrt(mean((predict(m, d[te, vars, drop = FALSE], type = "response", num.threads = CFG$cores)$predictions - d$y[te])^2))

select_blocks <- function(s, cand, sampled_idx, m, inr, n_add) {
  if (s == "random") return(sample(cand, n_add))
  ci <- inr[d$rblock[inr] %in% cand]; cb <- d$rblock[ci]
  if (s == "space") {
    dmin <- chord_km(get.knnx(X3[sampled_idx, , drop = FALSE], X3[ci, , drop = FALSE], k = 1)$nn.dist[, 1])
    picked <- character(0)
    for (j in seq_len(n_add)) {
      bs <- tapply(dmin, cb, mean); bs <- bs[setdiff(names(bs), picked)]
      b <- names(bs)[which.max(bs)]; picked <- c(picked, b)
      dmin <- pmin(dmin, chord_km(get.knnx(X3[ci[cb == b], , drop = FALSE], X3[ci, , drop = FALSE], k = 1)$nn.dist[, 1]))
    }
    return(picked)
  }
  sc <- switch(s,
    novelty     = get.knnx(Zw[sampled_idx, , drop = FALSE], Zw[ci, , drop = FALSE], k = 1)$nn.dist[, 1],
    uncertainty = { q <- predict(m, d[ci, vars, drop = FALSE], type = "quantiles", quantiles = c(0.1, 0.9), num.threads = CFG$cores)$predictions; q[, 2] - q[, 1] },
    oracle      = abs(predict(m, d[ci, vars, drop = FALSE], type = "response", num.threads = CFG$cores)$predictions - d$y[ci]))
  bs <- tapply(sc, cb, mean)
  names(sort(bs, decreasing = TRUE))[seq_len(n_add)]
}

f_runs <- file.path(DIRS$models, "targeted_sampling_runs_toc.rds")
if (!file.exists(f_runs)) {
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
      m0 <- fit(c(outr, inr[d$rblock[inr] %in% start_b])); r0 <- test_rmse(m0, te)
      for (s in STRATS) {
        set.seed(CFG$seed + 100 * r + 10 * rep + match(s, STRATS))
        samp <- start_b; m <- m0
        runs[[length(runs) + 1]] <- data.frame(region = RN[r], rep = rep, strategy = s, pct_pool = 20, n_blocks = length(samp), rmse = r0)
        for (target in c(0.4, 0.6)) {
          n_add <- round(target * length(pool)) - length(samp)
          if (n_add < 1) next
          sampled_idx <- c(outr, inr[d$rblock[inr] %in% samp])
          samp <- c(samp, select_blocks(s, setdiff(pool, samp), sampled_idx, m, inr, n_add))
          m <- fit(c(outr, inr[d$rblock[inr] %in% samp]))
          runs[[length(runs) + 1]] <- data.frame(region = RN[r], rep = rep, strategy = s, pct_pool = 100 * target, n_blocks = length(samp), rmse = test_rmse(m, te))
        }
      }
    }
    msg("region %d done (%.1f min; %d blocks: %d test, %d pool)", r, as.numeric(difftime(Sys.time(), t0, units = "mins")), length(blk), length(test_b), length(pool))
  }
  saveRDS(do.call(rbind, runs), f_runs)
}
R <- readRDS(f_runs)
write.csv(R, file.path(DIRS$tables, "targeted_sampling_runs_toc.csv"), row.names = FALSE)

# gain = RMSE reduction from the common 20% start; efficiency = gain / gain of random (same region, rep, step)
st <- R[R$pct_pool == 20, c("region", "rep", "strategy", "rmse")]; names(st)[4] <- "rmse_start"
G <- merge(R[R$pct_pool > 20, ], st, by = c("region", "rep", "strategy"))
G$gain <- G$rmse_start - G$rmse
rnd <- G[G$strategy == "random", c("region", "rep", "pct_pool", "gain")]; names(rnd)[4] <- "gain_random"
G <- merge(G, rnd, by = c("region", "rep", "pct_pool"))
G$gain_minus_random <- G$gain - G$gain_random
summ <- aggregate(cbind(rmse, gain, gain_minus_random) ~ region + strategy + pct_pool, G, mean)
write.csv(summ, file.path(DIRS$tables, "targeted_sampling_by_region_toc.csv"), row.names = FALSE)
set.seed(CFG$seed)
pooled <- do.call(rbind, lapply(split(G, list(G$strategy, G$pct_pool), drop = TRUE), function(x) {
  b <- replicate(2000, mean(sample(x$gain_minus_random, replace = TRUE)))
  data.frame(strategy = x$strategy[1], pct_pool = x$pct_pool[1], mean_gain = mean(x$gain), mean_gain_minus_random = mean(x$gain_minus_random),
             ci_lo = quantile(b, 0.025), ci_hi = quantile(b, 0.975), regions_better_than_random = sum(tapply(x$gain_minus_random, x$region, mean) > 0))
}))
write.csv(pooled, file.path(DIRS$tables, "targeted_sampling_pooled_toc.csv"), row.names = FALSE)
print(pooled, digits = 3, row.names = FALSE)

lc <- aggregate(rmse ~ region + strategy + pct_pool, R, mean)
lc$strategy <- factor(lc$strategy, levels = STRATS)
p <- ggplot(lc, aes(pct_pool, rmse, colour = strategy)) + geom_line() + geom_point(size = 1.2) + facet_wrap(~region, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = c(random = "grey45", space = "#1D6A73", novelty = "#7F77DD", uncertainty = "#C2410C", oracle = "#111111")) +
  scale_x_continuous(breaks = c(20, 40, 60)) +
  labs(x = "% of the region's candidate blocks sampled", y = "RMSE on fixed test blocks (log10 TOC)",
       title = "Directed vs random addition of data (oracle = upper bound, uses TOC values)") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(DIRS$figs, "targeted_sampling_toc.png"), p, width = 220, height = 130, units = "mm", dpi = 250, bg = "white")
msg("targeted sampling done")
