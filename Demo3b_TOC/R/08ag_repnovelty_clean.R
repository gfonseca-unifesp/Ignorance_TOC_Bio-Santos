# =============================================================================
# 08ag_repnovelty_clean.R — MSv9 (pre-submission review, point M1): does representative novelty
# still win when it cannot see the evaluation cells, and when it is applied greedily?
# =============================================================================
# 08ae scored each candidate block with Equation 1 summed over `inr`, the region's observation cells
# — which include the cells of the fixed TEST blocks. No response value was used, but the criterion
# could see where the evaluation cells lie in predictor space, which the other strategies could not.
# It also ranked blocks once and took the top n (one shot), while the paper describes a greedy rule.
#
# This script repeats the same experiment (same regions, blocks, test blocks, starts and seeds) with
# three evaluation sets and two selection modes:
#   obs_all     the region's observation cells, as published            (reference, reproduces 08ae)
#   obs_notest  the same minus every cell of the fixed test blocks      (no sight of the evaluation)
#   grid        a sample of the region's PREDICTION grid cells          (what a survey planner has)
#   oneshot     rank once, take the top n                               (as published)
#   greedy      take one block, update the distances, score again       (as Section 2.8 describes)
# and reports the gain against random from the same start, and the paired difference against the
# largest-error benchmark of 08h.
#   DEMO3B_ITER=iter2c_clean Rscript R/08ag_repnovelty_clean.R
# Writes outputs/comparison/m1_*; nothing is overwritten.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(FNN); library(terra); library(sf) })
set.seed(CFG$seed)
t_start <- Sys.time()
N_STARTS <- as.integer(Sys.getenv("DEMO3B_SAMPLING_STARTS", "10"))
N_GRID   <- as.integer(Sys.getenv("DEMO3B_GRID_EVAL", "3000"))   # evaluation cells per region, grid variant

## ---- the setup of 08h / 08ae, reproduced line for line ------------------------------
model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)
km <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)
d$region <- km$cluster
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
d$rblock <- paste(d$region, floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic",
        "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
imp <- pmax(model$finalModel$variable.importance[vars], 0); imp <- imp / sum(imp)
Zt  <- scale(as.matrix(d[, vars]))
ctr <- attr(Zt, "scaled:center"); scl <- attr(Zt, "scaled:scale")
Zw  <- sweep(Zt, 2, imp, "*")

fit <- function(idx) ranger(x = d[idx, vars, drop = FALSE], y = d$y[idx], num.trees = 300,
  mtry = min(bt$mtry, length(vars)), min.node.size = bt$min.node.size, splitrule = "variance",
  num.threads = CFG$cores, seed = CFG$seed)
test_rmse <- function(m, te) sqrt(mean((predict(m, d[te, vars, drop = FALSE], type = "response", num.threads = CFG$cores)$predictions - d$y[te])^2))

## ---- prediction-grid cells of each region (the "grid" evaluation set) ----------------
msg("sampling the prediction grid")
stk <- rast(path_toc_stack())[[vars]]
oc  <- which(!is.na(values(stk[[1]])[, 1]))
set.seed(CFG$seed)
samp_cells <- sample(oc, min(length(oc), 120000))
xy  <- xyFromCell(stk, samp_cells)
V   <- as.matrix(stk[samp_cells])
ok  <- complete.cases(V); V <- V[ok, , drop = FALSE]; xy <- xy[ok, , drop = FALSE]
XYZ  <- to_xyz(xy[, 1], xy[, 2])                       # region of a grid cell = nearest k-means centroid
Dc   <- sapply(seq_len(nrow(km$centers)), function(k) rowSums(sweep(XYZ, 2, km$centers[k, ], "-")^2))
greg <- max.col(-Dc)
Zg  <- sweep(sweep(sweep(V, 2, ctr, "-"), 2, scl, "/"), 2, imp, "*")
msg("grid cells per region: %s", paste(table(greg), collapse = " "))

## ---- representative novelty on a given evaluation set --------------------------------
# E: rows of the evaluation set in the weighted predictor space; sampled: rows already sampled
select_repnov <- function(cand, sampled_idx, inr, n_add, E, greedy) {
  ci <- inr[d$rblock[inr] %in% cand]; cb <- d$rblock[ci]
  dmin <- get.knnx(Zw[sampled_idx, , drop = FALSE], E, k = 1)$nn.dist[, 1]
  blocks <- unique(cb); picked <- character(0)
  score <- function(b) { dnew <- get.knnx(Zw[ci[cb == b], , drop = FALSE], E, k = 1)$nn.dist[, 1]; sum(pmax(dmin - dnew, 0)) }
  if (!greedy) {
    sc <- vapply(blocks, score, numeric(1))
    return(names(sort(sc, decreasing = TRUE))[seq_len(n_add)])
  }
  for (j in seq_len(n_add)) {
    left <- setdiff(blocks, picked)
    sc <- vapply(left, score, numeric(1))
    b <- left[which.max(sc)]; picked <- c(picked, b)
    dmin <- pmin(dmin, get.knnx(Zw[ci[cb == b], , drop = FALSE], E, k = 1)$nn.dist[, 1])
  }
  picked
}

VARIANTS <- list(
  list(name = "obs_all_oneshot",    eval = "obs_all",    greedy = FALSE),
  list(name = "obs_notest_oneshot", eval = "obs_notest", greedy = FALSE),
  list(name = "obs_notest_greedy",  eval = "obs_notest", greedy = TRUE),
  list(name = "grid_greedy",        eval = "grid",       greedy = TRUE))

old <- read.csv(file.path(DIRS$tables, "targeted_sampling_runs_toc.csv"))     # 08h: random and benchmark
runs <- list()
for (r in sort(unique(d$region))) {
  t0 <- Sys.time()
  inr <- which(d$region == r); outr <- which(d$region != r)
  blk <- unique(d$rblock[inr]); set.seed(CFG$seed + r)
  test_b <- sample(blk, round(0.2 * length(blk))); pool <- setdiff(blk, test_b)
  te <- inr[d$rblock[inr] %in% test_b]
  notest <- inr[!(d$rblock[inr] %in% test_b)]
  gi <- which(greg == r)
  if (length(gi) > N_GRID) { set.seed(CFG$seed + r); gi <- sample(gi, N_GRID) }
  EVAL <- list(obs_all = Zw[inr, , drop = FALSE], obs_notest = Zw[notest, , drop = FALSE], grid = Zg[gi, , drop = FALSE])
  for (rep in seq_len(N_STARTS)) {
    set.seed(CFG$seed + 10 * r + rep)
    start_b <- sample(pool, max(1, round(0.2 * length(pool))))
    r0 <- old$rmse[old$region == RN[r] & old$rep == rep & old$strategy == "random" & old$pct_pool == 20]
    stopifnot(length(r0) == 1)
    for (v in VARIANTS) {
      set.seed(CFG$seed + 100 * r + 10 * rep + 6)
      samp <- start_b
      runs[[length(runs) + 1]] <- data.frame(region = RN[r], rep = rep, strategy = v$name, pct_pool = 20,
                                             n_blocks = length(samp), rmse = r0)
      for (target in c(0.4, 0.6)) {
        n_add <- round(target * length(pool)) - length(samp)
        if (n_add < 1) next
        sampled_idx <- c(outr, inr[d$rblock[inr] %in% samp])
        samp <- c(samp, select_repnov(setdiff(pool, samp), sampled_idx, inr, n_add, EVAL[[v$eval]], v$greedy))
        m <- fit(c(outr, inr[d$rblock[inr] %in% samp]))
        runs[[length(runs) + 1]] <- data.frame(region = RN[r], rep = rep, strategy = v$name, pct_pool = 100 * target,
                                               n_blocks = length(samp), rmse = test_rmse(m, te))
      }
    }
  }
  msg("region %d done (%.1f min)", r, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
NEW <- do.call(rbind, runs)
write.csv(NEW, file.path(DIRS$compare, "m1_repnovelty_clean_runs_toc.csv"), row.names = FALSE)

## ---- gains against random, and the paired difference against the benchmark ------------
f54 <- read.csv(file.path(DIRS$compare, "f54_repnovelty_runs_toc.csv"))
R <- rbind(old[, names(NEW)], f54[, names(NEW)], NEW)
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
             mean_gain = mean(x$gain), mean_gain_random = mean(x$gain_random),
             ratio_vs_random = mean(x$gain) / mean(x$gain_random),
             mean_gain_minus_random = mean(x$gain_minus_random),
             ci_lo = quantile(b, 0.025), ci_hi = quantile(b, 0.975),
             regions_better_than_random = sum(tapply(x$gain_minus_random, x$region, mean) > 0))
}))
write.csv(pooled, file.path(DIRS$compare, "m1_pooled_toc.csv"), row.names = FALSE)
print(pooled[order(pooled$pct_pool, -pooled$mean_gain), ], digits = 3, row.names = FALSE)

# paired against the largest-error benchmark, run by run
bm <- G[G$strategy == "oracle", c("region", "rep", "pct_pool", "gain")]; names(bm)[4] <- "gain_benchmark"
P <- merge(G[G$strategy != "oracle", ], bm, by = c("region", "rep", "pct_pool"))
P$minus_benchmark <- P$gain - P$gain_benchmark
paired <- do.call(rbind, lapply(split(P, list(P$strategy, P$pct_pool), drop = TRUE), function(x) {
  set.seed(CFG$seed)
  b <- replicate(2000, mean(sample(x$minus_benchmark, replace = TRUE)))
  data.frame(strategy = x$strategy[1], pct_pool = x$pct_pool[1], n_runs = nrow(x),
             mean_minus_benchmark = mean(x$minus_benchmark), ci_lo = quantile(b, 0.025), ci_hi = quantile(b, 0.975),
             runs_better_than_benchmark = mean(x$minus_benchmark > 0)) }))
write.csv(paired, file.path(DIRS$compare, "m1_paired_vs_benchmark_toc.csv"), row.names = FALSE)
print(paired[order(paired$pct_pool, -paired$mean_minus_benchmark), ], digits = 3, row.names = FALSE)

## ---- effort equivalence: how much random sampling matches the guided gain --------------
# the random learning curve is known at 20, 40 and 60% of the pool; linear interpolation in between
cur <- aggregate(gain ~ pct_pool, G[G$strategy == "random", ], mean)
eq <- do.call(rbind, lapply(unique(G$strategy), function(s) do.call(rbind, lapply(c(40, 60), function(p) {
  g <- mean(G$gain[G$strategy == s & G$pct_pool == p])
  x <- approx(c(0, cur$gain), c(20, cur$pct_pool), xout = g, rule = 2)$y
  data.frame(strategy = s, pct_pool = p, gain = g, random_pct_pool_for_same_gain = x)
}))))
write.csv(eq, file.path(DIRS$compare, "m1_effort_equivalence_toc.csv"), row.names = FALSE)
print(eq[order(eq$pct_pool, -eq$gain), ], digits = 3, row.names = FALSE)
msg("M1 done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
