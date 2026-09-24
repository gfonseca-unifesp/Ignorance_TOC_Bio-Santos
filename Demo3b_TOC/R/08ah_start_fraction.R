# =============================================================================
# 08ah_start_fraction.R — MSv9 (pre-submission review, point M4): does the advantage of directed
# collection depend on how much is already known?
# =============================================================================
# The paper says that the strategy which reduces ignorance depends on how much is known, and supports
# it with the contrast between the two cases — two systems that differ in many other ways. This script
# tests the claim INSIDE Case 1, by repeating the retrospective test from four starting fractions of
# the candidate pool (5, 10, 20 and 40%) and comparing representative novelty with random addition
# from the same start. Representative novelty is used in its clean form (evaluation set without the
# test blocks, greedy selection; 08ag).
#   DEMO3B_ITER=iter2c_clean Rscript R/08ah_start_fraction.R
# Writes outputs/comparison/m4_*; nothing is overwritten.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(FNN); library(sf) })
set.seed(CFG$seed)
t_start <- Sys.time()
N_STARTS <- as.integer(Sys.getenv("DEMO3B_SAMPLING_STARTS", "10"))
FRACS <- c(0.05, 0.10, 0.20, 0.40)
ADD   <- c(0.20, 0.40)          # percentage points of the pool added to the start

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
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
  num.threads = CFG$cores, seed = CFG$seed)
test_rmse <- function(m, te) sqrt(mean((predict(m, d[te, vars, drop = FALSE], type = "response", num.threads = CFG$cores)$predictions - d$y[te])^2))

select_repnov <- function(cand, sampled_idx, inr, n_add, E) {      # greedy, as 08ag
  ci <- inr[d$rblock[inr] %in% cand]; cb <- d$rblock[ci]
  dmin <- get.knnx(Zw[sampled_idx, , drop = FALSE], E, k = 1)$nn.dist[, 1]
  blocks <- unique(cb); picked <- character(0)
  for (j in seq_len(n_add)) {
    left <- setdiff(blocks, picked)
    sc <- vapply(left, function(b) {
      dnew <- get.knnx(Zw[ci[cb == b], , drop = FALSE], E, k = 1)$nn.dist[, 1]; sum(pmax(dmin - dnew, 0)) }, numeric(1))
    b <- left[which.max(sc)]; picked <- c(picked, b)
    dmin <- pmin(dmin, get.knnx(Zw[ci[cb == b], , drop = FALSE], E, k = 1)$nn.dist[, 1])
  }
  picked
}

runs <- list()
for (r in sort(unique(d$region))) {
  t0 <- Sys.time()
  inr <- which(d$region == r); outr <- which(d$region != r)
  blk <- unique(d$rblock[inr]); set.seed(CFG$seed + r)
  test_b <- sample(blk, round(0.2 * length(blk))); pool <- setdiff(blk, test_b)
  te <- inr[d$rblock[inr] %in% test_b]
  E  <- Zw[inr[!(d$rblock[inr] %in% test_b)], , drop = FALSE]     # evaluation set without the test blocks
  for (f in FRACS) {
    n_start <- max(1, round(f * length(pool)))
    if (n_start + round(ADD[1] * length(pool)) > length(pool)) next
    for (rep in seq_len(N_STARTS)) {
      set.seed(CFG$seed + 1000 * round(100 * f) + 10 * r + rep)
      start_b <- sample(pool, n_start)
      m0 <- fit(c(outr, inr[d$rblock[inr] %in% start_b])); r0 <- test_rmse(m0, te)
      for (s in c("random", "repnovelty")) {
        samp <- start_b
        runs[[length(runs) + 1]] <- data.frame(region = RN[r], rep = rep, start_frac = f, strategy = s,
                                               added_pp = 0, n_blocks = length(samp), rmse = r0)
        for (a in ADD) {
          n_add <- round((f + a) * length(pool)) - length(samp)
          if (n_add < 1 || n_add > length(setdiff(pool, samp))) next
          set.seed(CFG$seed + 1000 * round(100 * f) + 100 * r + 10 * rep + match(s, c("random", "repnovelty")))
          samp <- c(samp, if (s == "random") sample(setdiff(pool, samp), n_add)
                          else select_repnov(setdiff(pool, samp), c(outr, inr[d$rblock[inr] %in% samp]), inr, n_add, E))
          m <- fit(c(outr, inr[d$rblock[inr] %in% samp]))
          runs[[length(runs) + 1]] <- data.frame(region = RN[r], rep = rep, start_frac = f, strategy = s,
                                                 added_pp = 100 * a, n_blocks = length(samp), rmse = test_rmse(m, te))
        }
      }
    }
  }
  msg("region %d done (%.1f min)", r, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
R <- do.call(rbind, runs)
write.csv(R, file.path(DIRS$compare, "m4_start_fraction_runs_toc.csv"), row.names = FALSE)

st <- R[R$added_pp == 0, c("region", "rep", "start_frac", "strategy", "rmse")]; names(st)[5] <- "rmse_start"
G <- merge(R[R$added_pp > 0, ], st, by = c("region", "rep", "start_frac", "strategy"))
G$gain <- G$rmse_start - G$rmse
rnd <- G[G$strategy == "random", c("region", "rep", "start_frac", "added_pp", "gain")]; names(rnd)[5] <- "gain_random"
G <- merge(G[G$strategy != "random", ], rnd, by = c("region", "rep", "start_frac", "added_pp"))
G$gain_minus_random <- G$gain - G$gain_random
set.seed(CFG$seed)
pooled <- do.call(rbind, lapply(split(G, list(G$start_frac, G$added_pp), drop = TRUE), function(x) {
  b <- replicate(2000, mean(sample(x$gain_minus_random, replace = TRUE)))
  data.frame(start_frac = x$start_frac[1], added_pp = x$added_pp[1], n_runs = nrow(x),
             mean_gain_guided = mean(x$gain), mean_gain_random = mean(x$gain_random),
             ratio_vs_random = mean(x$gain) / mean(x$gain_random),
             mean_gain_minus_random = mean(x$gain_minus_random),
             ci_lo = quantile(b, 0.025), ci_hi = quantile(b, 0.975),
             regions_better_than_random = sum(tapply(x$gain_minus_random, x$region, mean) > 0)) }))
write.csv(pooled, file.path(DIRS$compare, "m4_start_fraction_pooled_toc.csv"), row.names = FALSE)
print(pooled[order(pooled$added_pp, pooled$start_frac), ], digits = 3, row.names = FALSE)
msg("M4 done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
