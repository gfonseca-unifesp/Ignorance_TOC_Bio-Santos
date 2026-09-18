# =============================================================================
# 08ab_region_sensitivity.R — F5.1 of ROADMAP_MSv6: does the decision rule survive the definition
# of the regions? (pending item 5 of the MSv5 review package)
# =============================================================================
# The eight regions of 08g come from k-means on 3-D chord coordinates of the observation cells, with
# one seed. The decision rule (collect / regionalize / new drivers / stop) is read region by region, so
# it inherits whatever that clustering decided. Here the whole test is repeated with
#   k = 6, k = 10, and three k-means seeds at k = 8,
# using exactly the test of 08g (same 300-km blocks, same folds, same model settings), but without the
# G25 step, which the decision does not use.
#   DEMO3B_ITER=iter2c_clean Rscript R/08ab_region_sensitivity.R
# Writes outputs/comparison/f51_region_sensitivity_*.csv; nothing is overwritten.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger) })
t_start <- Sys.time()

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
XYZ <- to_xyz(d$x_lon, d$y_lat)
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))

fit <- function(idx) ranger(x = d[idx, vars, drop = FALSE], y = d$y[idx], num.trees = 300,
  mtry = min(bt$mtry, length(vars)), min.node.size = bt$min.node.size, splitrule = "variance",
  num.threads = CFG$cores, seed = CFG$seed)
prd <- function(m, idx) predict(m, d[idx, vars, drop = FALSE], num.threads = CFG$cores)$predictions
rmse <- function(p, o) sqrt(mean((p - o)^2))

# the reference clustering: k = 8 with the project seed, exactly as 05 / 07 / 08g / 08h / 08aa
set.seed(CFG$seed)
ref_region <- kmeans(XYZ, centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic",
        "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")

## ---- the 08g test and decision, for one clustering -----------------------------------
run_config <- function(label, region) {
  d$region <- region
  d$rblock <- paste(d$region, floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
  out <- list()
  for (r in sort(unique(d$region))) {
    t0 <- Sys.time()
    inr <- which(d$region == r); outr <- which(d$region != r)
    blk <- unique(d$rblock[inr]); k <- min(5, length(blk))
    set.seed(CFG$seed + r); bfold <- setNames(sample(rep_len(seq_len(k), length(blk))), blk)
    fold <- bfold[d$rblock[inr]]
    mT <- fit(outr); pts <- list()
    for (kk in seq_len(k)) {
      te <- inr[fold == kk]; trr <- inr[fold != kk]
      if (length(te) < 5 || length(trr) < 20) next
      res <- data.frame(block = d$rblock[te], obs = d$y[te],
                        T = prd(mT, te), G = prd(fit(c(outr, trr)), te), R = prd(fit(trr), te))
      tb <- unique(d$rblock[trr])
      res$G50 <- rowMeans(sapply(1:2, function(draw) { set.seed(CFG$seed + 100 * r + 10 * kk + draw)
        keep <- sample(tb, max(1, round(length(tb) * 0.5)))
        prd(fit(c(outr, trr[d$rblock[trr] %in% keep])), te) }))
      pts[[length(pts) + 1]] <- res
    }
    x <- do.call(rbind, pts)
    if (is.null(x) || nrow(x) < 10) { msg("  %s region %d skipped (too few test points)", label, r); next }
    floor_r <- median(d$y_sd[d$region == r], na.rm = TRUE)
    MOD <- c("T", "G50", "G", "R"); est <- sapply(MOD, function(m) rmse(x[[m]], x$obs))
    bl <- unique(x$block); bi <- split(seq_len(nrow(x)), x$block)
    set.seed(CFG$seed + r)
    bs <- t(replicate(1000, { i <- unlist(bi[sample(bl, length(bl), replace = TRUE)], use.names = FALSE)
      s <- sapply(MOD, function(m) rmse(x[[m]][i], x$obs[i]))
      c(TG = s[["T"]] - s[["G"]], GR = s[["G"]] - s[["R"]], G50G = s[["G50"]] - s[["G"]]) }))
    ci <- apply(bs, 2, quantile, c(0.025, 0.975))
    reducible <- 100 * (1 - floor_r^2 / est[["G"]]^2)
    decision <- if (reducible < 25) "stop: near noise floor (L5)" else
                if (ci[1, "GR"] > 0) "regionalize the model (L2)" else
                if (ci[1, "TG"] > 0 && ci[1, "G50G"] > 0) "collect data and load into the global model (L3)" else
                "new drivers or finer resolution (L2/L4)"
    # which reference region does this one mostly overlap, and how much?
    tb_ref <- table(ref_region[inr]); ref_main <- as.integer(names(tb_ref)[which.max(tb_ref)])
    out[[length(out) + 1]] <- data.frame(
      config = label, region = r, n_cells = length(inr), n_test = nrow(x), n_blocks = length(bl),
      lon_centroid = round(mean(d$x_lon[inr]), 1), lat_centroid = round(mean(d$y_lat[inr]), 1),
      noise_floor = floor_r, RMSE_transfer_T = est[["T"]], RMSE_global_G = est[["G"]], RMSE_regional_R = est[["R"]],
      T_minus_G = est[["T"]] - est[["G"]], T_minus_G_lo = ci[1, "TG"], T_minus_G_hi = ci[2, "TG"],
      G_minus_R_lo = ci[1, "GR"], G50_minus_G_lo = ci[1, "G50G"], reducible_pct_G = reducible,
      decision = decision, ref_region = RN[ref_main],
      overlap_with_ref_pct = round(100 * max(tb_ref) / length(inr), 1))
    msg("  %s region %d done (%.1f min): %s", label, r, as.numeric(difftime(Sys.time(), t0, units = "mins")), decision)
  }
  do.call(rbind, out)
}

CONFIGS <- list(list(label = "k8 seed 20260915 (reference)", k = 8, seed = CFG$seed),
                list(label = "k6 seed 20260915",  k = 6,  seed = CFG$seed),
                list(label = "k10 seed 20260915", k = 10, seed = CFG$seed),
                list(label = "k8 seed 20260916",  k = 8,  seed = CFG$seed + 1),
                list(label = "k8 seed 20260917",  k = 8,  seed = CFG$seed + 2))
f_all <- file.path(DIRS$compare, "f51_region_sensitivity_by_region.csv")
res <- list()
for (cfg in CONFIGS) {
  msg("configuration: %s", cfg$label)
  set.seed(cfg$seed)
  reg <- kmeans(XYZ, centers = cfg$k, nstart = 25, iter.max = 100)$cluster
  res[[cfg$label]] <- run_config(cfg$label, reg)
}
ALL <- do.call(rbind, res)
write.csv(ALL, f_all, row.names = FALSE)

## ---- stability table -------------------------------------------------------------------
ref <- res[[CONFIGS[[1]]$label]]
dec_of_ref <- setNames(ref$decision, ref$ref_region)
stab <- do.call(rbind, lapply(res, function(x) {
  same <- x$decision == dec_of_ref[x$ref_region]          # same call as the reference region it overlaps most
  data.frame(config = x$config[1], n_regions = nrow(x),
             pct_regions_same_decision = round(100 * mean(same), 1),
             pct_cells_same_decision = round(100 * sum(x$n_cells[same]) / sum(x$n_cells), 1),
             n_collect = sum(grepl("^collect", x$decision)), n_regionalize = sum(grepl("^regionalize", x$decision)),
             n_new_drivers = sum(grepl("^new drivers", x$decision)), n_stop = sum(grepl("^stop", x$decision)),
             T_minus_G_min = round(min(x$T_minus_G), 3), T_minus_G_max = round(max(x$T_minus_G), 3),
             T_minus_G_median = round(median(x$T_minus_G), 3),
             pct_regions_TG_CI_above_zero = round(100 * mean(x$T_minus_G_lo > 0), 1))
}))
write.csv(stab, file.path(DIRS$compare, "f51_region_sensitivity_stability.csv"), row.names = FALSE)
print(stab, row.names = FALSE)
msg("F5.1 done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
