# =============================================================================
# 08ac_stratum_interaction.R — F5.2 of ROADMAP_MSv6: does the RELATIONSHIP change between strata?
# (pending item 6 of the MSv5 review package; the last untested step of the structure test)
# =============================================================================
# The structure test (08w, 08x) showed that depth zone and latitude band structure the MAGNITUDE of the
# error and do not help as predictors, which sent them to the error model of the ignorance map rather
# than to the mean model. One step was never exercised: if the RELATIONSHIP between predictors and TOC
# itself changed between strata, a separate model per stratum would be warranted.
#
# Test, under leave-region-out (the 8 regions of 08g, so that the comparison is out of sample in space):
#   global     one model trained on all other regions
#   stratified one model per stratum, trained on the cells of that stratum in the other regions
#              (a fully saturated interaction: every split may differ between strata)
#   plus_flag  the global model with the stratum added as a predictor (interaction only where the
#              forest chooses to split on it) — reported for context, not used for the decision
# Strata: the 4 depth zones and the 3 latitude bands of 07_calibration.R / 08y.
#
# CRITERION DECLARED BEFORE RUNNING (ROADMAP F5.2): a separate model per stratum is warranted only if
# the out-of-sample gain over the global model exceeds 2% of its RMSE AND the bootstrap CI of the gain
# excludes zero. Anything less means one model, with the strata in the error model.
#   DEMO3B_ITER=iter2c_clean Rscript R/08ac_stratum_interaction.R
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger) })
t_start <- Sys.time()
GAIN_THRESHOLD_PCT <- 2                       # declared before running

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
d$block <- paste(floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))      # 300-km blocks, for the bootstrap
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic",
        "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")

zone <- function(z) as.character(cut(z, c(-Inf, 200, 1000, 3000, Inf),
                                     labels = c("shelf <=200 m", "slope 200-1000 m", "1000-3000 m", "abyss >3000 m")))
band <- function(lat) as.character(cut(abs(lat), c(0, 30, 60, 90), include.lowest = TRUE,
                                       labels = c("tropics 0-30", "mid-latitudes 30-60", "high latitudes 60-90")))
d$depth_zone <- zone(d$depth); d$latitude_band <- band(d$y_lat)

fit <- function(idx, extra = NULL) {
  X <- d[idx, vars, drop = FALSE]
  if (!is.null(extra)) X[[extra]] <- factor(d[[extra]][idx])
  ranger(x = X, y = d$y[idx], num.trees = 300, mtry = min(bt$mtry, ncol(X)),
         min.node.size = bt$min.node.size, splitrule = "variance", num.threads = CFG$cores, seed = CFG$seed)
}
prd <- function(m, idx, extra = NULL) {
  X <- d[idx, vars, drop = FALSE]
  if (!is.null(extra)) X[[extra]] <- factor(d[[extra]][idx], levels = levels(factor(d[[extra]])))
  predict(m, X, num.threads = CFG$cores)$predictions
}

## ---- leave-region-out predictions, three ways -----------------------------------------
res <- list()
for (pn in c("depth_zone", "latitude_band")) {
  for (r in sort(unique(d$region))) {
    t0 <- Sys.time()
    te <- which(d$region == r); tr <- which(d$region != r)
    p_glob <- prd(fit(tr), te)
    p_flag <- prd(fit(tr, extra = pn), te, extra = pn)
    p_strat <- rep(NA_real_, length(te))
    for (s in unique(d[[pn]][te])) {
      itr <- tr[d[[pn]][tr] == s]; ite <- which(d[[pn]][te] == s)
      if (length(itr) < 100) next                     # too few cells of that stratum outside the region
      p_strat[ite] <- prd(fit(itr), te[ite])
    }
    res[[length(res) + 1]] <- data.frame(partition = pn, region = RN[r], block = d$block[te],
                                         stratum = d[[pn]][te], obs = d$y[te],
                                         global = p_glob, stratified = p_strat, plus_flag = p_flag)
    msg("%s | region %d done (%.1f min; %d test cells, %d strata)", pn, r,
        as.numeric(difftime(Sys.time(), t0, units = "mins")), length(te), length(unique(d[[pn]][te])))
  }
}
P <- do.call(rbind, res)
P <- P[!is.na(P$stratified), ]
write.csv(P, file.path(DIRS$compare, "f52_stratum_interaction_points.csv"), row.names = FALSE)

## ---- gain of the stratified model over the global one -----------------------------------
rmse <- function(p, o) sqrt(mean((p - o)^2))
summarise <- function(x, label) {
  bl <- unique(x$block); bi <- split(seq_len(nrow(x)), x$block)
  set.seed(CFG$seed)
  bs <- t(replicate(1000, { i <- unlist(bi[sample(bl, length(bl), replace = TRUE)], use.names = FALSE)
    g <- rmse(x$global[i], x$obs[i])
    c(strat = 100 * (g - rmse(x$stratified[i], x$obs[i])) / g,
      flag  = 100 * (g - rmse(x$plus_flag[i], x$obs[i])) / g) }))
  g <- rmse(x$global, x$obs)
  data.frame(level = label, n = nrow(x), RMSE_global = g,
             RMSE_stratified = rmse(x$stratified, x$obs), RMSE_plus_flag = rmse(x$plus_flag, x$obs),
             gain_stratified_pct = 100 * (g - rmse(x$stratified, x$obs)) / g,
             gain_lo = quantile(bs[, "strat"], 0.025), gain_hi = quantile(bs[, "strat"], 0.975),
             gain_plus_flag_pct = 100 * (g - rmse(x$plus_flag, x$obs)) / g)
}
tab <- do.call(rbind, c(
  lapply(split(P, P$partition), function(x) summarise(x, paste("all strata:", x$partition[1]))),
  lapply(split(P, list(P$partition, P$stratum), drop = TRUE), function(x)
    summarise(x, paste0(x$partition[1], ": ", x$stratum[1])))))
tab$warrants_separate_model <- tab$gain_stratified_pct > GAIN_THRESHOLD_PCT & tab$gain_lo > 0
rownames(tab) <- NULL
write.csv(tab, file.path(DIRS$compare, "f52_stratum_interaction_gain.csv"), row.names = FALSE)
print(tab[, c("level", "n", "RMSE_global", "RMSE_stratified", "gain_stratified_pct", "gain_lo", "gain_hi",
              "gain_plus_flag_pct", "warrants_separate_model")], digits = 3, row.names = FALSE)
msg("criterion: gain > %d%% with CI above zero -> separate models warranted in %d of %d levels",
    GAIN_THRESHOLD_PCT, sum(tab$warrants_separate_model), nrow(tab))
msg("F5.2 done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
