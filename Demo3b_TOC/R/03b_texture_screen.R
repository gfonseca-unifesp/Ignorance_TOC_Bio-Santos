# =============================================================================
# 03b_texture_screen.R — does sediment texture reduce error? (test of H3 before iteration 2)
# =============================================================================
# Greedy nested forward selection (03_model.R, mtry = 2, one variable at a time) added no
# texture variable. That is a weak test: with mtry = 2 a single new column is rarely tried,
# and the lithology fractions act as a composition. Here the hypothesis is tested directly,
# with paired comparisons on identical folds:
#   A. spatial-block CV (same folds as 03_model.R): iteration-1 predictors vs + each texture
#      variable, + texture groups (grain size, porosity, lithology), + all texture;
#   B. leave-region-out (the 8 regions of iteration 1, where the ignorance map failed):
#      RMSE, NW Atlantic RMSE and the sediment-type bias component of the ignorance budget.
Sys.setenv(DEMO3B_ITER = "iter2_texture")
source("R/00_config.R")
suppressPackageStartupMessages({ library(caret); library(ranger); library(CAST) })
set.seed(CFG$seed)

obs <- read.csv(file.path(DIRS$data, "toc_cells.csv"))
stk <- rast(path_toc_stack())
d <- cbind(obs, stk[obs$cell][, CANDIDATES])
d <- d[complete.cases(d[, c("y", CANDIDATES)]), ]
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
bs <- CFG$cv_block_km * 1000
d$block <- paste(floor(exy[, 1] / bs), floor(exy[, 2] / bs))
folds <- CreateSpacetimeFolds(d, spacevar = "block", k = CFG$cv_k, seed = CFG$seed)

base <- readRDS(file.path(ROOT, "outputs", "iter1_baseline", "models", "ffs_toc.rds"))$selectedvars
groups <- list(grain = c("gs_d50", "gs_d16", "litho_gs"), porosity = "porosity",
               lithology = paste0("litho_t", 1:6), all_texture = TEXTURE_VARS)
sets <- c(list(baseline = base), setNames(lapply(TEXTURE_VARS, function(v) c(base, v)), paste0("+", TEXTURE_VARS)),
          setNames(lapply(groups, function(g) c(base, g)), paste0("+", names(groups))))
NT <- 300

# --- A. spatial-block CV, per-fold RMSE -------------------------------------------------------
fit_fold <- function(v, trn, tst, mtry) {
  m <- ranger(x = d[trn, v], y = d$y[trn], num.trees = NT, mtry = mtry, min.node.size = 5,
              num.threads = CFG$cores, seed = CFG$seed)
  predict(m, d[tst, v], num.threads = CFG$cores)$predictions
}
resA <- do.call(rbind, lapply(names(sets), function(s) {
  v <- sets[[s]]; p <- length(v)
  do.call(rbind, lapply(unique(c(2, floor(sqrt(p)), ceiling(p / 3))), function(mt) {
    oof <- rep(NA_real_, nrow(d))
    for (i in seq_along(folds$index)) oof[folds$indexOut[[i]]] <- fit_fold(v, folds$index[[i]], folds$indexOut[[i]], mt)
    fr <- vapply(folds$indexOut, function(k) sqrt(mean((oof[k] - d$y[k])^2)), 0)
    data.frame(set = s, n_vars = p, mtry = mt, fold = seq_along(fr), RMSE_fold = fr,
               RMSE_pooled = sqrt(mean((oof - d$y)^2)), R2_pooled = 1 - sum((oof - d$y)^2) / sum((d$y - mean(d$y))^2))
  }))
}))
# best mtry per set, then paired differences against the baseline (at its own best mtry)
best <- aggregate(RMSE_pooled ~ set + mtry, resA, `[`, 1)
best <- best[order(best$set, best$RMSE_pooled), ]; best <- best[!duplicated(best$set), ]
resA_b <- merge(resA, best[, c("set", "mtry")])
b0 <- resA_b[resA_b$set == "baseline", c("fold", "RMSE_fold")]; names(b0)[2] <- "RMSE_base"
cmpA <- merge(resA_b, b0, by = "fold")
sumA <- do.call(rbind, lapply(split(cmpA, cmpA$set), function(x) data.frame(
  set = x$set[1], n_vars = x$n_vars[1], mtry = x$mtry[1], R2_pooled = x$R2_pooled[1], RMSE_pooled = x$RMSE_pooled[1],
  dRMSE_mean_fold = mean(x$RMSE_fold - x$RMSE_base), folds_improved = sum(x$RMSE_fold < x$RMSE_base))))
sumA <- sumA[order(sumA$RMSE_pooled), ]
write.csv(resA, file.path(DIRS$tables, "texture_screen_spatialCV_folds.csv"), row.names = FALSE)
write.csv(sumA, file.path(DIRS$tables, "texture_screen_spatialCV.csv"), row.names = FALSE)
print(sumA, digits = 3)

# --- B. leave-region-out: the test where iteration 1 failed ------------------------------------------
reg <- read.csv(file.path(ROOT, "outputs", "iter1_baseline", "tables", "calibration_outer_points_toc.csv"))[, c("region", "lon", "lat")]
key <- function(lon, lat) paste(round(lon, 3), round(lat, 3))
d$region <- reg$region[match(key(d$x_lon, d$y_lat), key(reg$lon, reg$lat))]
msg("cells matched to iteration-1 regions: %d of %d", sum(!is.na(d$region)), nrow(d))
rc <- read.csv(file.path(ROOT, "outputs", "iter1_baseline", "tables", "transfer_regions_toc.csv"))
nwa <- rc$region[which.min((rc$lon + 65)^2 + (rc$lat - 40)^2)]
lith <- with(d, { v <- as.matrix(d[, paste0("litho_t", 1:6)]); s <- rowSums(v)
  ifelse(is.na(s) | s < 0.5, "unclassified", paste0("type ", max.col(v, ties.method = "first"))) })
setsB <- sets[c("baseline", "+grain", "+porosity", "+lithology", "+all_texture")]
dd <- d[!is.na(d$region), ]; lith <- lith[!is.na(d$region)]
resB <- do.call(rbind, lapply(names(setsB), function(s) {
  v <- setsB[[s]]; mt <- best$mtry[best$set == s]; pr <- rep(NA_real_, nrow(dd))
  for (r in sort(unique(dd$region))) {
    tr <- dd$region != r
    m <- ranger(x = dd[tr, v], y = dd$y[tr], num.trees = NT, mtry = mt, min.node.size = 5, num.threads = CFG$cores, seed = CFG$seed)
    pr[!tr] <- predict(m, dd[!tr, v], num.threads = CFG$cores)$predictions
  }
  e <- pr - dd$y; rg <- ave(e, dd$region); rl <- ave(e, dd$region, lith)
  within <- mean((e - rl)^2); nf <- min(median(dd$y_sd, na.rm = TRUE)^2, within)
  data.frame(set = s, mtry = mt, RMSE_regionCV = sqrt(mean(e^2)),
             RMSE_NW_Atlantic = sqrt(mean(e[dd$region == nwa]^2)),
             mse_total = mean(e^2), mse_regional_offset = mean(rg^2), mse_sediment_type_bias = mean((rl - rg)^2),
             mse_unstructured_above_noise = within - nf, mse_noise = nf,
             rho_abs_error_litho_gs_NWA = cor(abs(e[dd$region == nwa]), dd$litho_gs[dd$region == nwa], method = "spearman"))
}))
write.csv(resB, file.path(DIRS$tables, "texture_screen_regionCV.csv"), row.names = FALSE)
print(resB, digits = 3)
msg("texture screen done (NW Atlantic = region %d)", nwa)
