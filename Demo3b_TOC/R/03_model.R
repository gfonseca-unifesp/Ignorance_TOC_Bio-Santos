# =============================================================================
# 03_model.R — random-forest regression of log10 TOC with spatial CV and ffs
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(caret); library(ranger); library(CAST); library(FNN) })
set.seed(CFG$seed)

obs <- read.csv(file.path(DIRS$data, "toc_cells.csv"))
stk <- rast(path_toc_stack())
d <- cbind(obs, stk[obs$cell][, CANDIDATES])
n0 <- nrow(d)
d <- d[complete.cases(d[, c("y", CANDIDATES)]), ]
msg("training cells: %d (dropped %d with incomplete predictors)", nrow(d), n0 - nrow(d))

# --- spatial folds ----------------------------------------------------------------------
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
bs <- CFG$cv_block_km * 1000
d$block <- paste(floor(exy[, 1] / bs), floor(exy[, 2] / bs))
folds <- CreateSpacetimeFolds(d, spacevar = "block", k = CFG$cv_k, seed = CFG$seed)
d$fold <- NA_integer_; for (i in seq_along(folds$indexOut)) d$fold[folds$indexOut[[i]]] <- i
msg("blocks: %d; fold sizes: %s", length(unique(d$block)), paste(tabulate(d$fold), collapse = " / "))
ctrl <- trainControl(method = "cv", index = folds$index, indexOut = folds$indexOut, savePredictions = "final", allowParallel = FALSE)

# --- forward feature selection -------------------------------------------------------------
f_ffs <- file.path(DIRS$models, "ffs_toc.rds")
base_it  <- if (!is.null(ITERS[[ITER]]$base)) ITERS[[ITER]]$base else "iter1_baseline"   # iteration whose selection is the starting set
base_ffs <- file.path(ROOT, "outputs", base_it, "models", "ffs_toc.rds")
if (!file.exists(f_ffs)) {
  t0 <- Sys.time()
  if (ITER != "iter1_baseline" && file.exists(base_ffs)) {
    # Hypothesis-driven integration: iteration-1 predictors + the pre-declared block of the iteration
    # (ITERS[[ITER]]$block; ANALYTIC_LOG), without data-driven selection. For texture, greedy nested
    # selection added no variable, and a paired screen (03b_texture_screen.R) showed no spatial-CV gain
    # but a modest leave-region-out gain. Adding the declared block avoids choosing the intervention
    # by its own test result.
    start <- readRDS(base_ffs)$selectedvars
    added <- setdiff(intersect(ITERS[[ITER]]$block, CANDIDATES), start)
    if (length(ITERS[[ITER]]$block)) stopifnot(length(added) > 0)   # no block = re-use the base selection (iter1c_clean)
    steps <- data.frame(step = 0:1, added = c("(iteration-1 predictors)", if (length(added)) paste(added, collapse = "; ") else "(none: base selection re-used)"))
    ffsmod <- list(selectedvars = c(start, added), steps = steps, type = "iter1_plus_declared_block")
  } else {
    ffsmod <- ffs(d[, CANDIDATES], d$y, method = "ranger", metric = "RMSE", trControl = ctrl,
                  tuneGrid = data.frame(mtry = 2, splitrule = "variance", min.node.size = 5),
                  num.trees = CFG$ffs_trees, seed = CFG$seed, verbose = FALSE)
  }
  msg("variable selection done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(ffsmod, f_ffs)
}
ffsmod <- readRDS(f_ffs)
sel <- ffsmod$selectedvars
msg("selected: %s", paste(sel, collapse = ", "))
if (!is.null(ffsmod$steps)) {
  write.csv(ffsmod$steps, file.path(DIRS$tables, "ffs_selection_toc.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(step = seq_along(sel), variable = sel, cv_RMSE = c(NA, ffsmod$selectedvars_perf)[seq_along(sel)]),
            file.path(DIRS$tables, "ffs_selection_toc.csv"), row.names = FALSE)
}

# --- final model and random-CV twin ------------------------------------------------------------
p <- length(sel)
grid <- expand.grid(mtry = sort(unique(pmin(p, c(2, floor(sqrt(p)), ceiling(p / 2), p)))),
                    splitrule = "variance", min.node.size = c(5, 10))
set.seed(CFG$seed)
model <- train(d[, sel], d$y, method = "ranger", metric = "RMSE", trControl = ctrl, tuneGrid = grid,
               num.trees = CFG$final_trees, importance = "permutation", num.threads = CFG$cores)
saveRDS(model, file.path(DIRS$models, "rf_spatialCV_toc.rds"))
set.seed(CFG$seed)
model_rand <- train(d[, sel], d$y, method = "ranger", metric = "RMSE",
                    trControl = trainControl(method = "cv", number = CFG$cv_k, savePredictions = "final"),
                    tuneGrid = model$bestTune, num.trees = CFG$final_trees, num.threads = CFG$cores)

perf <- function(m, label) {
  pr <- m$pred; e <- pr$pred - pr$obs
  back <- function(z) 10^z - CFG$toc_offset
  data.frame(cv = label, R2_log = 1 - sum(e^2) / sum((pr$obs - mean(pr$obs))^2), RMSE_log = sqrt(mean(e^2)), MAE_log = mean(abs(e)),
             R2_pct = cor(back(pr$pred), back(pr$obs))^2, MAE_pct = mean(abs(back(pr$pred) - back(pr$obs))))
}
nf <- median(d$y_sd, na.rm = TRUE)
cvperf <- rbind(perf(model, "spatial_blocks"), perf(model_rand, "random"))
cvperf$noise_floor_withincell_sd_log <- nf
write.csv(cvperf, file.path(DIRS$tables, "cv_performance_toc.csv"), row.names = FALSE)
print(cvperf, digits = 3)
vi <- varImp(model, scale = FALSE)$importance
write.csv(data.frame(variable = rownames(vi), importance = vi[, 1]), file.path(DIRS$tables, "variable_importance_toc.csv"), row.names = FALSE)
write.csv(d, file.path(DIRS$data_iter, "toc_training_table.csv"), row.names = FALSE)

# --- geographic distance diagnostic ------------------------------------------------------------------
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))
Xt <- to_xyz(d$x_lon, d$y_lat)
s2s <- chord_km(get.knn(Xt, k = 1)$nn.dist[, 1])
cvd <- unlist(lapply(seq_along(folds$indexOut), function(i)
  chord_km(get.knnx(Xt[folds$index[[i]], , drop = FALSE], Xt[folds$indexOut[[i]], , drop = FALSE], k = 1)$nn.dist[, 1])))
tmpl <- stk[["sst_mean"]]; oc <- which(!is.na(values(tmpl)[, 1]))
set.seed(CFG$seed); smp <- sample(oc, 20000, prob = values(cellSize(tmpl, unit = "km"))[oc, 1])
pxy <- xyFromCell(tmpl, smp)
p2s <- chord_km(get.knnx(Xt, to_xyz(pxy[, 1], pxy[, 2]), k = 1)$nn.dist[, 1])
gd <- rbind(data.frame(what = "sample-to-sample", dist_km = s2s), data.frame(what = "CV test-to-train", dist_km = cvd),
            data.frame(what = "prediction-to-sample (global ocean)", dist_km = p2s))
write.csv(gd, file.path(DIRS$tables, "geodist_toc.csv"), row.names = FALSE)
print(aggregate(dist_km ~ what, gd, function(z) c(median = median(z), q90 = unname(quantile(z, 0.9)))))
