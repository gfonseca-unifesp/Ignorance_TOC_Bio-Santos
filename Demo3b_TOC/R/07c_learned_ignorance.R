# =============================================================================
# 07c_learned_ignorance.R — H6: can a LEARNED error model locate error where the
# DI-based ignorance map cannot?
# =============================================================================
# DI-based error models (C0-C3, 07_calibration.R) assume that error follows dissimilarity.
# Here the error model is learned: a random forest predicts squared error from
#   - the environmental predictors of the TOC model (where does the model fail?),
#   - DI, distance to the nearest observation, predicted TOC level,
#   - local sampling density (observations within 500 km) and local heterogeneity
#     (SD of observed TOC among the 10 nearest observations).
# Same nested leave-one-region-out design as 07: for each withheld region, the TOC model, DI and
# the error model are built from the other regions only (inner block-CV and leave-region-out
# errors). Covariates of calibration rows use only their own CV-training data, as at prediction.
# Evaluation on the same withheld points as C0-C3: calibration metrics and localisation
# (Spearman) at point, 2, 5 and 10 deg cells and region scale, with bootstrap CIs.
source("R/00_config.R")
suppressPackageStartupMessages({ library(caret); library(ranger); library(CAST); library(FNN) })
set.seed(CFG$seed)

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
bt    <- model$bestTune
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz   <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster   # same as 05 / 07

fit_rf <- function(df, imp = "none") ranger(x = df[, vars, drop = FALSE], y = df$y, num.trees = 300,
  mtry = min(bt$mtry, length(vars)), min.node.size = bt$min.node.size, splitrule = "variance",
  importance = imp, num.threads = CFG$cores, seed = CFG$seed)

# local context of target rows, given only the observations available to the model that predicts them
context <- function(src, tgt) {
  k  <- min(200, nrow(src))
  nn <- get.knnx(to_xyz(src$x_lon, src$y_lat), to_xyz(tgt$x_lon, tgt$y_lat), k = k)
  dk <- nn$nn.dist; dk[] <- chord_km(as.vector(dk))   # keep matrix dimensions (pmin drops them)
  data.frame(dist_km = dk[, 1], ldens = log1p(rowSums(dk < 500)),
             lsd10 = apply(nn$nn.index[, 1:min(10, k), drop = FALSE], 1, function(i) sd(src$y[i])))
}

cv_errors <- function(df, fold_id, w) {
  ks <- sort(unique(fold_id))
  CVtest <- lapply(ks, function(k) which(fold_id == k)); CVtrain <- lapply(ks, function(k) which(fold_id != k))
  tdi <- suppressMessages(trainDI(train = df[, vars, drop = FALSE], variables = vars, weight = w,
                                  CVtest = CVtest, CVtrain = CVtrain, verbose = FALSE))
  do.call(rbind, lapply(seq_along(ks), function(i) {
    tr <- CVtrain[[i]]; te <- CVtest[[i]]
    p <- predict(fit_rf(df[tr, ]), df[te, vars, drop = FALSE], num.threads = CFG$cores)$predictions
    cbind(df[te, vars, drop = FALSE], data.frame(e2 = (p - df$y[te])^2, DI = tdi$trainDI[te], pred = p), context(df[tr, ], df[te, ]))
  }))
}

feats <- c(vars, "DI", "ldist", "pred", "ldens", "lsd10")
f_outer <- file.path(DIRS$models, "learned_ignorance_outer_toc.rds")
if (!file.exists(f_outer)) {
  outer <- list(); imp <- list()
  for (r in sort(unique(d$region))) {
    t0 <- Sys.time()
    trn <- d[d$region != r, ]; tst <- d[d$region == r, ]
    m_out <- fit_rf(trn, "permutation")
    w <- as.data.frame(t(pmax(m_out$variable.importance[vars], 0)))
    ks <- sort(unique(trn$fold))
    tdi_o <- suppressMessages(trainDI(train = trn[, vars], variables = vars, weight = w,
                                      CVtest = lapply(ks, function(k) which(trn$fold == k)),
                                      CVtrain = lapply(ks, function(k) which(trn$fold != k)), verbose = FALSE))
    a <- suppressMessages(aoa(tst[, vars], trainDI = tdi_o, verbose = FALSE))
    p <- predict(m_out, tst[, vars], num.threads = CFG$cores)$predictions
    cal <- rbind(cv_errors(trn, trn$fold, w), cv_errors(trn, trn$region, w))
    cal$ldist <- log10(cal$dist_km + 1)
    em <- ranger(x = cal[, feats], y = cal$e2, num.trees = 500, min.node.size = 100, importance = "permutation",
                 num.threads = CFG$cores, seed = CFG$seed)
    nd <- cbind(tst[, vars, drop = FALSE], data.frame(DI = a$DI, pred = p), context(trn, tst)); nd$ldist <- log10(nd$dist_km + 1)
    res <- data.frame(region = r, lon = tst$x_lon, lat = tst$y_lat, obs = tst$y, pred = p, e2 = (p - tst$y)^2,
                      exp_L = sqrt(pmax(predict(em, nd[, feats], num.threads = CFG$cores)$predictions, 1e-6)))
    outer[[length(outer) + 1]] <- res
    imp[[length(imp) + 1]] <- data.frame(region = r, feature = names(em$variable.importance), importance = unname(em$variable.importance))
    msg("outer region %d done (%.1f min): realised RMSE %.3f | learned expected %.3f | rho points %.2f", r,
        as.numeric(difftime(Sys.time(), t0, units = "mins")), sqrt(mean(res$e2)), sqrt(mean(res$exp_L^2)),
        cor(res$exp_L, sqrt(res$e2), method = "spearman"))
  }
  saveRDS(list(outer = outer, imp = imp), f_outer)
}
L <- readRDS(f_outer)
ev <- do.call(rbind, L$outer)

# same withheld points, DI-based models from 07
key <- function(lon, lat) paste(round(lon, 3), round(lat, 3))
c07 <- read.csv(file.path(DIRS$tables, "calibration_outer_points_toc.csv"))
ev <- cbind(ev, c07[match(key(ev$lon, ev$lat), key(c07$lon, c07$lat)), c("exp_C0", "exp_C2", "exp_C3")])
msg("points matched to 07: %d of %d", sum(!is.na(ev$exp_C2)), nrow(ev))
write.csv(ev, file.path(DIRS$tables, "learned_ignorance_outer_points_toc.csv"), row.names = FALSE)

evaluate <- function(ex, e2, label) {
  b <- cut(rank(ex, ties.method = "first"), 10, labels = FALSE)
  bb <- data.frame(expected = tapply(ex^2, b, function(z) sqrt(mean(z))), realised = tapply(e2, b, function(z) sqrt(mean(z))))
  lmb <- lm(realised ~ expected, bb)
  data.frame(method = label, slope = unname(coef(lmb)[2]), MACE = mean(abs(bb$realised - bb$expected)),
             ratio_realised_expected = sqrt(mean(e2)) / sqrt(mean(ex^2)),
             spearman_bins = cor(bb$expected, bb$realised, method = "spearman"),
             spearman_points = cor(ex, sqrt(e2), method = "spearman"), coverage90 = mean(sqrt(e2) <= 1.645 * ex))
}
M <- c(L = "exp_L", C0 = "exp_C0", C2 = "exp_C2", C3 = "exp_C3")
ok <- complete.cases(ev[, M])
metrics <- do.call(rbind, lapply(names(M), function(m) evaluate(ev[[M[[m]]]][ok], ev$e2[ok], m)))
write.csv(metrics, file.path(DIRS$tables, "learned_ignorance_metrics_toc.csv"), row.names = FALSE)
print(metrics, digits = 3)

set.seed(CFG$seed)
boot_rho <- function(a, b, B = 500) { n <- length(a)
  quantile(replicate(B, { i <- sample(n, n, replace = TRUE); suppressWarnings(cor(a[i], b[i], method = "spearman")) }), c(0.025, 0.975), na.rm = TRUE) }
e <- ev[ok, ]
sc <- do.call(rbind, lapply(names(M), function(m) { ex <- e[[M[[m]]]]
  one <- function(scale, real, expd) { ci <- boot_rho(real, expd)
    data.frame(method = m, scale = scale, n_units = length(real), spearman = cor(real, expd, method = "spearman"), ci_lo = ci[[1]], ci_hi = ci[[2]]) }
  out <- list(one("point", sqrt(e$e2), ex))
  for (s in c(2, 5, 10)) { g <- paste(e$region, floor(e$lon / s), floor(e$lat / s)); k <- table(g); keep <- names(k)[k >= 10]
    out[[length(out) + 1]] <- one(sprintf("%g deg", s), sqrt(tapply(e$e2, g, mean))[keep], sqrt(tapply(ex^2, g, mean))[keep]) }
  out[[length(out) + 1]] <- one("region (8)", sqrt(tapply(e$e2, e$region, mean)), sqrt(tapply(ex^2, e$region, mean)))
  do.call(rbind, out) }))
write.csv(sc, file.path(DIRS$tables, "learned_ignorance_by_scale_toc.csv"), row.names = FALSE)
sc$txt <- sprintf("%5.2f [%5.2f,%5.2f]", sc$spearman, sc$ci_lo, sc$ci_hi)
print(reshape(sc[, c("method", "scale", "txt")], idvar = "scale", timevar = "method", direction = "wide"), row.names = FALSE)

imp <- do.call(rbind, L$imp)
imp_m <- aggregate(importance ~ feature, imp, mean); imp_m <- imp_m[order(-imp_m$importance), ]
imp_m$share <- 100 * pmax(imp_m$importance, 0) / sum(pmax(imp_m$importance, 0))
write.csv(imp_m, file.path(DIRS$tables, "learned_ignorance_importance_toc.csv"), row.names = FALSE)
print(head(imp_m, 12), digits = 3, row.names = FALSE)
msg("learned ignorance done")
