# =============================================================================
# 08w_error_structure.R — where does the out-of-sample error have structure, and at which level of ignorance?
# =============================================================================
# Protocol in Ignorance_MS/error_structure_reasoning.md. Out-of-sample residuals of the spatial-block CV (1000-km
# blocks, caret savePredictions) of the clean chain (iter1c_clean, iter2c_clean) are tested against declared partitions:
#   region (8, k-means)       space, not in the model
#   latitude band (30 deg)    space, not in the model (correlated with temperature)
#   longitude sector (60 deg) space, not in the model
#   depth zone                environment; depth IS a predictor -> structure means form/interaction (L4)
#   lithology class           candidate variable: absent in iteration 1 (L2 named), present in iteration 2 (L4)
# Responses: bias = mean signed residual; magnitude = mean absolute residual. Replicate units = 300-km blocks within
# partition levels (>= 5 cells). Effect size eta2 = between-level share of variance among unit means; p from 2000
# permutations of level labels among units; Benjamini-Hochberg within iteration. Structured = eta2 >= 0.01 and q < 0.05.
# Out-of-sample gain (partitions not in the model): refit with the partition as a predictor on the same folds
# (ranger, 300 trees, tuned mtry/min.node.size), paired 300-km block bootstrap of the MSE difference.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(ggplot2) })
ITS <- c("iter1c_clean", "iter2c_clean"); NPERM <- 2000
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
tex <- rast(path_texture_stack())[[paste0("litho_t", 1:6)]]
eta2 <- function(y, g) { m <- tapply(y, g, mean); n <- tapply(y, g, length); sum(n * (m - mean(y))^2) / sum((y - mean(y))^2) }

tests <- list(); gains <- list(); lvl <- list()
for (it in ITS) {
  m <- readRDS(file.path(ROOT, "outputs", it, "models", "rf_spatialCV_toc.rds")); bt <- m$bestTune
  vars <- setdiff(names(m$trainingData), c(".outcome", ".weights"))
  d <- read.csv(file.path(ROOT, "outputs", it, "data", "toc_training_table.csv"))
  stopifnot(nrow(m$pred) == nrow(d))
  d$pred_cv <- NA_real_; d$pred_cv[m$pred$rowIndex] <- m$pred$pred
  d$e <- d$pred_cv - d$y; d$ae <- abs(d$e)
  set.seed(CFG$seed); d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster
  exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
  d$blk <- paste(floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
  lv <- as.matrix(extract(tex, cbind(d$x_lon, d$y_lat))); ls <- rowSums(lv)
  P <- list(
    region           = list(g = RN[d$region], type = "space", in_model = FALSE),
    latitude_band    = list(g = as.character(cut(d$y_lat, seq(-90, 90, 30), include.lowest = TRUE)), type = "space", in_model = FALSE),
    longitude_sector = list(g = as.character(cut(d$x_lon, seq(-180, 180, 60), include.lowest = TRUE)), type = "space", in_model = FALSE),
    depth_zone       = list(g = as.character(cut(d$depth, c(-Inf, 200, 1000, 3000, Inf), labels = c("shelf <=200 m", "slope 200-1000 m", "1000-3000 m", "abyss >3000 m"))),
                            type = "environment", in_model = "depth" %in% vars),
    lithology_class  = list(g = ifelse(is.na(ls) | ls < 0.5, "unclassified", paste0("type ", max.col(lv, ties.method = "first"))),
                            type = "candidate variable", in_model = any(LITHO_VARS %in% vars)))
  rows <- list()
  for (pn in names(P)) {
    g <- P[[pn]]$g
    U <- aggregate(cbind(bias = d$e, magnitude = d$ae) ~ blk + level, data = data.frame(d[, c("e", "ae", "blk")], level = g), FUN = mean)
    cnt <- table(paste(d$blk, g, sep = "|")); U$n <- as.numeric(cnt[paste(U$blk, U$level, sep = "|")])
    U <- U[U$n >= 5, ]
    for (resp in c("bias", "magnitude")) {
      obs <- eta2(U[[resp]], U$level); set.seed(CFG$seed)
      perm <- replicate(NPERM, eta2(U[[resp]], sample(U$level)))
      rows[[length(rows) + 1]] <- data.frame(iteration = it, partition = pn, type = P[[pn]]$type, in_model = P[[pn]]$in_model, response = resp,
        n_units = nrow(U), n_levels = length(unique(U$level)), eta2 = obs, f2 = obs / (1 - obs), p_perm = (1 + sum(perm >= obs)) / (NPERM + 1))
      lm_ <- tapply(U[[resp]], U$level, mean)
      lvl[[length(lvl) + 1]] <- data.frame(iteration = it, partition = pn, response = resp, level = names(lm_), unit_mean = as.numeric(lm_),
                                           n_units = as.numeric(table(U$level)[names(lm_)]))
    }
  }
  R <- do.call(rbind, rows); R$q_BH <- p.adjust(R$p_perm, "BH"); tests[[it]] <- R
  msg("%s: structure tests done", it)

  # out-of-sample gain from adding each partition not in the model (same folds)
  # two designs: spatial blocks (interpolation-like) and leave-region-out (transfer, where the lithology gain appeared);
  # partitions enter as 0/1 dummies so that a level absent from the training folds is handled
  fitcv <- function(X, fold) { p <- rep(NA_real_, nrow(d))
    for (k in sort(unique(fold))) { tr <- fold != k
      mod <- ranger(x = X[tr, , drop = FALSE], y = d$y[tr], num.trees = 300, mtry = min(bt$mtry, ncol(X)), min.node.size = bt$min.node.size,
                    splitrule = "variance", num.threads = CFG$cores, seed = CFG$seed)
      p[!tr] <- predict(mod, X[!tr, , drop = FALSE], num.threads = CFG$cores)$predictions }
    p }
  bl <- split(seq_len(nrow(d)), d$blk); set.seed(CFG$seed)
  IDX <- replicate(1000, unlist(bl[sample(length(bl), replace = TRUE)], use.names = FALSE), simplify = FALSE)
  for (design in c("spatial_blocks_1000km", "leave_region_out")) {
    fold <- if (design == "spatial_blocks_1000km") d$fold else d$region
    e0 <- fitcv(d[, vars], fold) - d$y
    for (pn in names(P)[!vapply(P, function(z) isTRUE(z$in_model), TRUE)]) {
      if (design == "leave_region_out" && pn == "region") next   # a withheld region has no level to learn from
      Dm <- model.matrix(~ g - 1, data.frame(g = factor(P[[pn]]$g))); colnames(Dm) <- paste0("part_", seq_len(ncol(Dm)))
      e1 <- fitcv(cbind(d[, vars], Dm), fold) - d$y
      dm <- mean(e1^2) - mean(e0^2); bs <- vapply(IDX, function(i) mean(e1[i]^2) - mean(e0[i]^2), 1)
      gains[[length(gains) + 1]] <- data.frame(iteration = it, design = design, partition = pn, rmse_without = sqrt(mean(e0^2)), rmse_with = sqrt(mean(e1^2)),
        mse_pct_change = 100 * dm / mean(e0^2), ci_lo_pct = 100 * quantile(bs, 0.025) / mean(e0^2), ci_hi_pct = 100 * quantile(bs, 0.975) / mean(e0^2))
      msg("%s [%s] + %s: RMSE %.4f -> %.4f", it, design, pn, sqrt(mean(e0^2)), sqrt(mean(e1^2)))
    }
  }
}
T <- do.call(rbind, tests); Gn <- do.call(rbind, gains); L <- do.call(rbind, lvl)
T$structured <- T$eta2 >= 0.01 & T$q_BH < 0.05
T$level_of_ignorance <- ifelse(!T$structured, "no structure", ifelse(T$in_model, "L4 form or interaction",
                         ifelse(T$type == "candidate variable", "L2 named driver", "L2 unnamed (space)")))
write.csv(T, file.path(DIRS$compare, "error_structure_tests.csv"), row.names = FALSE)
write.csv(Gn, file.path(DIRS$compare, "error_structure_gain.csv"), row.names = FALSE)
write.csv(L, file.path(DIRS$compare, "error_structure_level_means.csv"), row.names = FALSE)
print(T, digits = 3, row.names = FALSE); print(Gn, digits = 3, row.names = FALSE)

T$lab <- factor(T$partition, levels = rev(names(P)))
T$it_lab <- ifelse(T$iteration == "iter1c_clean", "iteration 1 (base)", "iteration 2 (+ lithology)")
fig <- ggplot(T, aes(eta2, lab, colour = level_of_ignorance, shape = response)) + geom_vline(xintercept = 0.01, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) + facet_wrap(~it_lab) + scale_x_sqrt() +
  scale_colour_manual(values = c(`no structure` = "grey60", `L4 form or interaction` = "#7F77DD", `L2 named driver` = "#C2410C", `L2 unnamed (space)` = "#1D6A73"), name = NULL) +
  labs(x = "eta2: share of variance among block means explained by the partition (sqrt scale; dashed = 0.01)", y = NULL,
       title = "Seafloor TOC: structure of the out-of-sample error", subtitle = "bias = signed residual, magnitude = absolute residual; structured = eta2 >= 0.01 and BH-corrected permutation p < 0.05") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(DIRS$compare, "fig_error_structure_toc.png"), fig, width = 200, height = 110, units = "mm", dpi = 250, bg = "white")
msg("error structure (TOC) done")
