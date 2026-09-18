# =============================================================================
# 07b_diagnose.R — why did the iteration-1 ignorance map fail to locate error?
# =============================================================================
# Run on iteration 1 outputs BEFORE texture enters the model. Candidate causes, each
# mapped to an ignorance level of the manuscript:
#   H1 extrapolation / sampling gap (L3): error should follow DI and distance to data.
#   H2 irreducible noise (L5): error should approach the within-cell noise floor.
#   H3 missing essential variable (L2, conceptual model incomplete): residuals should be
#      explained by a variable absent from the predictors (candidate: sediment texture).
#   H4 local heterogeneity at sub-grid scale (L5 at this resolution, or L3 if resolvable):
#      error should follow the spatial variability of TOC itself.
# Tests use the leave-region-out residuals of iteration 1 (05_transfer.R):
#   (a) Spearman of |error| with DI, distance, texture variables and local heterogeneity;
#   (b) out-of-region explained variance of the signed residual by a random forest on
#       texture only vs baseline predictors only (leave-region-out, so no leakage across regions);
#   (c) the regions where the map failed most: texture of withheld vs training data.

Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1_baseline"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(FNN); library(ggplot2); library(patchwork) })
set.seed(CFG$seed)

tt  <- read.csv(file.path(DIRS$tables, "transfer_points_toc.csv"))
cal <- read.csv(file.path(DIRS$tables, "calibration_per_region_toc.csv"))
stk <- rast(path_texture_stack())
X   <- stk[tt$cell]
tt  <- cbind(tt, X[, c(TEXTURE_VARS, CANDIDATES_BASE)])
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))

# local heterogeneity: SD of observed log10 TOC among the 10 nearest observations (within the same region)
tt$local_sd <- NA_real_
for (r in unique(tt$region)) {
  i <- which(tt$region == r); P <- to_xyz(tt$lon[i], tt$lat[i])
  nn <- get.knn(P, k = min(10, length(i) - 1))$nn.index
  tt$local_sd[i] <- vapply(seq_along(i), function(j) sd(tt$obs[i][c(j, nn[j, ])]), 0)
}
noise_floor <- median(read.csv(file.path(DIRS$data, "toc_cells.csv"))$y_sd, na.rm = TRUE)

# (a) correlations
vars_a <- c("DInorm", "dist_km_proxy", "local_sd", TEXTURE_VARS)
tt$dist_km_proxy <- NA_real_
d_train <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
for (r in unique(tt$region)) {
  i <- which(tt$region == r); j <- which(!d_train$cell %in% tt$cell[i])
  tt$dist_km_proxy[i] <- chord_km(get.knnx(to_xyz(d_train$x_lon[j], d_train$y_lat[j]), to_xyz(tt$lon[i], tt$lat[i]), k = 1)$nn.dist[, 1])
}
cor_tab <- do.call(rbind, lapply(vars_a, function(v) {
  ok <- is.finite(tt[[v]]) & is.finite(tt$abs_err)
  c1 <- suppressWarnings(cor.test(tt[[v]][ok], tt$abs_err[ok], method = "spearman", exact = FALSE))
  c2 <- suppressWarnings(cor.test(tt[[v]][ok], tt$err[ok], method = "spearman", exact = FALSE))
  data.frame(variable = v, rho_abs_error = unname(c1$estimate), p_abs = c1$p.value, rho_signed_error = unname(c2$estimate), p_signed = c2$p.value)
}))
cor_tab$hypothesis <- c(DInorm = "H1", dist_km_proxy = "H1", local_sd = "H4", setNames(rep("H3", length(TEXTURE_VARS)), TEXTURE_VARS))[cor_tab$variable]
write.csv(cor_tab, file.path(DIRS$tables, "diagnosis_correlations_toc.csv"), row.names = FALSE)
print(cor_tab, digits = 3)

# (b) out-of-region explained variance of residuals
oor_r2 <- function(vars) {
  pred <- rep(NA_real_, nrow(tt))
  for (r in unique(tt$region)) {
    tr <- tt$region != r & complete.cases(tt[, vars]); te <- tt$region == r & complete.cases(tt[, vars])
    m <- ranger(x = tt[tr, vars, drop = FALSE], y = tt$err[tr], num.trees = 300, num.threads = CFG$cores, seed = CFG$seed)
    pred[te] <- predict(m, tt[te, vars, drop = FALSE])$predictions
  }
  ok <- !is.na(pred)
  1 - sum((tt$err[ok] - pred[ok])^2) / sum((tt$err[ok] - mean(tt$err[ok]))^2)
}
res_tab <- data.frame(predictor_set = c("baseline predictors (already in the model)", "texture only", "texture + baseline"),
                      out_of_region_R2_of_residual = c(oor_r2(CANDIDATES_BASE), oor_r2(TEXTURE_VARS), oor_r2(c(CANDIDATES_BASE, TEXTURE_VARS))))
write.csv(res_tab, file.path(DIRS$tables, "diagnosis_residual_R2_toc.csv"), row.names = FALSE)
print(res_tab, digits = 3)

# (c) worst-calibrated regions: texture of withheld vs training cells
worst <- cal$region[order(-abs(log(cal$ratio_best)))][1:2]
reg_tab <- do.call(rbind, lapply(worst, function(r) {
  w <- tt[tt$region == r, ]; o <- tt[tt$region != r, ]
  do.call(rbind, lapply(TEXTURE_VARS, function(v) data.frame(region = r, ratio_realised_expected = cal$ratio_best[cal$region == r],
    variable = v, withheld_median = median(w[[v]], na.rm = TRUE), others_median = median(o[[v]], na.rm = TRUE),
    rho_abs_error_within_region = suppressWarnings(cor(w[[v]], w$abs_err, method = "spearman", use = "complete.obs")))))
}))
write.csv(reg_tab, file.path(DIRS$tables, "diagnosis_worst_regions_toc.csv"), row.names = FALSE)
print(reg_tab, digits = 3)

# figure
bin_plot <- function(v, lab) {
  b <- tt[is.finite(tt[[v]]), ]; b$bin <- cut(rank(b[[v]], ties.method = "first"), 10, labels = FALSE)
  s <- do.call(rbind, lapply(split(b, b$bin), function(z) data.frame(x = median(z[[v]]), rmse = sqrt(mean(z$err^2)), bias = mean(z$err))))
  ggplot(s, aes(x)) + geom_hline(yintercept = 0, colour = "grey70") + geom_hline(yintercept = noise_floor, linetype = 3, colour = "grey40") +
    geom_line(aes(y = rmse, colour = "RMSE")) + geom_point(aes(y = rmse, colour = "RMSE")) +
    geom_line(aes(y = bias, colour = "mean error (bias)")) + geom_point(aes(y = bias, colour = "mean error (bias)")) +
    scale_colour_manual(values = c("RMSE" = "#D7301F", "mean error (bias)" = "#3E5C76"), name = NULL) +
    labs(x = lab, y = "log10 TOC") + theme_classic(base_size = 8) + theme(legend.position = "bottom")
}
g <- (bin_plot("DInorm", "H1: DI / AOA threshold") | bin_plot("local_sd", "H4: local SD of observed TOC")) /
     (bin_plot("gs_d50", "H3: log10 grain size D50 (mm)") | bin_plot("porosity", "H3: porosity (%)")) +
  plot_annotation(title = "Iteration 1 diagnosis: error of withheld regions against candidate causes (deciles; dotted = noise floor)",
                  subtitle = sprintf("Out-of-region R2 of residuals: baseline %.2f | texture %.2f | texture + baseline %.2f",
                                     res_tab$out_of_region_R2_of_residual[1], res_tab$out_of_region_R2_of_residual[2], res_tab$out_of_region_R2_of_residual[3]))
ggsave(file.path(DIRS$figs, "toc_diagnosis_iter1.png"), g, width = 200, height = 160, units = "mm", dpi = 250, bg = "white")
msg("diagnosis written to %s", DIRS$tables)
