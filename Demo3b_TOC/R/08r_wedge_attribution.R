# =============================================================================
# 08r_wedge_attribution.R — which predictor draws the SE Pacific wedge into the iteration-2 prediction?
# =============================================================================
# 08q showed that smoothing the Bio-ORACLE bottom layers at 0.5 deg removes the sharp seam steps without loss of skill,
# but the wedge of low TOC remains. Attribution: predict the SE Pacific box with the full iteration-2 model (500 trees,
# all data) after replacing one predictor at a time by its 2-degree moving mean. The wedge is a mesoscale structure, so
# its source is the predictor whose smoothing removes most of the prediction's 0.5-2 deg structure
# (mean |prediction - 2-deg mean of prediction|), and whose smoothed panel no longer shows the wedge.
Sys.setenv(DEMO3B_ITER = "iter2b_supply_texture")
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(ggplot2); library(patchwork); library(tidyterra) })

m2 <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds")); bt <- m2$bestTune
V  <- setdiff(names(m2$trainingData), c(".outcome", ".weights"))
d  <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
mod <- ranger(x = d[, V], y = d$y, num.trees = 500, mtry = bt$mtry, min.node.size = bt$min.node.size, splitrule = "variance",
              num.threads = CFG$cores, seed = CFG$seed)
eB <- ext(-140, -90, -50, 0); eE <- ext(-143, -87, -53, 3)
S  <- crop(rast(path_toc_stack())[[V]], eE)
S2 <- rast(lapply(V, function(k) focal(S[[k]], w = 21, fun = "mean", na.rm = TRUE))); names(S2) <- V
X0 <- as.data.frame(S, na.rm = FALSE); X2 <- as.data.frame(S2, na.rm = FALSE)
for (k in V) { na <- is.na(X2[[k]]); X2[[k]][na] <- X0[[k]][na] }
ok <- complete.cases(X0)
inB <- rep(FALSE, ncell(S)); inB[cells(S, eB)] <- TRUE
pr <- function(X) { p <- rep(NA_real_, nrow(X)); p[ok] <- predict(mod, X[ok, V], num.threads = CFG$cores)$predictions; r <- S[[1]]; values(r) <- p; names(r) <- "pred"; r }
hf <- function(r) { h <- abs(r - focal(r, w = 21, fun = "mean", na.rm = TRUE)); mean(values(h)[inB, 1], na.rm = TRUE) }
r0 <- pr(X0); h0 <- hf(r0)
res <- list(); R <- list(full = r0)
for (k in V) {
  X <- X0; X[[k]] <- X2[[k]]; rk <- pr(X); R[[k]] <- rk
  res[[k]] <- data.frame(variable = k, mesoscale_structure = hf(rk), pct_removed = 100 * (1 - hf(rk) / h0),
                         mean_abs_change_pred = mean(abs(values(rk)[inB, 1] - values(r0)[inB, 1]), na.rm = TRUE))
}
res <- do.call(rbind, res); res <- res[order(-res$pct_removed), ]
write.csv(res, file.path(DIRS$compare, "wedge_attribution.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
top <- head(res$variable, 5)
pz <- function(r, title) ggplot() + geom_spatraster(data = crop(r, eB)) +
  scale_fill_viridis_c(limits = c(-1.1, 0.3), oob = scales::squish, na.value = "grey85", name = "log10 TOC") + coord_sf(expand = FALSE) +
  labs(title = title) + theme_minimal(base_size = 6.5) + theme(axis.text = element_blank(), plot.title = element_text(size = 7, face = "bold"), legend.key.width = unit(2, "mm"))
fig <- wrap_plots(c(list(pz(r0, "full iteration-2 prediction")), lapply(top, function(k) pz(R[[k]], sprintf("%s smoothed (2 deg): -%.0f%%", k, res$pct_removed[res$variable == k])))),
                  nrow = 2, guides = "collect")
ggsave(file.path(DIRS$compare, "fig_wedge_attribution.png"), fig, width = 200, height = 140, units = "mm", dpi = 250, bg = "white")
msg("wedge attribution done")
