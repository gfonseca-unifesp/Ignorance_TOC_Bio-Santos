# =============================================================================
# 08s_wedge_neutralise.R — the SE Pacific wedge predates supply and texture: which base predictor carries it?
# =============================================================================
# 08r: smoothing any single predictor at 2 deg leaves the wedge (it spans ~20 deg, so smoothing only blurs its edges).
# The wedge is already in the preliminary model (11 base predictors). Here each of them is held at its median over the
# box (all other predictors unchanged) and the box is re-predicted with the preliminary model (500 trees, all data).
# The predictor whose neutralisation erases the wedge carries it; pairs of bottom-water layers are tested as well.
Sys.setenv(DEMO3B_ITER = "iter1_baseline")
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(ggplot2); library(patchwork); library(tidyterra) })

m0 <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds")); bt <- m0$bestTune
V  <- setdiff(names(m0$trainingData), c(".outcome", ".weights"))
d  <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
mod <- ranger(x = d[, V], y = d$y, num.trees = 500, mtry = bt$mtry, min.node.size = bt$min.node.size, splitrule = "variance",
              num.threads = CFG$cores, seed = CFG$seed)
eB <- ext(-140, -90, -50, 0)
S  <- crop(rast(path_base_toc_stack())[[V]], eB)
X0 <- as.data.frame(S, na.rm = FALSE); ok <- complete.cases(X0)
pr <- function(X) { p <- rep(NA_real_, nrow(X)); p[ok] <- predict(mod, X[ok, V], num.threads = CFG$cores)$predictions; r <- S[[1]]; values(r) <- p; names(r) <- "pred"; r }
hold <- function(ks) { X <- X0; for (k in ks) X[[k]][ok] <- median(X0[[k]][ok]); pr(X) }
R <- c(list(`none (preliminary model)` = pr(X0)), setNames(lapply(V, hold), paste(V, "held")),
       list(`sbt + sbs held` = hold(c("sbt_mean", "sbs_mean")), `sbt + sbs + o2b held` = hold(c("sbt_mean", "sbs_mean", "o2b_mean")),
            `all bottom-water layers held` = hold(intersect(c("sbt_mean", "sbs_mean", "o2b_mean", "phyc_bot", "sws_bot"), V))))
rng <- as.numeric(quantile(values(R[[1]]), c(0.02, 0.98), na.rm = TRUE))
pz <- function(r, title) ggplot() + geom_spatraster(data = r) +
  scale_fill_viridis_c(limits = rng, oob = scales::squish, na.value = "grey85", name = "log10 TOC") + coord_sf(expand = FALSE) +
  labs(title = title) + theme_minimal(base_size = 6) + theme(axis.text = element_blank(), plot.title = element_text(size = 6.5, face = "bold"), legend.key.width = unit(2, "mm"))
fig <- wrap_plots(Map(pz, R, names(R)), ncol = 5, guides = "collect")
ggsave(file.path(DIRS$compare, "fig_wedge_neutralise.png"), fig, width = 240, height = 170, units = "mm", dpi = 220, bg = "white")
msg("wedge neutralisation done: %s", paste(V, collapse = ", "))
