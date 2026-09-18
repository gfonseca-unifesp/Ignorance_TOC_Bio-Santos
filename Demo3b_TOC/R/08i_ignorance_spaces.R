# =============================================================================
# 08i_ignorance_spaces.R — the same ignorance in two spaces: geographic and environmental
# =============================================================================
# H7 showed that data chosen for environmental novelty (DI logic) or model uncertainty reduce error
# about twice as fast as random, while geographic gap filling does not beat random. The ignorance map is
# therefore delivered in more than one space, each answering a different question:
#   (a) geographic space   — where on the seafloor are conditions least represented by the data?
#   (b) environmental space — which combinations of conditions are unsampled? (importance-weighted,
#       standardised predictors, as in the DI; first two principal components; ocean cells coloured by
#       mean DI / AOA threshold, observation cells overlaid)
#   (c) the test that decides which space guides new data (H7, 08h_targeted_sampling.R)
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1_baseline"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(tidyterra); library(rnaturalearth) })
set.seed(CFG$seed)

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
summ  <- read.csv(file.path(DIRS$tables, "toc_summary.csv"))
thr   <- as.numeric(summ$value[summ$metric == "AOA_DI_threshold"])
stk   <- rast(path_toc_stack())
DIn   <- rast(file.path(DIRS$rasters, "toc_DI.tif")) / thr; names(DIn) <- "DInorm"

# --- (b) environmental space ------------------------------------------------------------------------------------
imp <- pmax(model$finalModel$variable.importance[vars], 0); imp <- imp / sum(imp)
oc  <- which(!is.na(values(DIn)[, 1]))
area <- values(cellSize(DIn, unit = "km"))[oc, 1]
smp <- oc[order(rexp(length(oc)) / area)[1:150000]]   # area-weighted sampling without replacement (exponential keys; fast)
X   <- stk[[vars]][smp]; X <- as.matrix(X); ok <- complete.cases(X); X <- X[ok, ]; smp <- smp[ok]
mu  <- colMeans(d[, vars]); sdv <- apply(d[, vars], 2, sd)
wz  <- function(M) sweep(sweep(sweep(as.matrix(M), 2, mu, "-"), 2, sdv, "/"), 2, imp, "*")
Zo  <- wz(X); Zt <- wz(d[, vars])
pc  <- prcomp(Zo, center = TRUE, scale. = FALSE)
ve  <- round(100 * pc$sdev[1:2]^2 / sum(pc$sdev^2))
env <- data.frame(predict(pc, Zo)[, 1:2], DI = values(DIn)[smp, 1])
trn <- data.frame(predict(pc, Zt)[, 1:2])
ld  <- data.frame(var = rownames(pc$rotation), PC1 = pc$rotation[, 1], PC2 = pc$rotation[, 2])
ld  <- ld[order(-(ld$PC1^2 + ld$PC2^2)), ][1:6, ]
sc  <- 0.8 * max(abs(range(env$PC1, env$PC2))) / max(sqrt(ld$PC1^2 + ld$PC2^2))
pb <- ggplot(env, aes(PC1, PC2)) +
  stat_summary_2d(aes(z = pmin(DI, 2)), fun = mean, bins = 80) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 2), name = "mean DI / AOA threshold\n(ocean cells)") +
  geom_point(data = trn, colour = "#00B4D8", size = 0.12, alpha = 0.25) +
  geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1 * sc, yend = PC2 * sc), inherit.aes = FALSE,
               arrow = arrow(length = unit(1.2, "mm")), colour = "grey20", linewidth = 0.3) +
  geom_text(data = ld, aes(x = PC1 * sc * 1.12, y = PC2 * sc * 1.12, label = var), inherit.aes = FALSE, size = 2.3, colour = "grey10") +
  labs(x = sprintf("environmental axis 1 (%d%%)", ve[1]), y = sprintf("environmental axis 2 (%d%%)", ve[2]),
       title = "b  Environmental space: which conditions are unsampled?",
       subtitle = "ocean cells (area-weighted sample) coloured by dissimilarity to the data; cyan = observation cells") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "right")

# --- (a) geographic space ----------------------------------------------------------------------------------------
PROJ <- "+proj=eqearth +datum=WGS84 +units=m"
tp <- project(DIn, PROJ, res = 25000); dp <- project(DIn, tp, method = "bilinear")
land <- st_transform(ne_countries(scale = 50, returnclass = "sf"), PROJ)
obs_sf <- st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), PROJ)
pa <- ggplot() + geom_spatraster(data = dp) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 2), oob = scales::squish, na.value = NA, name = "DI / AOA threshold") +
  geom_sf(data = land, fill = "grey82", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = obs_sf, size = 0.02, colour = "#00B4D8", alpha = 0.35, inherit.aes = FALSE) +
  coord_sf(crs = PROJ, expand = FALSE, datum = NA) + theme_void(base_size = 8) +
  labs(title = "a  Geographic space: where are conditions least represented by the data?",
       subtitle = "colour = dissimilarity of each cell's conditions to the training data (DI / AOA threshold); cyan dots = observation cells (training data)") +
  theme(plot.title = element_text(face = "bold", size = 9), legend.position = "right")

# --- (c) which space guides new data (H7) --------------------------------------------------------------------------
h7 <- read.csv(file.path(DIRS$tables, "targeted_sampling_pooled_toc.csv"))
h7 <- h7[h7$strategy != "random", ]
lab <- c(novelty = "environmental novelty (DI)", uncertainty = "model uncertainty (QRF)", space = "geographic gap filling", oracle = "largest current error (greedy)")
h7$strategy_lab <- factor(lab[h7$strategy], levels = rev(lab))
h7$effort <- factor(paste0("sampled to ", h7$pct_pool, "% of candidates"), levels = paste0("sampled to ", c(40, 60), "% of candidates"))
pc3 <- ggplot(h7, aes(mean_gain_minus_random, strategy_lab, colour = effort)) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = c("#1D6A73", "#C2410C"), name = NULL) +
  labs(x = "extra RMSE reduction vs random addition of data (log10 TOC; 95% CI)", y = NULL,
       title = "c  Which space should guide new data?", subtitle = "retrospective test within 8 regions, 3 random starts each") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")

fig <- pa / (pb | pc3) + plot_layout(heights = c(1.05, 1))
ggsave(file.path(DIRS$compare, "fig_ignorance_spaces.png"), fig, width = 240, height = 230, units = "mm", dpi = 250, bg = "white")
write.csv(data.frame(axis = c("PC1", "PC2"), variance_pct = ve), file.path(DIRS$tables, "environmental_space_pca_toc.csv"), row.names = FALSE)
msg("ignorance spaces figure written")
