# =============================================================================
# 08b_ridge_check.R — is the increase in expected error (iteration 2 - iteration 1)
# associated with mid-ocean ridges (and trenches)?
# =============================================================================
# Ridge/trench proxy from the bathymetry already in the predictor stack (no external data):
#   relief = (mean depth in a ~5 deg window) - depth; > 0 means shallower than surroundings.
#   ridge / elevated relief : 1500-4000 m and relief > 500 m
#   trench                  : > 5500 m and relief < -1000 m
# Tests: area enrichment of "red" cells (delta > 0.01) by class; delta vs relief in deep water;
# mechanism (change in DI / threshold, texture values, sampling of each class by observations).
Sys.setenv(DEMO3B_ITER = "iter2_texture")
source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(tidyterra); library(patchwork) })

f2 <- 2   # work at 0.2 deg
dlt <- aggregate(rast(file.path(DIRS$compare, "delta_expected_RMSE_iter2_minus_iter1.tif")), f2, "mean", na.rm = TRUE)
stk <- rast(path_toc_stack())
dep <- stk[["depth"]]
if (mean(values(dep), na.rm = TRUE) < 0) dep <- -dep
coarse <- aggregate(dep, fact = 5, fun = "mean", na.rm = TRUE)                   # 0.5 deg
bg <- resample(focal(coarse, w = 11, fun = "mean", na.rm = TRUE), dep, method = "bilinear")  # ~5.5 deg window
relief <- bg - dep
agg <- function(r) resample(aggregate(r, f2, "mean", na.rm = TRUE), dlt)
dep2 <- agg(dep); rel2 <- agg(relief); tri2 <- agg(stk[["tri"]])
cls <- ifel(dep2 < 1500, 1, ifel(dep2 > 5500 & rel2 < -1000, 4, ifel(dep2 <= 4000 & rel2 > 500, 3, 2)))
CL <- c("shelf & slope (<1500 m)", "abyssal plains & basins", "ridges & elevated relief (1500-4000 m, >500 m above surroundings)", "trenches (>5500 m, >1000 m below surroundings)")

thr <- function(it) { s <- read.csv(file.path(ROOT, "outputs", it, "tables", "toc_summary.csv")); as.numeric(s$value[s$metric == "AOA_DI_threshold"]) }
rr <- function(it, f) agg(rast(file.path(ROOT, "outputs", it, "rasters", f)))
dDI <- rr("iter2_texture", "toc_DI.tif") / thr("iter2_texture") - rr("iter1_baseline", "toc_DI.tif") / thr("iter1_baseline")
lsum <- agg(sum(stk[[paste0("litho_t", 1:6)]]))
v <- data.frame(delta = values(dlt)[, 1], cls = values(cls)[, 1], relief = values(rel2)[, 1], depth = values(dep2)[, 1],
                tri = values(tri2)[, 1], dDI = values(dDI)[, 1], area = values(cellSize(dlt, unit = "km"))[, 1],
                dist = values(rr("iter2_texture", "toc_distance_nearest_observation_km.tif"))[, 1],
                gs_d50 = values(agg(stk[["gs_d50"]]))[, 1], porosity = values(agg(stk[["porosity"]]))[, 1],
                litho_unclassified = values(lsum)[, 1] < 0.5)
v <- v[complete.cases(v[, c("delta", "cls", "area")]), ]
v$red <- v$delta > 0.01
obs <- read.csv(file.path(DIRS$data, "toc_cells.csv"))
ocls <- extract(cls, cbind(obs$x_lon, obs$y_lat))[, 1]

wm <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)])
tab <- do.call(rbind, lapply(1:4, function(k) { x <- v[v$cls == k, ]
  data.frame(class = CL[k], pct_ocean_area = 100 * sum(x$area) / sum(v$area),
             pct_red_area = 100 * sum(x$area[x$red]) / sum(v$area[v$red]),
             enrichment_red = (sum(x$area[x$red]) / sum(v$area[v$red])) / (sum(x$area) / sum(v$area)),
             pct_of_class_red = 100 * sum(x$area[x$red]) / sum(x$area),
             mean_delta = wm(x$delta, x$area), mean_dDI_over_threshold = wm(x$dDI, x$area),
             median_dist_obs_km = median(x$dist, na.rm = TRUE),
             pct_observation_cells = 100 * mean(ocls == k, na.rm = TRUE),
             pct_litho_unclassified = 100 * wm(x$litho_unclassified, x$area),
             mean_gs_d50 = wm(x$gs_d50, x$area), mean_porosity = wm(x$porosity, x$area)) }))
write.csv(tab, file.path(DIRS$compare, "ridge_check_by_class.csv"), row.names = FALSE)
print(tab, digits = 3)

deep <- v[v$depth > 1500, ]
set.seed(CFG$seed); s <- deep[sample(nrow(deep), min(2e5, nrow(deep))), ]
sp <- function(a, b) cor(a, b, method = "spearman", use = "complete.obs")
stats <- data.frame(test = c("rho(delta, relief), deep water", "rho(delta, TRI), deep water", "rho(delta, change in DI/threshold), deep water",
                             "rho(delta, distance to observations), deep water", "rho(delta, depth), deep water"),
                    rho = c(sp(s$delta, s$relief), sp(s$delta, s$tri), sp(s$delta, s$dDI), sp(s$delta, s$dist), sp(s$delta, s$depth)))
write.csv(stats, file.path(DIRS$compare, "ridge_check_stats.csv"), row.names = FALSE)
print(stats, digits = 3)

# delta by relief bin in deep water
deep$relief_bin <- cut(deep$relief, c(-Inf, -1500, -1000, -500, -250, 0, 250, 500, 750, 1000, 1500, Inf))
byb <- do.call(rbind, lapply(split(deep, deep$relief_bin), function(x) if (nrow(x))
  data.frame(relief_bin = x$relief_bin[1], mean_delta = wm(x$delta, x$area), pct_red = 100 * wm(x$red, x$area), n = nrow(x))))
write.csv(byb, file.path(DIRS$compare, "ridge_check_by_relief_bin.csv"), row.names = FALSE)

# --- figure: Atlantic zoom + relief response --------------------------------------------------------------
box <- ext(-70, 20, -60, 65)
d_c <- crop(dlt, box); c_c <- crop(cls, box)
pa <- ggplot() + geom_spatraster(data = d_c) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-0.03, 0.03), oob = scales::squish, na.value = "grey80", name = "change in\nexpected RMSE") +
  coord_sf(expand = FALSE) + labs(title = "a  Change in expected error (iteration 2 - 1)") + theme_minimal(base_size = 8)
pb <- ggplot() + geom_spatraster(data = as.factor(c_c)) +
  scale_fill_manual(values = c("1" = "#D9D9D9", "2" = "#F2F2F2", "3" = "#D95F02", "4" = "#1B9E77"), labels = c("shelf & slope", "abyssal", "ridges & elevated relief", "trenches"),
                    na.value = "grey60", name = NULL, na.translate = FALSE) +
  coord_sf(expand = FALSE) + labs(title = "b  Bathymetric classes (relief proxy)") + theme_minimal(base_size = 8)
pc <- ggplot(byb, aes(relief_bin, mean_delta)) + geom_hline(yintercept = 0, colour = "grey50") + geom_col(fill = "#B2182B") +
  labs(x = "relief: m shallower than surroundings (deep water > 1500 m)", y = "mean change in expected RMSE",
       title = "c  Change in expected error vs relief") + theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig <- (pa | pb) / pc + plot_layout(heights = c(1.6, 1))
ggsave(file.path(DIRS$compare, "fig_ridge_check.png"), fig, width = 200, height = 220, units = "mm", dpi = 250, bg = "white")
msg("ridge check done")
