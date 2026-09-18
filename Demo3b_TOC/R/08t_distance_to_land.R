# =============================================================================
# 08t_distance_to_land.R — the SE Pacific wedge comes from distance to the nearest coast (08s). Fix and test.
# =============================================================================
# dist_coast_km (02_predictors.R) is the distance to the nearest land cell, including tiny oceanic islands. In the open
# ocean its field is a patchwork of island "catchments" with straight boundaries, and the random forest turns that
# geometry into TOC structure (08s: holding dist_coast_km at its median erases the wedge; holding the bottom-water
# layers does not). Candidate fix: distance to land masses >= 25,000 km2 (continents and large islands), computed the
# same way (0.2 deg, geodesic, bilinear to 0.1 deg, log10). Tests, with the design of 08q:
#   V9  iteration 2 without any distance to land      V10 iteration 2 with distance to large land masses
#   V11 iteration 1 with distance to large land masses (is the texture gain unchanged? V10 vs V11)
# Leave-one-region-out predictions as 07, budget bootstrap, and SE Pacific predictions of the full-data models.
Sys.setenv(DEMO3B_ITER = "iter2b_supply_texture")
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(ggplot2); library(patchwork); library(tidyterra) })

base  <- rast(path_base_toc_stack()); ocean <- !is.na(base[["sst_mean"]])
f_dl <- file.path(CACHE, "predictor_dist_large_land_0.1deg.tif")
if (!file.exists(f_dl)) {
  land2 <- aggregate(ifel(ocean, NA, 1), fact = 2, fun = "max", na.rm = TRUE)
  pt <- patches(land2, directions = 8)
  za <- zonal(cellSize(land2, unit = "km"), pt, fun = "sum")
  keep <- za[[1]][za[[2]] >= 25000]
  big <- ifel(pt %in% keep, 1, NA)
  dl <- resample(distance(big) / 1000, base[["sst_mean"]], method = "bilinear")
  dl <- mask(dl, ocean, maskvalues = FALSE)
  lo <- max(global(dl, function(z) quantile(z[z > 0], 0.001, na.rm = TRUE))[[1]], 1e-6)
  dl <- log10(dl + lo); names(dl) <- "dist_large_land"
  writeRaster(dl, f_dl, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  msg("distance to large land masses written: %d of %d land patches kept", length(keep), nrow(za))
}
dl <- rast(f_dl)

m2 <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds")); bt2 <- m2$bestTune
m1 <- readRDS(file.path(ROOT, "outputs", "iter1b_supply", "models", "rf_spatialCV_toc.rds")); bt1 <- m1$bestTune
V2B <- setdiff(names(m2$trainingData), c(".outcome", ".weights")); V1B <- setdiff(names(m1$trainingData), c(".outcome", ".weights"))
d <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
d$dist_large_land <- extract(dl, cbind(d$x_lon, d$y_lat))[, 1]
na <- is.na(d$dist_large_land); d$dist_large_land[na] <- d$dist_coast_km[na]
msg("training cells: rho(dist_coast_km, dist_large_land) = %.3f; %d cells where they differ by > 0.1 log10",
    cor(d$dist_coast_km, d$dist_large_land, method = "spearman"), sum(abs(d$dist_coast_km - d$dist_large_land) > 0.1))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster

swap <- function(v) replace(v, v == "dist_coast_km", "dist_large_land")
VAR <- list(V0_full = V2B, V8_iteration1_no_texture = V1B,
            V9_no_distance_to_land = setdiff(V2B, "dist_coast_km"), V10_distance_to_large_land = swap(V2B),
            V11_iteration1_distance_to_large_land = swap(V1B))
hp <- function(nm) if (nm %in% c("V8_iteration1_no_texture", "V11_iteration1_distance_to_large_land")) list(mtry = bt1$mtry, node = bt1$min.node.size) else
  list(mtry = max(2, round(length(VAR[[nm]]) * bt2$mtry / length(V2B))), node = bt2$min.node.size)
fit <- function(df, nm, trees = 300) { v <- VAR[[nm]]; h <- hp(nm)
  ranger(x = df[, v, drop = FALSE], y = df$y, num.trees = trees, mtry = min(h$mtry, length(v)), min.node.size = h$node,
         splitrule = "variance", num.threads = CFG$cores, seed = CFG$seed) }

f_lro <- file.path(DIRS$models, "distance_to_land_lro.rds")
if (!file.exists(f_lro)) {
  P0 <- readRDS(file.path(DIRS$models, "layer_sensitivity_lro.rds"))
  stopifnot(identical(P0$x_lon, d$x_lon), identical(P0$region, d$region))
  P <- P0[, c("x_lon", "y_lat", "y", "region", "V0_full", "V8_iteration1_no_texture")]
  for (nm in c("V9_no_distance_to_land", "V10_distance_to_large_land", "V11_iteration1_distance_to_large_land")) {
    p <- rep(NA_real_, nrow(d))
    for (r in sort(unique(d$region))) { tr <- d$region != r; p[!tr] <- predict(fit(d[tr, ], nm), d[!tr, VAR[[nm]], drop = FALSE], num.threads = CFG$cores)$predictions }
    P[[nm]] <- p; msg("%s: withheld RMSE %.4f", nm, sqrt(mean((p - d$y)^2)))
  }
  saveRDS(P, f_lro)
}
P <- readRDS(f_lro)

noise <- median(read.csv(file.path(DIRS$data, "toc_cells.csv"))$y_sd, na.rm = TRUE)
tex <- rast(path_texture_stack())[[paste0("litho_t", 1:6)]]
lv <- as.matrix(extract(tex, cbind(P$x_lon, P$y_lat))); ls <- rowSums(lv)
P$lith <- ifelse(is.na(ls) | ls < 0.5, "unclassified", paste0("type ", max.col(lv, ties.method = "first")))
exy <- st_coordinates(st_transform(st_as_sf(P, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
P$block <- paste(P$region, floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
budget <- function(err, region, lith) {
  reg <- ave(err, region); regl <- ave(err, region, lith); within <- mean((err - regl)^2); nf <- min(noise^2, within)
  c(total = mean(err^2), regional_offset = mean(reg^2), sediment_bias = mean((regl - reg)^2), unstructured_above_noise = within - nf)
}
idx_by_block <- split(seq_len(nrow(P)), P$block)
blocks_by_region <- split(names(idx_by_block), sub(" .*", "", names(idx_by_block)))
set.seed(CFG$seed)
IDX <- replicate(1000, unlist(lapply(blocks_by_region, function(bl) unlist(idx_by_block[sample(bl, length(bl), replace = TRUE)], use.names = FALSE)), use.names = FALSE), simplify = FALSE)
NM <- names(VAR); E <- sapply(NM, function(nm) P[[nm]] - P$y)
BS <- setNames(lapply(NM, function(nm) t(vapply(IDX, function(i) budget(E[i, nm], P$region[i], P$lith[i]), numeric(4)))), NM)
cmp <- function(a, b) { ea <- budget(E[, a], P$region, P$lith); eb <- budget(E[, b], P$region, P$lith)
  ci <- apply(BS[[b]] - BS[[a]], 2, quantile, c(0.025, 0.975))
  data.frame(reference = a, variant = b, component = names(ea), value_reference = ea, value_variant = eb, difference = eb - ea,
             ci_lo = ci[1, ], ci_hi = ci[2, ], pct_change = 100 * (eb - ea) / ea, ci_excludes_zero = ci[1, ] > 0 | ci[2, ] < 0) }
bud <- rbind(cmp("V0_full", "V9_no_distance_to_land"), cmp("V0_full", "V10_distance_to_large_land"),
             cmp("V8_iteration1_no_texture", "V11_iteration1_distance_to_large_land"),
             cmp("V11_iteration1_distance_to_large_land", "V10_distance_to_large_land"), cmp("V8_iteration1_no_texture", "V0_full"))
write.csv(bud, file.path(DIRS$compare, "distance_to_land_budget.csv"), row.names = FALSE)
print(bud[bud$component %in% c("total", "regional_offset"), ], digits = 3, row.names = FALSE)
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
reg <- do.call(rbind, lapply(sort(unique(P$region)), function(r) { j <- P$region == r; data.frame(region = RN[r], n = sum(j), t(sqrt(colMeans(E[j, , drop = FALSE]^2)))) }))
write.csv(reg, file.path(DIRS$compare, "distance_to_land_regions.csv"), row.names = FALSE)
print(reg, digits = 3, row.names = FALSE)

# SE Pacific predictions of full-data models
eB <- ext(-140, -90, -50, 0)
S <- c(crop(rast(path_toc_stack()), eB), crop(dl, eB))
X <- as.data.frame(S, na.rm = FALSE)
pz <- function(nm) { v <- VAR[[nm]]; ok <- complete.cases(X[, v]); p <- rep(NA_real_, nrow(X))
  p[ok] <- predict(fit(d, nm, trees = 500), X[ok, v, drop = FALSE], num.threads = CFG$cores)$predictions
  r <- S[[1]]; values(r) <- p
  ggplot() + geom_spatraster(data = r) + scale_fill_viridis_c(limits = c(-1.1, 0.3), oob = scales::squish, na.value = "grey85", name = "log10 TOC") +
    coord_sf(expand = FALSE) + labs(title = nm) + theme_minimal(base_size = 6.5) +
    theme(axis.text = element_blank(), plot.title = element_text(size = 6.5, face = "bold"), legend.key.width = unit(2, "mm")) }
pdl <- ggplot() + geom_spatraster(data = crop(dl, eB)) + scale_fill_viridis_c(na.value = "grey85", name = "log10 km") + coord_sf(expand = FALSE) +
  labs(title = "distance to land masses >= 25,000 km2") + theme_minimal(base_size = 6.5) + theme(axis.text = element_blank(), plot.title = element_text(size = 6.5, face = "bold"))
pdc <- ggplot() + geom_spatraster(data = crop(rast(path_toc_stack())[["dist_coast_km"]], eB)) + scale_fill_viridis_c(na.value = "grey85", name = "log10 km") +
  coord_sf(expand = FALSE) + labs(title = "distance to nearest coast (any island)") + theme_minimal(base_size = 6.5) + theme(axis.text = element_blank(), plot.title = element_text(size = 6.5, face = "bold"))
fig <- wrap_plots(list(pdc, pdl, pz("V0_full"), pz("V9_no_distance_to_land"), pz("V10_distance_to_large_land"), pz("V11_iteration1_distance_to_large_land")), ncol = 3)
ggsave(file.path(DIRS$compare, "fig_distance_to_land.png"), fig, width = 220, height = 160, units = "mm", dpi = 250, bg = "white")
msg("distance to land test done")
