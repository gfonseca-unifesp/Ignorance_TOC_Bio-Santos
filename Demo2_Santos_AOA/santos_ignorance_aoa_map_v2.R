## ================================================================
## Demo 2b (v2) — Santos Basin: ignorance field combined with the AOA
## ================================================================
## Run after santos_aoa.R (v1 figure and script are kept unchanged):
##   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_ignorance_aoa_map_v2.R
##
## Changes from v1
##  1. Background: biological-tier DI / threshold shown unsmoothed. The only
##     operation is filling isolated empty raster cells.
##  2. Soft AOA edge for the environmental tier, instead of a hard line at 1:
##     inside (< 0.8), transition (0.8-1.2), extrapolation (> 1.2).
##     The Q75 + 1.5 IQR threshold is a convention; the band shows that.
##  3. The class map is generalised with a sieve (patches < SIEVE_CELLS cells merged into
##     their neighbours). A sensitivity table and figure (0-320 cells) show the area of each
##     class does not hinge on that choice. Statistics use the unsieved classes.
##  4. Diagnostic of the rectangular blocks in the background: do their edges coincide
##     with the most frequent Long/Lat split points of the environmental random forests?
## ================================================================

.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(caret); library(randomForest); library(terra); library(sf)
  library(ggplot2); library(patchwork); library(rnaturalearth)
})
set.seed(20260914)

ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")  # F0.2: ver .Renviron.example
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs")
RES  <- 0.02                          # raster resolution (deg); grid spacing is ~2 km
SIEVE_CELLS <- 40                     # generalisation used in the map
SIEVE_GRID  <- c(0, 10, 20, 40, 80, 160, 320)
EDGE <- c(0.8, 1.2)                   # soft-edge band of DI / threshold
BIO <- c("Meio_N", "Meio_Nema", "Meio_Cop", "Meio_Kino", "Meio_Poly", "Meio_S",
         "Macro_N", "Macro_Annel", "Macro_Ploy", "Macro_Crust", "Macro_Moll", "Macro_S")
SURVEY <- c(`1` = "2019 (Survey 1)", `2` = "2021 (Survey 2)")
CLASS_LAB <- c("inside AOA", sprintf("transition (DI/threshold %.1f-%.1f)", EDGE[1], EDGE[2]),
               sprintf("extrapolation (DI/threshold > %.1f)", EDGE[2]))

## ---- 1. station ignorance (S2) ---------------------------------------------------
sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
W <- sd1$Wat_Sed_Org_Geo_mean2; G <- sd1$Geo_mean; B <- sd1[["Biod_mean_Numeric_2024-05-20_log10"]]
rfW <- attr(W, "rf"); rfG <- attr(G, "rf"); rm(sd1); invisible(gc())
stopifnot(identical(rownames(B), rownames(W)), identical(rownames(G), rownames(W)))

O <- as.matrix(B[, BIO]); P <- O * NA
for (v in BIO) {
  m  <- rfW[[paste0(v, "~Wat_Sed_Org_Geo_mean2")]]$m
  pr <- m$pred
  for (tn in names(m$bestTune)) pr <- pr[pr[[tn]] == m$bestTune[[tn]], ]
  oof <- tapply(pr$pred, pr$rowIndex, mean)
  P[rownames(m$trainingData)[as.integer(names(oof))], v] <- oof
  te <- setdiff(rownames(W), rownames(m$trainingData))
  P[te, v] <- predict(m, W[te, predictors(m), drop = FALSE])
}
stopifnot(!anyNA(P))
st <- read.csv(file.path(OUT, "santos_aoa_stations.csv"))
stopifnot(identical(st$station, rownames(W)))
st$ignorance <- rowMeans(sweep(abs(O - P), 2, apply(O, 2, sd), "/"))
st$survey <- SURVEY[as.character(st$Camp)]

## ---- 2. grid layers -> rasters ---------------------------------------------------
x <- read.csv(file.path(OUT, "santos_aoa_consensus.csv"))
tmpl <- rast(ext(min(x$Long) - RES, max(x$Long) + RES, min(x$Lat) - RES, max(x$Lat) + RES),
             resolution = RES, crs = "EPSG:4326")
to_rast <- function(d, field) {
  r <- rasterize(as.matrix(d[, c("Long", "Lat")]), tmpl, values = d[[field]], fun = mean)
  focal(r, 3, "mean", na.policy = "only", na.rm = TRUE)          # fills isolated empty cells only
}
classify_edge <- function(r) classify(r, rbind(c(-Inf, EDGE[1], 0), c(EDGE[1], EDGE[2], 1), c(EDGE[2], Inf, 2)),
                                      include.lowest = TRUE, right = FALSE)
sieve_class <- function(cls, n) {
  if (n == 0) return(cls)
  filled <- ifel(is.na(cls), 0, cls)
  mask(sieve(filled, threshold = n, directions = 8), cls)
}
layers <- list(); sens <- list()
for (cp in c("1", "2")) {
  d   <- x[x$Camp == as.integer(cp), ]
  bio <- to_rast(d, "DInorm_mean"); env <- to_rast(d, "env_DInorm_mean")
  cls <- classify_edge(env)
  cls_s <- sieve_class(cls, SIEVE_CELLS)
  area <- cellSize(cls, unit = "km")
  for (n in SIEVE_GRID) {
    cn <- sieve_class(cls, n)
    tot <- global(ifel(is.na(cn), NA, area), "sum", na.rm = TRUE)[[1]]
    sens[[length(sens) + 1]] <- data.frame(survey = SURVEY[cp], sieve_cells = n,
      pct_transition    = 100 * global(ifel(cn == 1, area, NA), "sum", na.rm = TRUE)[[1]] / tot,
      pct_extrapolation = 100 * global(ifel(cn == 2, area, NA), "sum", na.rm = TRUE)[[1]] / tot)
  }
  s <- c(bio, env, cls, cls_s)
  names(s) <- c("bio_DInorm", "env_DInorm", "env_edge_class", "env_edge_class_sieved")
  writeRaster(s, file.path(OUT, sprintf("santos_ignorance_aoa_layers_v2_camp%s.tif", cp)), overwrite = TRUE,
              gdal = "COMPRESS=DEFLATE")
  layers[[cp]] <- s
}
sens <- do.call(rbind, sens)
write.csv(sens, file.path(OUT, "santos_aoa_sieve_sensitivity_v2.csv"), row.names = FALSE)

## ---- 3. block diagnostic: RF split points on Long / Lat ------------------------------
env_models <- setNames(lapply(rfG, `[[`, "m"), sub("~Geo_mean$", "", names(rfG)))
imp_share <- t(sapply(env_models, function(m) {
  v <- varImp(m, scale = FALSE)$importance; s <- pmax(setNames(v[, 1], rownames(v)), 0)
  s[c("Long", "Lat", "Camp.1", "Depth")] / sum(s) }))
imp_share[is.na(imp_share)] <- 0
write.csv(data.frame(env_variable = rownames(imp_share), round(imp_share, 3)),
          file.path(OUT, "santos_env_models_importance_share_v2.csv"), row.names = FALSE)
top_var <- rownames(imp_share)[which.max(imp_share[, "Long"] + imp_share[, "Lat"])]
fm <- env_models[[top_var]]$finalModel
splits <- do.call(rbind, lapply(seq_len(fm$ntree), function(k) {
  t <- randomForest::getTree(fm, k, labelVar = TRUE)
  t <- t[!is.na(t$`split var`) & t$`split var` %in% c("Long", "Lat"), c("split var", "split point")]
  if (nrow(t)) t else NULL }))
names(splits) <- c("var", "point")
splits$var <- as.character(splits$var)                       # drop empty factor levels (Camp.1, Depth)
top_splits <- do.call(rbind, lapply(split(splits, splits$var), function(s) {
  tb <- head(sort(table(round(s$point, 3)), decreasing = TRUE), 8)
  data.frame(var = s$var[1], point = as.numeric(names(tb)), n_trees_nodes = as.integer(tb)) }))
write.csv(top_splits, file.path(OUT, "santos_rf_top_splits_v2.csv"), row.names = FALSE)

d1 <- x[x$Camp == 1, ]
EG1 <- readRDS(file.path(ROOT, "Savepoint_part2.rds"))$saved_data$Env_predictions_grid[x$Camp == 1, top_var]
d1$val <- EG1
rv <- to_rast(d1, "val"); rv_df <- as.data.frame(rv, xy = TRUE, na.rm = TRUE); names(rv_df)[3] <- "val"
# straight step edges: column/row boundaries carrying >= MIN_RUN strong value jumps (top 5% of jumps)
MIN_RUN <- 15                                              # ~30 km of edge
mv <- as.matrix(rv, wide = TRUE)
dx <- abs(mv[, -1, drop = FALSE] - mv[, -ncol(mv), drop = FALSE])   # between column j and j+1
dy <- abs(mv[-1, , drop = FALSE] - mv[-nrow(mv), , drop = FALSE])   # between row i and i+1
step_thr <- as.numeric(quantile(c(dx, dy), 0.95, na.rm = TRUE))
xs <- xFromCol(rv, which(colSums(dx > step_thr, na.rm = TRUE) >= MIN_RUN)) + RES / 2
ys <- yFromRow(rv, which(rowSums(dy > step_thr, na.rm = TRUE) >= MIN_RUN)) - RES / 2
# Permutation test: split-use density at the step edges vs at random positions.
# score(z) = number of forest nodes splitting on that coordinate within +/- 1.5 cells of z.
# Observed mean score over the detected step lines is compared with 2000 sets of the same
# number of positions drawn uniformly over the raster's extent on that axis.
N_PERM <- 2000
score <- function(z, pts) vapply(z, function(q) sum(abs(pts - q) <= 1.5 * RES), numeric(1))
perm_test <- function(edges, pts, lo, hi) {
  if (!length(edges)) return(c(ratio = NA, p = NA))
  obs  <- mean(score(edges, pts))
  null <- replicate(N_PERM, mean(score(runif(length(edges), lo, hi), pts)))
  c(ratio = obs / mean(null), p = (1 + sum(null >= obs)) / (N_PERM + 1))
}
set.seed(20260914)
tL <- perm_test(xs, splits$point[splits$var == "Long"], xmin(rv), xmax(rv))
tA <- perm_test(ys, splits$point[splits$var == "Lat"],  ymin(rv), ymax(rv))
blk <- data.frame(variable = top_var,
                  long_lat_importance_share = round(sum(imp_share[top_var, c("Long", "Lat")]), 3),
                  median_long_lat_share_44_models = round(median(imp_share[, "Long"] + imp_share[, "Lat"]), 3),
                  n_vertical_step_lines = length(xs), long_split_density_ratio = round(tL[["ratio"]], 2), long_perm_p = signif(tL[["p"]], 3),
                  n_horizontal_step_lines = length(ys), lat_split_density_ratio = round(tA[["ratio"]], 2), lat_perm_p = signif(tA[["p"]], 3))
write.csv(blk, file.path(OUT, "santos_rf_block_check_v2.csv"), row.names = FALSE)
print(blk)

land <- ne_countries(scale = 50, returnclass = "sf")
map_base <- list(geom_sf(data = land, fill = "#E6DFCF", colour = "grey45", linewidth = 0.2, inherit.aes = FALSE),
                 coord_sf(xlim = c(-49, -40.6), ylim = c(-28.1, -22.6), expand = FALSE),
                 theme_bw(base_size = 8), theme(panel.grid = element_blank(), legend.position = "bottom",
                                                legend.key.width = unit(10, "mm")))
gb <- ggplot() + geom_raster(data = rv_df, aes(x, y, fill = val)) + scale_fill_viridis_c(name = top_var) +
  geom_vline(xintercept = top_splits$point[top_splits$var == "Long"], colour = "red", linewidth = 0.25, linetype = 2) +
  geom_hline(yintercept = top_splits$point[top_splits$var == "Lat"], colour = "red", linewidth = 0.25, linetype = 2) +
  map_base + labs(x = "Longitude", y = "Latitude",
    title = sprintf("Block check: predicted %s (2019) and the 8 most frequent Long/Lat split points of its random forest (red)", top_var),
    subtitle = sprintf(paste0("Long+Lat importance share %.0f%% (median over the 44 environmental models %.0f%%).\n",
                              "Split-use density at straight value steps vs random positions: vertical x%.2f (p = %.3g), horizontal x%.2f (p = %.3g)"),
                       100 * blk$long_lat_importance_share, 100 * blk$median_long_lat_share_44_models,
                       blk$long_split_density_ratio, blk$long_perm_p, blk$lat_split_density_ratio, blk$lat_perm_p))
ggsave(file.path(OUT, "santos_rf_block_check_v2.png"), gb, width = 210, height = 150, units = "mm", dpi = 250, bg = "white")

## ---- 4. statistics --------------------------------------------------------------------
sp <- function(a, b) { ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)); c(unname(ct$estimate), ct$p.value) }
r_all <- sp(st$ignorance, st$bio_DInorm_mean); r_te <- sp(st$ignorance[st$set == "test"], st$bio_DInorm_mean[st$set == "test"])
row_cls <- cut(x$env_DInorm_mean, c(-Inf, EDGE, Inf), labels = CLASS_LAB, right = FALSE)
S <- function(metric, value, description) data.frame(metric, value = round(value, 3), description)
summ <- rbind(
  S("ignorance_mean_2019", mean(st$ignorance[st$Camp == 1]), "reproduces S2 (0.39)"),
  S("ignorance_mean_2021", mean(st$ignorance[st$Camp == 2]), "reproduces S2 (0.46)"),
  S("spearman_ignorance_DI_all", r_all[1], "198 stations"), S("spearman_ignorance_DI_all_p", r_all[2], ""),
  S("spearman_ignorance_DI_test", r_te[1], "38 held-out stations"), S("spearman_ignorance_DI_test_p", r_te[2], ""),
  S("pct_grid_env_inside", 100 * mean(row_cls == CLASS_LAB[1]), "grid rows, unsieved"),
  S("pct_grid_env_transition", 100 * mean(row_cls == CLASS_LAB[2]), "grid rows, unsieved"),
  S("pct_grid_env_extrapolation", 100 * mean(row_cls == CLASS_LAB[3]), "grid rows, unsieved"),
  S("pct_deep_rows_extrapolation", 100 * mean(row_cls[x$Depth > 2400] == CLASS_LAB[3]), "cells deeper than 2400 m"),
  S("sieve_cells_in_map", SIEVE_CELLS, "patches smaller than this merged into neighbours (display only)"),
  S("sieve_range_pct_extrapolation_max_minus_min",
    diff(range(sens$pct_extrapolation[sens$sieve_cells > 0])), "percentage points across sieve 10-320 cells")
)
write.csv(summ, file.path(OUT, "santos_ignorance_aoa_summary_v2.csv"), row.names = FALSE)
print(summ, right = FALSE); print(sens, digits = 3)

## ---- 5. combined map (v2) ----------------------------------------------------------------
hatch_for <- function(poly, spacing = 0.06) {
  bb <- st_bbox(poly); w <- bb[["xmax"]] - bb[["xmin"]]; h <- bb[["ymax"]] - bb[["ymin"]]
  offs <- seq(-h, w, by = spacing)
  lines <- st_sfc(lapply(offs, function(o) st_linestring(rbind(c(bb[["xmin"]] + o, bb[["ymin"]]),
                                                             c(bb[["xmin"]] + o + h, bb[["ymax"]])))), crs = 4326)
  suppressWarnings(st_intersection(lines, st_union(poly)))
}
bio_df <- do.call(rbind, lapply(names(layers), function(cp) {
  d <- as.data.frame(layers[[cp]][["bio_DInorm"]], xy = TRUE, na.rm = TRUE); d$survey <- SURVEY[cp]; d }))
polys <- list(); hatches <- list()
for (cp in names(layers)) {
  cl <- layers[[cp]][["env_edge_class_sieved"]]
  for (k in 1:2) {
    p <- st_as_sf(as.polygons(ifel(cl == k, 1, NA), dissolve = TRUE))
    if (!nrow(p)) next
    p <- st_sf(geometry = st_union(p)); p$survey <- SURVEY[cp]; p$class <- CLASS_LAB[k + 1]
    polys[[length(polys) + 1]] <- p
    if (k == 2) { h <- st_sf(geometry = hatch_for(p)); h$survey <- SURVEY[cp]; hatches[[length(hatches) + 1]] <- h }
  }
}
polys <- do.call(rbind, polys); hatches <- do.call(rbind, hatches)
polys$class <- factor(polys$class, levels = CLASS_LAB[2:3])
ign_cap <- as.numeric(quantile(st$ignorance, 0.975))
lab <- data.frame(survey = SURVEY, txt = sprintf("mean ignorance = %.2f (n = 99)", tapply(st$ignorance, st$Camp, mean)))

g <- ggplot() +
  geom_sf(data = land, fill = "#E6DFCF", colour = "grey45", linewidth = 0.2) +
  geom_raster(data = bio_df, aes(x, y, fill = bio_DInorm)) +
  scale_fill_gradient(low = "#F4F6F8", high = "#3E5C76", limits = c(0, 1),
                      name = "Dissimilarity to sampled conditions\n(biological tier, DI / AOA threshold; unsmoothed)") +
  geom_sf(data = polys[polys$class == CLASS_LAB[2], ], fill = "grey20", alpha = 0.10, colour = NA) +
  geom_sf(data = polys[polys$class == CLASS_LAB[3], ], fill = "grey20", alpha = 0.20, colour = NA) +
  geom_sf(data = hatches, colour = "grey25", linewidth = 0.15) +
  geom_sf(data = polys, aes(linetype = class), fill = NA, colour = "grey20", linewidth = 0.3) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Environmental predictions",
                        guide = guide_legend(override.aes = list(fill = c("grey85", "grey65")), ncol = 1)) +
  geom_point(data = st, aes(Long, Lat, colour = pmin(ignorance, ign_cap)), size = 2.3) +
  geom_point(data = st, aes(Long, Lat), shape = 21, size = 2.3, colour = "grey15", stroke = 0.25) +
  scale_colour_distiller(palette = "YlOrRd", direction = 1,
                         name = "Station ignorance\n(mean standardised |obs - pred|, out-of-sample)") +
  geom_text(data = lab, aes(x = -48.85, y = -22.72, label = txt), hjust = 0, vjust = 1, size = 2.7) +
  facet_wrap(~survey) +
  coord_sf(xlim = c(-49, -40.6), ylim = c(-28.1, -22.6), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude",
       title = "Santos Basin: where the benthic models are ignorant, and where they extrapolate",
       # F1.2: era "(Suppl. Fig. S2)" — numeração antiga do suplemento
       subtitle = sprintf(paste0("Points: out-of-sample ignorance (mean over 12 indicators). Background: dissimilarity of the modelled environment to the training stations.\n",
                                 "Environmental-tier AOA edge as a band: dotted = transition (DI/threshold %.1f-%.1f), hatched = extrapolation (> %.1f); patches < %d cells merged for display.\n",
                                 "Stations: ignorance vs dissimilarity, Spearman rho = %.2f (held-out test stations %.2f)."),
                          EDGE[1], EDGE[2], EDGE[2], SIEVE_CELLS, r_all[1], r_te[1])) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.title.position = "top",
        legend.key.width = unit(9, "mm"), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"), strip.background = element_rect(fill = "grey92"))
ggsave(file.path(OUT, "santos_ignorance_aoa_map_v2.png"), g, width = 270, height = 165, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "santos_ignorance_aoa_map_v2.pdf"), g, width = 270, height = 165, units = "mm", bg = "white")

gs <- ggplot(reshape(sens, direction = "long", varying = c("pct_transition", "pct_extrapolation"), v.names = "pct",
                     timevar = "class", times = CLASS_LAB[2:3]),
             aes(sieve_cells, pct, colour = class, linetype = survey)) +
  geom_vline(xintercept = SIEVE_CELLS, colour = "grey60") + geom_line() + geom_point(size = 1) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 10), breaks = SIEVE_GRID) +
  scale_colour_manual(values = c("grey55", "grey15"), name = NULL) +
  labs(x = "Sieve: minimum patch size (raster cells; 1 cell ~ 4.9 km2)", y = "% of grid area",
       title = "Sensitivity of the environmental-tier AOA classes to map generalisation") +
  theme_classic(base_size = 8) + theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(OUT, "santos_aoa_sieve_sensitivity_v2.png"), gs, width = 140, height = 100, units = "mm", dpi = 250, bg = "white")
writeLines(capture.output(sessioninfo::session_info()), file.path(OUT, "sessionInfo_ignorance_map_v2.txt"))
cat("done\n")
