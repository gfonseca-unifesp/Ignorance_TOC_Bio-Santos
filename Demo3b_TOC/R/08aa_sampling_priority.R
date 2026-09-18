# =============================================================================
# 08aa_sampling_priority.R — F2 of ROADMAP_MSv6: where should the next TOC samples be taken?
# =============================================================================
# The ignorance map says where the model is ignorant. This turns it into an instruction: N ranked
# places where one new TOC sample would reduce ignorance most, read together with the per-region
# decision rule (08g/08l). Run on the final map:
#   DEMO3B_ITER=iter2c_clean Rscript R/08aa_sampling_priority.R
#
# What is ranked, and what that means
#   * Main ranking — REPRESENTATIVE NOVELTY. A site is worth more when it represents a large area of
#     unsampled conditions, not when it is a rare outlier. Gain of a candidate s = sum over its K_NN
#     environmental neighbours c of max(0, DI_c - d(c,s)/avg): how much the dissimilarity of the
#     neighbourhood drops if s is sampled. Greedy: pick the best, update the DI of every evaluation
#     cell, repeat. The AOA threshold is held fixed (it would change if the training set changed);
#     this is an approximation, declared as such.
#   * Sensitivity ranking — MODEL UNCERTAINTY (width of the 80% QRF interval), the strategy that won
#     the retrospective test in most regions (08h, F5.3). Reported with its overlap (Jaccard at 300 km)
#     with the main ranking, so the reader sees how much the guide depends on the criterion.
#   * The gain curve measures ENVIRONMENTAL ignorance (DI) and is, by construction, favourable to
#     novelty. Evidence that directed data reduce ERROR comes from the retrospective test (panel d).
# Nothing here is overwritten: every output is new and goes to outputs/comparison/.
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter2c_clean"))
source("R/00_config.R")
suppressPackageStartupMessages({
  library(ranger); library(FNN); library(ggplot2); library(patchwork); library(tidyterra); library(rnaturalearth); library(ggnewscale)
})
set.seed(CFG$seed)

N_SITES <- as.integer(Sys.getenv("DEMO3B_N_SITES", "50"))     # ranked sites to return
N_EVAL  <- as.integer(Sys.getenv("DEMO3B_N_EVAL",  "60000"))  # area-weighted ocean sample used to evaluate gain
N_CAND  <- as.integer(Sys.getenv("DEMO3B_N_CAND",  "20000"))  # candidate sites, drawn from the evaluation sample
K_NN    <- as.integer(Sys.getenv("DEMO3B_K_NN",    "400"))    # environmental neighbours that a candidate can represent
SEP_KM  <- 300                                                # minimum separation of the uncertainty ranking (08h block size)
N_RAND  <- 20                                                 # random site sets for the gain curve
t_start <- Sys.time()

## ---- 1. inputs ---------------------------------------------------------------------------------
model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
tdi   <- readRDS(file.path(DIRS$models, "trainDI_toc.rds"))
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
summ  <- read.csv(file.path(DIRS$tables, "toc_summary.csv"))
thr   <- as.numeric(summ$value[summ$metric == "AOA_DI_threshold"])
coll  <- read.csv(file.path(DIRS$tables, "collection_by_region_toc.csv"))   # F5.3: iter2c, 10 random starts
stk   <- rast(path_toc_stack())
r_DI   <- rast(file.path(DIRS$rasters, "toc_DI.tif"))
r_AOA  <- rast(file.path(DIRS$rasters, "toc_AOA.tif"))
r_pred <- rast(file.path(DIRS$rasters, "toc_prediction_log10.tif"))
r_err  <- rast(file.path(DIRS$rasters, "toc_ignorance_calibrated_expectedRMSE.tif"))
r_dist <- rast(file.path(DIRS$rasters, "toc_distance_nearest_observation_km.tif"))
vars <- tdi$variables
msg("iteration %s | %d predictors | AOA threshold %.4f | trainDist_avrgmean %.4f", ITER, length(vars), thr, tdi$trainDist_avrgmean)

# the DI space of CAST: centre and scale with the training parameters, weight by variable importance,
# and express distances as multiples of the mean nearest-neighbour distance within the training data
w   <- unlist(tdi$weight)[vars]
ctr <- tdi$scaleparam[["scaled:center"]][vars]
scl <- tdi$scaleparam[["scaled:scale"]][vars]
wz  <- function(M) sweep(sweep(sweep(as.matrix(M)[, vars, drop = FALSE], 2, ctr, "-"), 2, scl, "/"), 2, w, "*")
Zt  <- wz(d[, vars])
di_of <- function(Z) get.knnx(Zt, Z, k = 1)$nn.dist[, 1] / tdi$trainDist_avrgmean
# dissimilarity of every evaluation cell to a set of candidate sites, in the same units as the DI
di_to_sites <- function(Ze, Zs) get.knnx(Zs, Ze, k = 1)$nn.dist[, 1] / tdi$trainDist_avrgmean

## ---- 2. mandatory check: does the reproduction match the published DI raster? --------------------
oc <- which(!is.na(values(r_DI)[, 1]))
set.seed(CFG$seed); chk <- sample(oc, 2000)
Xc <- stk[[vars]][chk]; okc <- complete.cases(Xc); chk <- chk[okc]
di_ref <- values(r_DI)[chk, 1]; di_rep <- di_of(wz(Xc[okc, ]))
chk_cor <- cor(di_ref, di_rep); chk_max <- max(abs(di_ref - di_rep))
msg("DI check on %d cells: correlation %.8f, maximum absolute error %.3g", length(chk), chk_cor, chk_max)
if (!(chk_cor > 0.999)) stop("DI reproduction failed (correlation <= 0.999): the sampling priorities would not be those of the published map")

## ---- 3. evaluation sample and candidates --------------------------------------------------------
area <- values(cellSize(r_DI, unit = "km"))[oc, 1]
set.seed(CFG$seed)
ev <- oc[order(rexp(length(oc)) / area)[1:min(N_EVAL, length(oc))]]     # area-weighted, without replacement
Xe <- stk[[vars]][ev]; oke <- complete.cases(Xe); ev <- ev[oke]; Xe <- Xe[oke, ]
Ze <- wz(Xe)
DI0 <- values(r_DI)[ev, 1]                     # published DI of the evaluation cells
xy  <- xyFromCell(r_DI, ev)
msg("evaluation sample: %d ocean cells (area-weighted); mean DI/threshold %.3f", length(ev), mean(DI0) / thr)

set.seed(CFG$seed)
hi <- order(-DI0)[1:(N_CAND %/% 2)]                                     # half from the least represented conditions
ca <- sample(setdiff(seq_along(ev), hi), N_CAND - length(hi))           # half at random over the ocean
cand <- sort(c(hi, ca))
msg("candidates: %d (half highest DI, half at random)", length(cand))

## ---- 4. main ranking: representative novelty (greedy) -------------------------------------------
nn <- get.knnx(Ze, Ze[cand, , drop = FALSE], k = K_NN)                  # fixed geometry: neighbours of each candidate
nn_idx <- nn$nn.index; nn_dn <- nn$nn.dist / tdi$trainDist_avrgmean     # distances in DI units
rm(nn); invisible(gc())

DI <- DI0                                                              # DI of every evaluation cell, updated as sites are added
picked <- integer(0); gain_at <- numeric(0); di_at <- numeric(0)
for (k in seq_len(N_SITES)) {
  g <- rowSums(pmax(matrix(DI[nn_idx], nrow = nrow(nn_idx)) - nn_dn, 0))   # pmax(0, M) would drop the dimensions
  g[cand %in% picked] <- -Inf
  b <- which.max(g)
  picked <- c(picked, cand[b]); gain_at <- c(gain_at, g[b]); di_at <- c(di_at, DI[cand[b]])
  DI <- pmin(DI, di_to_sites(Ze, Ze[cand[b], , drop = FALSE]))
  if (k %% 10 == 0) msg("  novelty: %d sites; mean DI/threshold %.4f (start %.4f)", k, mean(DI) / thr, mean(DI0) / thr)
}
DI_novelty <- DI

## ---- 5. sensitivity ranking: model uncertainty (QRF, 300-km separation) --------------------------
bt <- model$bestTune
qrf <- ranger(x = d[, vars, drop = FALSE], y = d$y, num.trees = 300, mtry = min(bt$mtry, length(vars)),
              min.node.size = bt$min.node.size, splitrule = "variance", quantreg = TRUE,
              num.threads = CFG$cores, seed = CFG$seed)
qq <- predict(qrf, as.data.frame(Xe), type = "quantiles", quantiles = c(0.1, 0.9), num.threads = CFG$cores)$predictions
width80 <- qq[, 2] - qq[, 1]
to_xyz   <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))
X3e <- to_xyz(xy[, 1], xy[, 2])
ord <- cand[order(-width80[cand])]
pick_u <- integer(0)
for (i in ord) {
  if (length(pick_u) == N_SITES) break
  if (length(pick_u) && min(chord_km(get.knnx(X3e[pick_u, , drop = FALSE], X3e[i, , drop = FALSE], k = 1)$nn.dist[, 1])) < SEP_KM) next
  pick_u <- c(pick_u, i)
}
msg("uncertainty ranking: %d sites, minimum separation %d km", length(pick_u), SEP_KM)

# overlap of the two rankings: a site of one ranking counts as matched if the other has a site within 300 km
matched <- sum(chord_km(get.knnx(X3e[pick_u, , drop = FALSE], X3e[picked, , drop = FALSE], k = 1)$nn.dist[, 1]) < SEP_KM)
jacc <- matched / (length(picked) + length(pick_u) - matched)
msg("overlap novelty vs uncertainty: %d of %d sites matched within %d km (Jaccard %.2f)", matched, length(picked), SEP_KM, jacc)

## ---- 6. gain curve: mean DI / threshold of the ocean against the number of sites -----------------
curve_of <- function(site_idx) {
  DIc <- DI0; out <- numeric(length(site_idx) + 1); out[1] <- mean(DIc) / thr
  for (k in seq_along(site_idx)) {
    DIc <- pmin(DIc, di_to_sites(Ze, Ze[site_idx[k], , drop = FALSE]))
    out[k + 1] <- mean(DIc) / thr
  }
  out
}
cv_nov <- curve_of(picked)
cv_unc <- curve_of(pick_u)
set.seed(CFG$seed)
cv_rnd <- vapply(seq_len(N_RAND), function(i) curve_of(sample(cand, N_SITES)), numeric(N_SITES + 1))
curve <- data.frame(n_sites = 0:N_SITES,
                    novelty = cv_nov, uncertainty = cv_unc,
                    random_mean = rowMeans(cv_rnd), random_min = apply(cv_rnd, 1, min), random_max = apply(cv_rnd, 1, max))
write.csv(curve, file.path(DIRS$compare, "sampling_priority_curve_toc.csv"), row.names = FALSE)
red <- function(v, k) 100 * (v[1] - v[k + 1]) / v[1]
msg("mean DI/threshold %.4f -> %.4f with 10 sites and %.4f with 50 (novelty); random 50 sites -> %.4f",
    cv_nov[1], cv_nov[11], cv_nov[N_SITES + 1], curve$random_mean[N_SITES + 1])

## ---- 7. annotate the sites ----------------------------------------------------------------------
set.seed(CFG$seed)                                    # the 8 regions of 05 / 07 / 08g / 08h
km <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic",
        "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
region_of <- function(lon, lat) RN[apply(as.matrix(dist(rbind(to_xyz(lon, lat), km$centers)))[seq_along(lon), length(lon) + seq_len(nrow(km$centers)), drop = FALSE], 1, which.min)]

annotate <- function(idx, gain = NA, di_sel = NA) {
  cells <- ev[idx]; p <- xy[idx, , drop = FALSE]
  reg <- region_of(p[, 1], p[, 2])
  out <- data.frame(
    rank = seq_along(idx), lon = round(p[, 1], 3), lat = round(p[, 2], 3),
    depth_m = round(values(stk[["depth"]])[cells, 1]),
    toc_predicted_pct = round(10^values(r_pred)[cells, 1] - CFG$toc_offset, 3),
    DInorm_initial = round(DI0[idx] / thr, 3),
    DInorm_at_selection = round(di_sel / thr, 3),
    gain = signif(gain, 4),
    expected_RMSE_calibrated = round(values(r_err)[cells, 1], 3),
    distance_nearest_observation_km = round(values(r_dist)[cells, 1]),
    inside_AOA = values(r_AOA)[cells, 1] == 1,
    region = reg)
  out <- merge(out, coll[, c("region", "decision", "path", "collection")], by = "region", all.x = TRUE, sort = FALSE)
  out <- out[order(out$rank), ]
  # flags, by declared rule — flagged, never silently dropped
  out$flag_region_needs_new_drivers <- grepl("^3 ", out$path)
  out$flag_outside_region_support   <- out$distance_nearest_observation_km > 1000
  rownames(out) <- NULL
  out[, c("rank", "lon", "lat", "depth_m", "toc_predicted_pct", "DInorm_initial", "DInorm_at_selection", "gain",
          "expected_RMSE_calibrated", "distance_nearest_observation_km", "inside_AOA", "region", "path",
          "decision", "collection", "flag_region_needs_new_drivers", "flag_outside_region_support")]
}
sites <- annotate(picked, gain_at, di_at)
sites_u <- annotate(pick_u, NA, DI0[pick_u])
write.csv(sites, file.path(DIRS$compare, "sampling_priority_sites_toc.csv"), row.names = FALSE)
write.csv(sites_u, file.path(DIRS$compare, "sampling_priority_uncertainty_sites_toc.csv"), row.names = FALSE)

by_reg <- as.data.frame(table(region = factor(sites$region, levels = RN)), responseName = "n_sites_novelty")
by_reg$n_sites_uncertainty <- as.integer(table(factor(sites_u$region, levels = RN)))
by_reg <- merge(by_reg, coll[, c("region", "path", "collection", "T_minus_G", "advantage_vs_random")], by = "region", all.x = TRUE)
by_reg$mean_DInorm_of_sites <- round(tapply(sites$DInorm_initial, factor(sites$region, levels = RN), mean)[by_reg$region], 3)
# regions are k-means clusters of the OBSERVATION cells, so an unsampled ocean (the South Pacific gyre) is
# attached to the nearest centroid, which may be an ocean away. That is what the flag below counts.
by_reg$n_sites_outside_region_support <- as.integer(tapply(sites$flag_outside_region_support,
                                                           factor(sites$region, levels = RN), sum)[by_reg$region])
by_reg$median_distance_nearest_observation_km <- round(tapply(sites$distance_nearest_observation_km,
                                                              factor(sites$region, levels = RN), median)[by_reg$region])
write.csv(by_reg, file.path(DIRS$compare, "sampling_priority_by_region_toc.csv"), row.names = FALSE)
print(by_reg, row.names = FALSE)

## ---- 7b. are the first sites layer artefacts? (F2 acceptance; numbers quoted in the text) ----------
# Compared with random ocean cells: depth, distance to land, dissimilarity, and the share of the 8 neighbouring
# cells whose sample-based lithology differs materially (> 0.2 in any lithology layer) - i.e. lithology edges,
# where predictions were shown to step far from data (08d).
LIT <- c("litho_gs", paste0("litho_t", 1:6))
edge_frac <- function(cells) vapply(cells, function(cl) {
  adj <- adjacent(r_DI, cl, directions = 8)[1, ]
  v0 <- as.numeric(stk[[LIT]][cl]); vn <- as.matrix(stk[[LIT]][adj])
  mean(apply(vn, 1, function(r) any(abs(r - v0) > 0.2, na.rm = TRUE)), na.rm = TRUE)
}, numeric(1))
set.seed(CFG$seed); rnd_cells <- sample(oc, 500)
grp <- list("sites 1-10" = ev[picked[1:10]], "sites 1-50" = ev[picked], "500 random ocean cells" = rnd_cells)
art <- do.call(rbind, lapply(names(grp), function(g) { cl <- grp[[g]]
  data.frame(group = g, n = length(cl),
             depth_median_m = round(median(values(stk[["depth"]])[cl, 1], na.rm = TRUE)),
             distance_to_land_median_km = round(median(10^values(stk[["dist_coast_km"]])[cl, 1], na.rm = TRUE)),
             lithology_edge_share_pct = round(100 * mean(edge_frac(cl), na.rm = TRUE), 1),
             DInorm_median = round(median(values(r_DI)[cl, 1], na.rm = TRUE) / thr, 3)) }))
write.csv(art, file.path(DIRS$compare, "sampling_priority_artefact_check_toc.csv"), row.names = FALSE)
print(art, row.names = FALSE)

## ---- 8. environmental space (PCA of the DI space) — also the data behind Figure S10b -------------
pc <- prcomp(Ze, center = TRUE, scale. = FALSE)
ve <- round(100 * pc$sdev[1:2]^2 / sum(pc$sdev^2))
env <- data.frame(PC1 = predict(pc, Ze)[, 1], PC2 = predict(pc, Ze)[, 2], DInorm = DI0 / thr,
                  lon = round(xy[, 1], 3), lat = round(xy[, 2], 3))
trn <- data.frame(predict(pc, Zt)[, 1:2])
write.csv(env, gzfile(file.path(DIRS$compare, "environmental_space_sample_toc.csv.gz")), row.names = FALSE)
write.csv(data.frame(variable = rownames(pc$rotation), PC1 = pc$rotation[, 1], PC2 = pc$rotation[, 2],
                     weight = w[rownames(pc$rotation)], variance_pct = c(ve, rep(NA, length(vars) - 2))),
          file.path(DIRS$compare, "environmental_space_loadings_toc.csv"), row.names = FALSE)
site_env <- data.frame(PC1 = predict(pc, Ze[picked, , drop = FALSE])[, 1],
                       PC2 = predict(pc, Ze[picked, , drop = FALSE])[, 2], rank = seq_along(picked))

## ---- 9. figure ----------------------------------------------------------------------------------
PROJ <- "+proj=eqearth +datum=WGS84 +units=m"
DIn  <- r_DI / thr
tp   <- project(DIn, PROJ, res = 25000); dp <- project(DIn, tp, method = "bilinear")
land <- st_transform(ne_countries(scale = 50, returnclass = "sf"), PROJ)
obs_sf  <- st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), PROJ)
site_sf <- st_transform(st_as_sf(sites, coords = c("lon", "lat"), crs = 4326), PROJ)
pa <- ggplot() + geom_spatraster(data = dp) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 2), oob = scales::squish, na.value = NA,
                       name = "DI / AOA threshold") +
  geom_sf(data = land, fill = "grey82", colour = NA, inherit.aes = FALSE) +
  geom_sf(data = obs_sf, size = 0.02, colour = "grey25", alpha = 0.35, inherit.aes = FALSE) +
  geom_sf_text(data = site_sf[site_sf$rank <= 10, ], aes(label = rank), size = 2.1, colour = "#08306B",
               fontface = "bold", nudge_y = 6e5, inherit.aes = FALSE) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = site_sf, aes(fill = rank), size = 1.9, shape = 21, colour = "grey25", stroke = 0.35, inherit.aes = FALSE) +
  scale_fill_gradient(low = "#9ECAE1", high = "#08306B", name = "priority rank", trans = "reverse") +
  coord_sf(crs = PROJ, expand = FALSE, datum = NA) + theme_void(base_size = 8) +
  labs(title = sprintf("a  Where to sample next: %d ranked sites over the ignorance map", N_SITES),
       subtitle = "background: dissimilarity of each ocean cell to the training data; dark grey: observation cells; circles: ranked sites (1 = darkest)") +
  theme(plot.title = element_text(face = "bold", size = 9), legend.position = "right")

ld <- data.frame(var = rownames(pc$rotation), PC1 = pc$rotation[, 1], PC2 = pc$rotation[, 2])
ld <- ld[order(-(ld$PC1^2 + ld$PC2^2)), ][1:6, ]
sc <- 0.8 * max(abs(range(env$PC1, env$PC2))) / max(sqrt(ld$PC1^2 + ld$PC2^2))
pb <- ggplot(env, aes(PC1, PC2)) +
  stat_summary_2d(aes(z = pmin(DInorm, 2)), fun = mean, bins = 80) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 2), name = "mean DI / AOA threshold") +
  geom_point(data = trn, colour = "grey25", size = 0.12, alpha = 0.25) +        # dark grey, as in panel a
  geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1 * sc, yend = PC2 * sc), inherit.aes = FALSE,
               arrow = arrow(length = unit(1.2, "mm")), colour = "grey20", linewidth = 0.3) +
  geom_text(data = ld, aes(PC1 * sc, PC2 * sc, label = var), inherit.aes = FALSE, size = 2, vjust = -0.4, colour = "grey20") +
  geom_point(data = site_env, aes(PC1, PC2, colour = rank), size = 1.4, inherit.aes = FALSE) +
  scale_colour_gradient(low = "#9ECAE1", high = "#08306B", name = "priority rank", trans = "reverse") +
  labs(x = sprintf("environmental axis 1 (%d%%)", ve[1]), y = sprintf("environmental axis 2 (%d%%)", ve[2]),
       title = "b  The same sites in environmental space",
       subtitle = "ocean cells (area-weighted sample) by dissimilarity; dark grey: training cells") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "right")

pc3 <- ggplot(curve, aes(n_sites)) +
  geom_ribbon(aes(ymin = random_min, ymax = random_max), fill = "grey80") +
  geom_line(aes(y = random_mean, colour = "random (mean of 20 sets)")) +
  geom_line(aes(y = uncertainty, colour = "model uncertainty (QRF)")) +
  geom_line(aes(y = novelty, colour = "representative novelty")) +
  scale_colour_manual(values = c("representative novelty" = "#08306B", "model uncertainty (QRF)" = "#C2410C",
                                 "random (mean of 20 sets)" = "grey45"), name = NULL) +
  labs(x = "sites added", y = "mean DI / AOA threshold of the ocean (area-weighted)",
       title = "c  How fast ignorance falls", subtitle = "environmental ignorance only; the error test is panel d") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")

# F5.4 added the criterion this map ranks by (representative novelty) to the same retrospective test
f54 <- file.path(DIRS$compare, "f54_repnovelty_pooled_toc.csv")
h7 <- read.csv(if (file.exists(f54)) f54 else file.path(DIRS$tables, "targeted_sampling_pooled_toc.csv"))
h7 <- h7[h7$strategy != "random", ]
lab <- c(repnovelty = "representative novelty (this map)", novelty = "environmental novelty (DI)",
         uncertainty = "model uncertainty (QRF)", space = "geographic gap filling",
         oracle = "largest current error (upper bound)")
h7$strategy_lab <- factor(lab[h7$strategy], levels = rev(lab))
h7$effort <- factor(paste0("sampled to ", h7$pct_pool, "% of candidates"), levels = paste0("sampled to ", c(40, 60), "% of candidates"))
pd <- ggplot(h7, aes(mean_gain_minus_random, strategy_lab, colour = effort)) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = c("#1D6A73", "#C2410C"), name = NULL) +
  labs(x = "extra RMSE reduction vs random addition of data (log10 TOC; 95% CI)", y = NULL,
       title = "d  Does directed collection reduce error?",
       subtitle = "yes, and most for this map's criterion;
8 regions, 10 random starts each (iteration 2)") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")

fig <- pa / (pb | pc3 | pd) + plot_layout(heights = c(1.05, 1))
ggsave(file.path(DIRS$compare, "fig_sampling_priority_toc.png"), fig, width = 300, height = 210, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(DIRS$compare, "fig_sampling_priority_toc.pdf"), fig, width = 300, height = 210, units = "mm", bg = "white")

## ---- 10. the numbers the text may quote ---------------------------------------------------------
S <- function(m, v, d = "") data.frame(metric = m, value = v, description = d)
out <- rbind(
  S("DI_check_correlation", signif(chk_cor, 8), "reproduced DI vs toc_DI.tif, 2000 cells"),
  S("DI_check_max_abs_error", signif(chk_max, 3), ""),
  S("n_eval_cells", length(ev), "area-weighted ocean sample"),
  S("n_candidates", length(cand), "half highest DI, half random"),
  S("mean_DInorm_start", round(cv_nov[1], 4), "area-weighted mean DI / threshold before new sites"),
  S("mean_DInorm_novelty_10", round(cv_nov[11], 4), ""),
  S(sprintf("mean_DInorm_novelty_%d", N_SITES), round(cv_nov[N_SITES + 1], 4), ""),
  S(sprintf("mean_DInorm_uncertainty_%d", N_SITES), round(cv_unc[N_SITES + 1], 4), ""),
  S(sprintf("mean_DInorm_random_%d", N_SITES), round(curve$random_mean[N_SITES + 1], 4), sprintf("mean of %d random sets", N_RAND)),
  S("pct_reduction_novelty_10", round(red(cv_nov, 10), 2), "relative reduction of mean DI / threshold"),
  S(sprintf("pct_reduction_novelty_%d", N_SITES), round(red(cv_nov, N_SITES), 2), ""),
  S(sprintf("pct_reduction_uncertainty_%d", N_SITES), round(red(cv_unc, N_SITES), 2), ""),
  S(sprintf("pct_reduction_random_%d", N_SITES), round(100 * (curve$random_mean[1] - curve$random_mean[N_SITES + 1]) / curve$random_mean[1], 2), ""),
  S(sprintf("ratio_novelty_over_random_%d", N_SITES), round(red(cv_nov, N_SITES) / (100 * (curve$random_mean[1] - curve$random_mean[N_SITES + 1]) / curve$random_mean[1]), 2), ""),
  S("jaccard_novelty_uncertainty_300km", round(jacc, 3), sprintf("%d of %d sites matched", matched, N_SITES)),
  S("n_sites_flagged_new_drivers", sum(sites$flag_region_needs_new_drivers), "in regions whose decision is 'new drivers or finer resolution'"),
  S("n_sites_flagged_outside_region_support", sum(sites$flag_outside_region_support), "more than 1000 km from any observation"),
  S("n_sites_outside_AOA", sum(!sites$inside_AOA), "sites in cells outside the area of applicability"))
write.csv(out, file.path(DIRS$compare, "sampling_priority_summary_toc.csv"), row.names = FALSE)
print(out, row.names = FALSE)
msg("sampling priority done in %.1f min -> %s", as.numeric(difftime(Sys.time(), t_start, units = "mins")), DIRS$compare)
