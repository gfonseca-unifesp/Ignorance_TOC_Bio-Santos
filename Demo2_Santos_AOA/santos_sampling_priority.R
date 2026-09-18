## ================================================================
## santos_sampling_priority.R — F3 of ROADMAP_MSv6: where, and when, to sample next in the Santos Basin
## ================================================================
## Run from the project root (needs .Renviron or IGNORANCE_ROOT):
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_sampling_priority.R
##
## The audit (santos_audit*.R) gave three results that this map must respect, and must not blur:
##   T3/T3b  with ~16 starting stations, no guide (novelty, uncertainty) beat random stations;
##   T4b     at equal effort, a new time interval pays about twice as much as new stations
##           (0.077 against 0.033; new stations in a new interval 0.083);
##   v2 map  18.6% of the basin is an environmental extrapolation (DI / threshold > 1.2).
## So the map carries TWO layers with DIFFERENT objectives, and the text must keep them apart:
##   layer A — new stations inside the stratified design, drawn AT RANDOM and NOT ranked. This is the
##             strategy the audit supports for reducing error at this level of ignorance.
##   layer B — ranked sites whose objective is to shrink the extrapolated area, i.e. to extend the
##             applicability of the published models. This objective was NEVER tested against error:
##             it is a design recommendation, not a result.
## Time is not a point on the map: the recommendation is a new survey covering the zones, with the
## existing stations plus layer A (reported in the figure footnote with the T4b gains).
##
## Nothing is overwritten: all outputs are new files in Demo2_Santos_AOA/outputs/.
suppressPackageStartupMessages({
  library(CAST); library(caret); library(randomForest); library(FNN)
  library(terra); library(sf); library(ggplot2); library(patchwork); library(rnaturalearth)
})
SEED <- 20260914                     # the seed of santos_aoa.R (applicability and maps)
set.seed(SEED)
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs")
msg  <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")

N_A     <- as.integer(Sys.getenv("SANTOS_N_A", "30"))    # layer A: random stations within the design
N_B     <- as.integer(Sys.getenv("SANTOS_N_B", "20"))    # layer B: ranked sites that extend applicability
M_CAND  <- 2000                                          # layer B candidates evaluated per step
K_NB    <- 250                                           # neighbours used to approximate a candidate's gain
N_RAND  <- 20                                            # random site sets for the curve
SEP_A   <- 20                                            # km: environmental change between surveys equals ~20 km (T4c)
EDGE    <- c(0.8, 1.2)                                   # soft AOA edge of the v2 map
ZONES   <- list("shelf (<= 200 m)" = c(0, 200), "slope (200-1000 m)" = c(200, 1000), "deep (1000-2400 m)" = c(1000, 2400))
ENV_PRED <- c("Long", "Lat", "Camp.1", "Depth")
t_start <- Sys.time()

## ---- 1. the published environmental models, and their applicability -------------------
sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
G   <- sd1$Geo_mean                       # 198 station rows: Long, Lat, Camp.1, Depth
GB  <- sd1$Grid_bat                       # grid: Long, Lat, Depth, Camp.1 (two campaigns)
rfG <- attr(G, "rf"); rm(sd1); invisible(gc())
env_models <- setNames(lapply(rfG, `[[`, "m"), sub("~Geo_mean$", "", names(rfG)))
msg("environmental models: %d | grid rows: %d (%d cells x 2 campaigns)", length(env_models), nrow(GB),
    nrow(unique(GB[, c("Long", "Lat")])))

cv_sets <- function(m) {                  # station-grouped folds, exactly as santos_aoa.R
  rn <- rownames(m$trainingData); stopifnot(all(rn %in% rownames(G)))
  st <- paste(G[rn, "Long"], G[rn, "Lat"])
  f  <- CreateSpacetimeFolds(data.frame(st = st), spacevar = "st", k = 5, seed = SEED)
  list(CVtest = f$indexOut, CVtrain = f$index)
}
cp2 <- which(GB$Camp.1 == 2)              # the most recent campaign
grid <- GB[cp2, ENV_PRED]
tdis <- lapply(env_models, function(m) { cvs <- cv_sets(m)
  suppressMessages(trainDI(m, CVtest = cvs$CVtest, CVtrain = cvs$CVtrain, verbose = FALSE)) })
msg("trainDI rebuilt for %d environmental models (station-grouped folds)", length(tdis))

# Every model uses the same four predictors, so a model only differs by its centring, scaling, weights,
# mean training distance and threshold. Distances are therefore a per-model rescaling of the raw
# predictor differences: a_mj = w_mj / (scale_mj * trainDist_avrgmean_m * threshold_m), and the
# normalised dissimilarity DI/threshold is the Euclidean distance in the space scaled by a_m.
coef_of <- function(tdi) {
  v <- tdi$variables; w <- unlist(tdi$weight)[v]
  a <- w / (tdi$scaleparam[["scaled:scale"]][v] * tdi$trainDist_avrgmean * tdi$threshold)
  out <- setNames(numeric(length(ENV_PRED)), ENV_PRED); out[v] <- a; out
}
A  <- t(vapply(tdis, coef_of, numeric(length(ENV_PRED))))      # 44 x 4
A2 <- A^2
Xg <- as.matrix(grid[, ENV_PRED])                             # grid cells x 4
Xt <- lapply(tdis, function(tdi) as.matrix(tdi$train[, ENV_PRED, drop = FALSE]))

# normalised DI of every grid cell for every model: nearest training station in that model's space
dinorm_grid <- function(i) {
  s <- A[i, ]
  get.knnx(sweep(Xt[[i]], 2, s, "*"), sweep(Xg, 2, s, "*"), k = 1)$nn.dist[, 1]
}
DIN <- vapply(seq_len(nrow(A)), dinorm_grid, numeric(nrow(Xg)))   # cells x 44, already DI / threshold
DIN0 <- DIN
dn_mean <- rowMeans(DIN)

## ---- 2. mandatory check against the published consensus ------------------------------
cons <- read.csv(file.path(OUT, "santos_aoa_consensus.csv"))
cons2 <- cons[cons$Camp == 2, ]
stopifnot(nrow(cons2) == nrow(Xg), max(abs(cons2$Long - grid$Long)) < 1e-4, max(abs(cons2$Lat - grid$Lat)) < 1e-4)   # the CSV is written with 7 significant digits
chk_cor <- cor(cons2$env_DInorm_mean, dn_mean); chk_max <- max(abs(cons2$env_DInorm_mean - dn_mean))
msg("check vs santos_aoa_consensus.csv (Camp 2): correlation %.8f, maximum absolute error %.3g", chk_cor, chk_max)
if (!(chk_cor > 0.999)) stop("the reproduced environmental DI does not match the published one")

pct_extrap <- function(v) 100 * mean(v > EDGE[2])
msg("basin in environmental extrapolation (Camp 2): %.2f%% | in transition: %.2f%%",
    pct_extrap(dn_mean), 100 * mean(dn_mean > EDGE[1] & dn_mean <= EDGE[2]))

## ---- 3. layer A — random stations within the stratified design -----------------------
stations <- unique(G[, c("Long", "Lat", "Depth")])
km_between <- function(a, b) {                      # great-circle distance, km (small basin, spherical law)
  r <- 6371; p <- pi / 180
  outer(seq_len(nrow(a)), seq_len(nrow(b)), function(i, j)
    r * acos(pmin(1, sin(a[i, 2] * p) * sin(b[j, 2] * p) + cos(a[i, 2] * p) * cos(b[j, 2] * p) * cos((a[i, 1] - b[j, 1]) * p))))
}
elig <- which(dn_mean <= EDGE[2] & grid$Depth <= 2400 & grid$Depth > 0)
set.seed(SEED)
A_rows <- list()
for (z in names(ZONES)) {
  lim <- ZONES[[z]]
  pool <- elig[grid$Depth[elig] > lim[1] & grid$Depth[elig] <= lim[2]]
  pool <- sample(pool)                              # random order: layer A is NOT ranked
  taken <- as.matrix(stations[, c("Long", "Lat")]); got <- integer(0)
  for (i in pool) {
    if (length(got) == N_A %/% length(ZONES)) break
    p <- as.matrix(grid[i, c("Long", "Lat")])
    if (min(km_between(p, taken)) < SEP_A) next
    got <- c(got, i); taken <- rbind(taken, p)
  }
  code <- c("S", "L", "D")[match(z, names(ZONES))]        # shelf / slope / deep
  A_rows[[z]] <- data.frame(station = sprintf("A%s%d", code, seq_along(got)), zone = z,
                            lon = round(grid$Long[got], 4), lat = round(grid$Lat[got], 4),
                            depth_m = round(grid$Depth[got]), DInorm_env = round(dn_mean[got], 3))
  msg("layer A: %s -> %d stations (pool %d cells)", z, length(got), length(pool))
}
layerA <- do.call(rbind, A_rows); rownames(layerA) <- NULL
write.csv(layerA, file.path(OUT, "santos_priority_A_stations.csv"), row.names = FALSE)

## ---- 4. layer B — ranked sites that extend applicability ------------------------------
# objective: the smooth area above the transition edge, sum over cells of max(0, mean DI/threshold - 0.8)
objective <- function(M) sum(pmax(rowMeans(M) - EDGE[1], 0))
upd_all <- function(M, i) {                          # exact update of all 44 DIs over the whole grid
  d2 <- sweep(Xg, 2, Xg[i, ], "-")^2
  pmin(M, sqrt(d2 %*% t(A2)))
}
cand_pool <- which(dn_mean > EDGE[1])
Sxyz <- scale(Xg[, c("Long", "Lat", "Depth")])       # geographic + depth space, scaled (used only to pick neighbours)
msg("layer B candidate pool: %d cells with DI/threshold > %.1f", length(cand_pool), EDGE[1])
nbrs <- NULL
pickB <- integer(0); gainB <- numeric(0); dinB <- numeric(0)
for (k in seq_len(N_B)) {
  set.seed(SEED + k)
  cand <- sample(cand_pool, min(M_CAND, length(cand_pool)))
  nb <- get.knnx(Sxyz, Sxyz[cand, , drop = FALSE], k = K_NB)$nn.index
  base <- rowMeans(DIN)
  gains <- vapply(seq_along(cand), function(j) {
    ci <- nb[j, ]
    d2 <- sweep(Xg[ci, , drop = FALSE], 2, Xg[cand[j], ], "-")^2
    newDIN <- pmin(DIN[ci, , drop = FALSE], sqrt(d2 %*% t(A2)))
    sum(pmax(base[ci] - EDGE[1], 0)) - sum(pmax(rowMeans(newDIN) - EDGE[1], 0))
  }, numeric(1))
  b <- cand[which.max(gains)]
  pickB <- c(pickB, b); gainB <- c(gainB, max(gains)); dinB <- c(dinB, base[b])
  DIN <- upd_all(DIN, b)
  if (k %% 5 == 0) msg("  layer B: %d sites; extrapolated area %.2f%% (start %.2f%%)", k,
                       pct_extrap(rowMeans(DIN)), pct_extrap(dn_mean))
}
DIN_B <- DIN

## ---- 5. curve: extrapolated area against the number of layer-B sites ------------------
curve_of <- function(idx) {
  M <- DIN0; out <- numeric(length(idx) + 1); out[1] <- pct_extrap(rowMeans(M))
  for (k in seq_along(idx)) { M <- upd_all(M, idx[k]); out[k + 1] <- pct_extrap(rowMeans(M)) }
  out
}
cv_greedy <- curve_of(pickB)
set.seed(SEED)
cv_rand <- vapply(seq_len(N_RAND), function(i) curve_of(sample(cand_pool, N_B)), numeric(N_B + 1))
curveB <- data.frame(n_sites = 0:N_B, greedy = cv_greedy, random_mean = rowMeans(cv_rand),
                     random_min = apply(cv_rand, 1, min), random_max = apply(cv_rand, 1, max))
write.csv(curveB, file.path(OUT, "santos_priority_curve.csv"), row.names = FALSE)
msg("extrapolated area: %.2f%% -> %.2f%% with %d ranked sites (random sites: %.2f%%)",
    cv_greedy[1], cv_greedy[N_B + 1], N_B, curveB$random_mean[N_B + 1])

layerB <- data.frame(rank = seq_along(pickB), lon = round(grid$Long[pickB], 4), lat = round(grid$Lat[pickB], 4),
                     depth_m = round(grid$Depth[pickB]), DInorm_initial = round(dn_mean[pickB], 3),
                     DInorm_at_selection = round(dinB, 3), gain_smooth_area = signif(gainB, 4),
                     pct_extrapolation_after = round(cv_greedy[-1], 3))
write.csv(layerB, file.path(OUT, "santos_priority_B_sites.csv"), row.names = FALSE)
print(head(layerB, 10), row.names = FALSE)

## ---- 6. figure -------------------------------------------------------------------------
RES <- 0.02
tmpl <- rast(ext(min(grid$Long) - RES, max(grid$Long) + RES, min(grid$Lat) - RES, max(grid$Lat) + RES),
             resolution = RES, crs = "EPSG:4326")
r_din <- rasterize(as.matrix(grid[, c("Long", "Lat")]), tmpl, values = dn_mean, fun = mean)
r_din <- focal(r_din, 3, "mean", na.policy = "only", na.rm = TRUE)      # fill isolated empty cells only
df <- as.data.frame(r_din, xy = TRUE, na.rm = TRUE); names(df)[3] <- "DInorm"
land <- ne_countries(scale = 10, country = "Brazil", returnclass = "sf")
labB <- layerB[layerB$rank <= 10, ]
T4B <- c(new_stations_same_survey = 0.033, same_stations_other_survey = 0.077, new_stations_other_survey = 0.083)

pa <- ggplot() +
  geom_raster(data = df, aes(x, y, fill = pmin(DInorm, 2))) +
  scale_fill_gradient(low = "#F4F6F8", high = "#3E5C76", limits = c(0, 2), oob = scales::squish,
                      name = "Environmental tier:\ndissimilarity to sampled\nconditions (DI / threshold)") +
  geom_contour(data = df, aes(x, y, z = DInorm), breaks = EDGE, colour = "grey25", linewidth = 0.25, linetype = "dashed") +
  geom_sf(data = land, fill = "#E6DFCF", colour = "grey45", linewidth = 0.2) +
  geom_point(data = stations, aes(Long, Lat), colour = "grey45", size = 1.1) +
  geom_point(data = transform(layerA, zone = factor(zone, levels = names(ZONES))),
             aes(lon, lat, shape = zone), colour = "#1D6A73", size = 2, stroke = 0.6, fill = NA) +
  scale_shape_manual(values = c("shelf (<= 200 m)" = 25, "slope (200-1000 m)" = 23, "deep (1000-2400 m)" = 24),
                     breaks = names(ZONES), name = "Layer A — random stations within\nthe stratified design (T3/T3b:\nno guide beat random)") +
  geom_point(data = layerB, aes(lon, lat, colour = rank), size = 2.2) +
  scale_colour_gradient(low = "#FCA5A5", high = "#7F1D1D", trans = "reverse",     # rank 1 = darkest
                        name = "Layer B — ranked sites that\nextend applicability (untested\nagainst error)") +
  geom_text(data = labB, aes(lon, lat, label = rank), size = 2.1, colour = "black", fontface = "bold", nudge_y = 0.12) +
  coord_sf(xlim = range(grid$Long), ylim = range(grid$Lat), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude",
       title = "Santos Basin: where to add stations, and where the models stop being applicable",
       subtitle = sprintf(paste0("Dashed contours: transition (%.1f) and extrapolation (%.1f) edges; grey dots: the 99 existing stations.\n",
                                 "When to sample: at equal effort a new time interval pays about twice as much as new stations\n",
                                 "(T4b gains %.3f same survey, %.3f other survey, %.3f both), so layers A and B are meant to be\n",
                                 "occupied in a NEW survey, together with the existing stations."),
                          EDGE[1], EDGE[2], T4B[1], T4B[2], T4B[3])) +
  theme_bw(base_size = 8) + theme(legend.position = "right", panel.grid = element_blank(),
                                  plot.title = element_text(face = "bold"))

pb <- ggplot(curveB, aes(n_sites)) +
  geom_ribbon(aes(ymin = random_min, ymax = random_max), fill = "grey85") +
  geom_line(aes(y = random_mean, colour = "random sites among the same candidates")) +
  geom_line(aes(y = greedy, colour = "layer B (ranked)")) +
  scale_colour_manual(values = c("layer B (ranked)" = "#7F1D1D", "random sites among the same candidates" = "grey45"), name = NULL) +
  labs(x = "layer-B sites added", y = "% of the basin in environmental extrapolation",
       title = "How much extrapolation each added site removes",
       subtitle = "objective: applicability of the published models, not prediction error") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")

fig <- pa / pb + plot_layout(heights = c(2.1, 1))
ggsave(file.path(OUT, "fig_santos_sampling_priority.png"), fig, width = 210, height = 240, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_santos_sampling_priority.pdf"), fig, width = 210, height = 240, units = "mm", bg = "white")

## ---- 7. the numbers the text may quote --------------------------------------------------
S <- function(m, v, d = "") data.frame(metric = m, value = v, description = d)
summ <- rbind(
  S("check_correlation_env_DInorm", signif(chk_cor, 8), "reproduced vs santos_aoa_consensus.csv, Camp 2"),
  S("check_max_abs_error", signif(chk_max, 3), ""),
  S("pct_extrapolation_camp2_start", round(cv_greedy[1], 2), "mean DI/threshold > 1.2, Camp 2 grid"),
  S("pct_extrapolation_after_B", round(cv_greedy[N_B + 1], 2), sprintf("%d ranked sites", N_B)),
  S("pct_extrapolation_after_random", round(curveB$random_mean[N_B + 1], 2), sprintf("mean of %d random sets", N_RAND)),
  S("n_layerA_stations", nrow(layerA), "random stations within the stratified design"),
  S("n_layerB_sites", nrow(layerB), "ranked sites"),
  S("min_separation_layerA_km", SEP_A, "T4c: environmental change between surveys equals ~20 km"),
  S("T4b_gain_new_stations_same_survey", T4B[1], ""),
  S("T4b_gain_same_stations_other_survey", T4B[2], ""),
  S("T4b_gain_new_stations_other_survey", T4B[3], ""))
write.csv(summ, file.path(OUT, "santos_priority_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
writeLines(capture.output(sessioninfo::session_info()), file.path(OUT, "sessionInfo_sampling_priority.txt"))
msg("done in %.1f min -> %s", as.numeric(difftime(Sys.time(), t_start, units = "mins")), OUT)
