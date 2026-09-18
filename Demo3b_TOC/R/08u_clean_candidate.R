# =============================================================================
# 08u_clean_candidate.R — a cleaner iteration 2: fewer modelled layers, fewer artefacts, same knowledge?
# =============================================================================
# Combines the three corrections supported by 08q/08t:
#   texture from sample-based lithology only (no porosity, no grain size; V4 of 08q)
#   distance to land masses >= 25,000 km2 instead of distance to any island (V10 of 08t; removes the SE Pacific wedge)
#   Bio-ORACLE bottom layers smoothed in deep water (V7 of 08q; removes seam steps)
# V12 = the clean candidate. Compared with the full iteration 2 (V0), iteration 1 (V8) and iteration 1 with the corrected
# distance (V11): leave-one-region-out predictions as 07, budget bootstrap as 08k, lithology step ratio far from data
# and bottom-seam step ratio in the SE Pacific (as 08q), and SE Pacific predictions.
Sys.setenv(DEMO3B_ITER = "iter2b_supply_texture")
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(ggplot2); library(patchwork); library(tidyterra) })

m2 <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds")); bt2 <- m2$bestTune
m1 <- readRDS(file.path(ROOT, "outputs", "iter1b_supply", "models", "rf_spatialCV_toc.rds")); bt1 <- m1$bestTune
V2B <- setdiff(names(m2$trainingData), c(".outcome", ".weights")); V1B <- setdiff(names(m1$trainingData), c(".outcome", ".weights"))
BOT <- c("sbt_mean", "sbs_mean", "o2b_mean", "phyc_bot", "sws_bot")
sm <- rast(file.path(CACHE, "predictors_toc_bottom_smoothed_deep_0.1deg.tif"))
dl <- rast(file.path(CACHE, "predictor_dist_large_land_0.1deg.tif"))
stk <- rast(path_toc_stack())
fill_new <- function(D) {
  for (k in BOT) { s <- D[[paste0(k, "_sm")]]; na <- is.na(s); s[na] <- D[[k]][na]; D[[paste0(k, "_sm")]] <- s }
  na <- is.na(D$dist_large_land); D$dist_large_land[na] <- D$dist_coast_km[na]; D }
d <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
d <- fill_new(cbind(d, as.data.frame(extract(sm, cbind(d$x_lon, d$y_lat))), dist_large_land = extract(dl, cbind(d$x_lon, d$y_lat))[, 1]))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster

V12 <- setdiff(V2B, c("porosity", "gs_d50", "gs_d16"))
V12 <- replace(V12, V12 == "dist_coast_km", "dist_large_land")
V12 <- c(setdiff(V12, BOT), paste0(BOT, "_sm"))
VAR <- list(V0_full = V2B, V8_iteration1_no_texture = V1B, V12_clean_candidate = V12)
fit <- function(df, nm, trees = 300) { v <- VAR[[nm]]
  h <- if (nm == "V8_iteration1_no_texture") list(mtry = bt1$mtry, node = bt1$min.node.size) else list(mtry = max(2, round(length(v) * bt2$mtry / length(V2B))), node = bt2$min.node.size)
  ranger(x = df[, v, drop = FALSE], y = df$y, num.trees = trees, mtry = min(h$mtry, length(v)), min.node.size = h$node,
         splitrule = "variance", num.threads = CFG$cores, seed = CFG$seed) }

# --- leave-one-region-out and budget ------------------------------------------------------------------------------------
P <- readRDS(file.path(DIRS$models, "distance_to_land_lro.rds"))
stopifnot(identical(P$x_lon, d$x_lon), identical(P$region, d$region))
p <- rep(NA_real_, nrow(d))
for (r in sort(unique(d$region))) { tr <- d$region != r; p[!tr] <- predict(fit(d[tr, ], "V12_clean_candidate"), d[!tr, V12, drop = FALSE], num.threads = CFG$cores)$predictions }
P$V12_clean_candidate <- p
msg("V12 clean candidate (%d predictors): withheld RMSE %.4f", length(V12), sqrt(mean((p - d$y)^2)))
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
NM <- c("V0_full", "V8_iteration1_no_texture", "V11_iteration1_distance_to_large_land", "V12_clean_candidate")
E <- sapply(NM, function(nm) P[[nm]] - P$y)
BS <- setNames(lapply(NM, function(nm) t(vapply(IDX, function(i) budget(E[i, nm], P$region[i], P$lith[i]), numeric(4)))), NM)
cmp <- function(a, b) { ea <- budget(E[, a], P$region, P$lith); eb <- budget(E[, b], P$region, P$lith)
  ci <- apply(BS[[b]] - BS[[a]], 2, quantile, c(0.025, 0.975))
  data.frame(reference = a, variant = b, component = names(ea), value_reference = ea, value_variant = eb, difference = eb - ea,
             ci_lo = ci[1, ], ci_hi = ci[2, ], pct_change = 100 * (eb - ea) / ea, ci_excludes_zero = ci[1, ] > 0 | ci[2, ] < 0) }
bud <- rbind(cmp("V0_full", "V12_clean_candidate"), cmp("V8_iteration1_no_texture", "V12_clean_candidate"),
             cmp("V11_iteration1_distance_to_large_land", "V12_clean_candidate"))
write.csv(bud, file.path(DIRS$compare, "clean_candidate_budget.csv"), row.names = FALSE)
print(bud, digits = 3, row.names = FALSE)
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
reg <- do.call(rbind, lapply(sort(unique(P$region)), function(r) { j <- P$region == r; data.frame(region = RN[r], n = sum(j), t(sqrt(colMeans(E[j, , drop = FALSE]^2)))) }))
write.csv(reg, file.path(DIRS$compare, "clean_candidate_regions.csv"), row.names = FALSE)
print(reg, digits = 3, row.names = FALSE)
saveRDS(P, file.path(DIRS$models, "clean_candidate_lro.rds"))

# --- artefacts (same sampled pairs and box as 08q) ------------------------------------------------------------------------
dep <- stk[["depth"]]; dom <- which.max(stk[[paste0("litho_t", 1:6)]])
dkm <- rast(file.path(DIRS$rasters, "toc_distance_nearest_observation_km.tif"))
Dm <- as.matrix(dom, wide = TRUE); Km <- as.matrix(dkm, wide = TRUE); Zm <- as.matrix(dep, wide = TRUE); nr <- nrow(Dm); nc <- ncol(Dm)
set.seed(CFG$seed + 7)
kh <- sample.int(nr * (nc - 1), 1.5e6); kv <- sample.int(nr * nc, 1.5e6); kv <- kv[kv %% nr != 0]
pr <- rbind(data.frame(k1 = kh, k2 = kh + nr), data.frame(k1 = kv, k2 = kv + 1))
pr$cross <- Dm[pr$k1] != Dm[pr$k2]; pr$km <- pmax(Km[pr$k1], Km[pr$k2]); pr$z <- pmin(Zm[pr$k1], Zm[pr$k2])
pr <- pr[complete.cases(pr), ]; rm(Dm, Km, Zm); invisible(gc())
pr$deep_far <- pr$z > 1500 & pr$km > 500
set.seed(CFG$seed + 8)
i_df <- which(pr$deep_far); i_ot <- which(!pr$deep_far)
pr <- pr[c(if (length(i_df) > 150000) sample(i_df, 150000) else i_df, sample(i_ot, min(150000, length(i_ot)))), ]
lin2cell <- function(k) { r <- (k - 1) %% nr + 1; cc <- (k - 1) %/% nr + 1; (r - 1) * nc + cc }
pr$c1 <- lin2cell(pr$k1); pr$c2 <- lin2cell(pr$k2)
eB <- ext(-140, -90, -50, 0)
cb <- cells(stk[[1]], eB); if (is.matrix(cb)) cb <- cb[, 1]
rcb <- rowColFromCell(stk, cb); hb <- cb[rcb[, 2] < max(rcb[, 2])]; sp <- data.frame(c1 = hb, c2 = hb + 1)
allc <- unique(c(pr$c1, pr$c2, cb, sp$c2))
ND <- fill_new(cbind(as.data.frame(extract(stk, allc)), as.data.frame(extract(sm, allc)), dist_large_land = extract(dl, allc)[, 1]))
PRED <- lapply(c(V0_full = "V0_full", V12_clean_candidate = "V12_clean_candidate"), function(nm) {
  ok <- complete.cases(ND[, VAR[[nm]]]); q <- rep(NA_real_, nrow(ND))
  q[ok] <- predict(fit(d, nm, trees = 500), ND[ok, VAR[[nm]], drop = FALSE], num.threads = CFG$cores)$predictions; q })
m1i <- match(pr$c1, allc); m2i <- match(pr$c2, allc); s1 <- match(sp$c1, allc); s2 <- match(sp$c2, allc)
nstep <- function(k) { a <- abs(ND[[k]][s1] - ND[[k]][s2]); a / median(a[a > 0], na.rm = TRUE) }
seam <- pmax(nstep("sbt_mean"), nstep("sbs_mean"), nstep("o2b_mean"), na.rm = TRUE) > 10
ratio <- function(x, across, s) mean(x[s & across], na.rm = TRUE) / mean(x[s & !across], na.rm = TRUE)
art <- do.call(rbind, lapply(names(PRED), function(nm) { q <- PRED[[nm]]; a <- abs(q[m1i] - q[m2i]); b <- abs(q[s1] - q[s2])
  data.frame(variant = nm, step_ratio_lithology_deep_far = ratio(a, pr$cross, pr$deep_far),
             step_ratio_bottom_seams_SE_Pacific = ratio(b, seam, rep(TRUE, length(b)))) }))
write.csv(art, file.path(DIRS$compare, "clean_candidate_artefacts.csv"), row.names = FALSE)
print(art, digits = 3, row.names = FALSE)
tmpl <- init(stk[[1]], NA)
zr <- function(q) { r <- tmpl; r[cb] <- q[match(cb, allc)]; crop(r, eB) }
pz <- function(r, title, sc) ggplot() + geom_spatraster(data = r) + sc + coord_sf(expand = FALSE) + labs(title = title, fill = NULL) +
  theme_minimal(base_size = 6.5) + theme(axis.text = element_blank(), plot.title = element_text(size = 7, face = "bold"), legend.key.width = unit(2, "mm"))
sct <- scale_fill_viridis_c(limits = c(-1.1, 0.3), oob = scales::squish, na.value = "grey85")
fig <- pz(zr(PRED$V0_full), "full iteration 2 (log10 TOC)", sct) | pz(zr(PRED$V12_clean_candidate), "clean candidate (log10 TOC)", sct) |
  pz(zr(PRED$V12_clean_candidate - PRED$V0_full), "clean minus full",
     scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-0.3, 0.3), oob = scales::squish, na.value = "grey85"))
ggsave(file.path(DIRS$compare, "fig_clean_candidate_SEPacific.png"), fig, width = 220, height = 80, units = "mm", dpi = 250, bg = "white")
msg("clean candidate done")
