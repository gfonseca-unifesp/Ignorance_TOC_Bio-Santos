# =============================================================================
# 08q_layer_sensitivity.R — do the conclusions depend on integrated layers of uncertain provenance?
# =============================================================================
# Every integrated layer is itself a product: porosity (machine-learning map, Martin et al. 2015), grain size
# (gridded NGDC samples), lithology (sample-based classification, Dutkiewicz et al. 2015), satellite surface POC
# and its Martin-curve flux, and Bio-ORACLE v3 bottom layers (ocean-model output interpolated to the seafloor,
# with seams; 08p). Sensitivity refits of iteration 2 (26 predictors):
#   V0 full                        V1 without porosity             V2 without grain size (D50, D16)
#   V3 without lithology (8)       V4 without porosity and grain size (texture from sample-based lithology only)
#   V5 without satellite supply (poc_surf, poc_flux)               V6 without near-bottom phytoplankton
#   V7 bottom layers smoothed in deep water (0.5-deg mean where depth > 1500 m: sbt, sbs, o2b, phyc, sws)
#   V8 iteration 1 (16 predictors, no texture): the reference for the texture gain
# For each: leave-one-region-out predictions exactly as 07_calibration.R (ranger, 300 trees, same seed, same regions;
# mtry = the tuned share of predictors, min.node.size of the tuned model), the ignorance budget, and a stratified
# block bootstrap of the difference to V0 and to V8 (1000 resamples, 300-km blocks within regions; as 08k).
# Artefacts, from full-data models (500 trees) predicted on sampled cells:
#   step ratio across lithology boundaries in deep water > 500 km from data (as 08o), and
#   step ratio at bottom-layer seams in the SE Pacific box (neighbour pairs where bottom T, S or O2 jump > 10x their
#   median cell-to-cell change; as 08p).
Sys.setenv(DEMO3B_ITER = "iter2b_supply_texture")
source("R/00_config.R")
suppressPackageStartupMessages({ library(ranger); library(ggplot2); library(patchwork); library(tidyterra) })

m2  <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds")); bt2 <- m2$bestTune
m1  <- readRDS(file.path(ROOT, "outputs", "iter1b_supply", "models", "rf_spatialCV_toc.rds")); bt1 <- m1$bestTune
V2B <- setdiff(names(m2$trainingData), c(".outcome", ".weights"))
V1B <- setdiff(names(m1$trainingData), c(".outcome", ".weights"))
d   <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
stopifnot(all(c(V1B, V2B) %in% names(d)))
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster   # same as 07

# --- smoothed bottom layers (deep water only) ------------------------------------------------------------------------
BOT <- c("sbt_mean", "sbs_mean", "o2b_mean", "phyc_bot", "sws_bot")
stk <- rast(path_toc_stack()); dep <- stk[["depth"]]
f_sm <- file.path(CACHE, "predictors_toc_bottom_smoothed_deep_0.1deg.tif")
if (!file.exists(f_sm)) {
  sm <- rast(lapply(BOT, function(k) { r <- stk[[k]]; ifel(dep > 1500, focal(r, w = 5, fun = "mean", na.rm = TRUE), r) }))
  names(sm) <- paste0(BOT, "_sm")
  writeRaster(sm, f_sm, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  msg("smoothed bottom layers written")
}
sm <- rast(f_sm)
fill_sm <- function(D) { for (k in BOT) { s <- D[[paste0(k, "_sm")]]; na <- is.na(s); s[na] <- D[[k]][na]; D[[paste0(k, "_sm")]] <- s }; D }
d <- fill_sm(cbind(d, as.data.frame(extract(sm, cbind(d$x_lon, d$y_lat)))))

POR <- "porosity"; GS <- c("gs_d50", "gs_d16"); LITH <- c("litho_gs", paste0("litho_t", 1:6)); SAT <- c("poc_surf", "poc_flux")
VAR <- list(
  V0_full                    = V2B,
  V1_no_porosity             = setdiff(V2B, POR),
  V2_no_grain_size           = setdiff(V2B, GS),
  V3_no_lithology            = setdiff(V2B, LITH),
  V4_lithology_only_texture  = setdiff(V2B, c(POR, GS)),
  V5_no_satellite_supply     = setdiff(V2B, SAT),
  V6_no_bottom_phytoplankton = setdiff(V2B, "phyc_bot"),
  V7_bottom_layers_smoothed  = c(setdiff(V2B, BOT), paste0(BOT, "_sm")),
  V8_iteration1_no_texture   = V1B)
hp <- function(nm) {
  if (nm == "V8_iteration1_no_texture") return(list(mtry = bt1$mtry, node = bt1$min.node.size))
  list(mtry = max(2, round(length(VAR[[nm]]) * bt2$mtry / length(V2B))), node = bt2$min.node.size)
}
fit <- function(df, nm, trees = 300) { v <- VAR[[nm]]; h <- hp(nm)
  ranger(x = df[, v, drop = FALSE], y = df$y, num.trees = trees, mtry = min(h$mtry, length(v)), min.node.size = h$node,
         splitrule = "variance", num.threads = CFG$cores, seed = CFG$seed) }

# --- (1) leave-one-region-out ----------------------------------------------------------------------------------------
f_lro <- file.path(DIRS$models, "layer_sensitivity_lro.rds")
if (!file.exists(f_lro)) {
  P <- d[, c("x_lon", "y_lat", "y", "region")]
  for (nm in names(VAR)) {
    t0 <- Sys.time(); p <- rep(NA_real_, nrow(d))
    for (r in sort(unique(d$region))) { tr <- d$region != r; p[!tr] <- predict(fit(d[tr, ], nm), d[!tr, VAR[[nm]], drop = FALSE], num.threads = CFG$cores)$predictions }
    P[[nm]] <- p
    msg("%s: %d predictors, withheld RMSE %.4f (%.1f min)", nm, length(VAR[[nm]]), sqrt(mean((p - d$y)^2)), as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  saveRDS(P, f_lro)
}
P <- readRDS(f_lro)
key <- function(lon, lat) paste(round(lon, 3), round(lat, 3))
chk <- function(it, nm) { o <- read.csv(file.path(ROOT, "outputs", it, "tables", "calibration_outer_points_toc.csv"))
  z <- merge(data.frame(k = key(P$x_lon, P$y_lat), a = P[[nm]]), data.frame(k = key(o$lon, o$lat), b = o$pred), by = "k")
  msg("reproduction %s vs %s (07): n = %d, r = %.4f, max |diff| = %.4f", nm, it, nrow(z), cor(z$a, z$b), max(abs(z$a - z$b))) }
chk("iter2b_supply_texture", "V0_full"); chk("iter1b_supply", "V8_iteration1_no_texture")

# --- (2) budget and bootstrap ---------------------------------------------------------------------------------------
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
E <- sapply(names(VAR), function(nm) P[[nm]] - P$y)
BS <- lapply(names(VAR), function(nm) t(vapply(IDX, function(i) budget(E[i, nm], P$region[i], P$lith[i]), numeric(4))))   # same resamples for all variants
names(BS) <- names(VAR)
cmp <- function(a, b) { ea <- budget(E[, a], P$region, P$lith); eb <- budget(E[, b], P$region, P$lith)
  ci <- apply(BS[[b]] - BS[[a]], 2, quantile, c(0.025, 0.975))
  data.frame(reference = a, variant = b, component = names(ea), value_reference = ea, value_variant = eb, difference = eb - ea,
             ci_lo = ci[1, ], ci_hi = ci[2, ], pct_change = 100 * (eb - ea) / ea, ci_excludes_zero = ci[1, ] > 0 | ci[2, ] < 0) }
bud <- rbind(do.call(rbind, lapply(setdiff(names(VAR), "V0_full"), function(b) cmp("V0_full", b))),
             do.call(rbind, lapply(setdiff(names(VAR), "V8_iteration1_no_texture"), function(b) cmp("V8_iteration1_no_texture", b))))
write.csv(bud, file.path(DIRS$compare, "layer_sensitivity_budget.csv"), row.names = FALSE)
print(bud[bud$component == "total", ], digits = 3, row.names = FALSE)
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
reg <- do.call(rbind, lapply(sort(unique(P$region)), function(r) { j <- P$region == r
  data.frame(region = RN[r], n = sum(j), t(sqrt(colMeans(E[j, , drop = FALSE]^2)))) }))
write.csv(reg, file.path(DIRS$compare, "layer_sensitivity_regions.csv"), row.names = FALSE)
print(reg, digits = 3, row.names = FALSE)

# --- (3) artefacts ----------------------------------------------------------------------------------------------------
dom <- which.max(stk[[paste0("litho_t", 1:6)]])
dkm <- rast(file.path(DIRS$rasters, "toc_distance_nearest_observation_km.tif"))
Dm <- as.matrix(dom, wide = TRUE); Km <- as.matrix(dkm, wide = TRUE); Zm <- as.matrix(dep, wide = TRUE)
nr <- nrow(Dm); nc <- ncol(Dm)
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
cb <- cells(stk[[1]], eB); if (is.matrix(cb)) cb <- cb[, 1]   # vector in recent terra versions
rcb <- rowColFromCell(stk, cb)
hb <- cb[rcb[, 2] < max(rcb[, 2])]; sp <- data.frame(c1 = hb, c2 = hb + 1)

allc <- unique(c(pr$c1, pr$c2, cb, sp$c2))
ND <- fill_sm(cbind(as.data.frame(extract(stk, allc)), as.data.frame(extract(sm, allc))))
f_pred <- file.path(DIRS$models, "layer_sensitivity_cellpred.rds")
if (!file.exists(f_pred)) {
  PRED <- list()
  for (nm in names(VAR)) {
    t0 <- Sys.time(); mfull <- fit(d, nm, trees = 500); ok <- complete.cases(ND[, VAR[[nm]]]); p <- rep(NA_real_, nrow(ND))
    p[ok] <- predict(mfull, ND[ok, VAR[[nm]], drop = FALSE], num.threads = CFG$cores)$predictions
    PRED[[nm]] <- p; msg("cell predictions %s (%.1f min)", nm, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  saveRDS(PRED, f_pred)
}
PRED <- readRDS(f_pred)
m1i <- match(pr$c1, allc); m2i <- match(pr$c2, allc)
s1 <- match(sp$c1, allc); s2 <- match(sp$c2, allc)
nstep <- function(k) { a <- abs(ND[[k]][s1] - ND[[k]][s2]); a / median(a[a > 0], na.rm = TRUE) }
seam <- pmax(nstep("sbt_mean"), nstep("sbs_mean"), nstep("o2b_mean"), na.rm = TRUE) > 10
art <- do.call(rbind, lapply(names(VAR), function(nm) { p <- PRED[[nm]]
  a <- abs(p[m1i] - p[m2i]); b <- abs(p[s1] - p[s2])
  ratio <- function(x, across, s) mean(x[s & across], na.rm = TRUE) / mean(x[s & !across], na.rm = TRUE)
  data.frame(variant = nm, step_ratio_lithology_deep_far = ratio(a, pr$cross, pr$deep_far),
             step_ratio_lithology_other_cells = ratio(a, pr$cross, !pr$deep_far),
             step_ratio_bottom_seams_SE_Pacific = ratio(b, seam, rep(TRUE, length(b))), seam_pairs = sum(seam, na.rm = TRUE)) }))
write.csv(art, file.path(DIRS$compare, "layer_sensitivity_artefacts.csv"), row.names = FALSE)
print(art, digits = 3, row.names = FALSE)

# --- figure ----------------------------------------------------------------------------------------------------------
LAB <- c(V0_full = "V0 full iteration 2", V1_no_porosity = "V1 without porosity", V2_no_grain_size = "V2 without grain size",
         V3_no_lithology = "V3 without lithology", V4_lithology_only_texture = "V4 lithology as the only texture",
         V5_no_satellite_supply = "V5 without satellite POC and flux", V6_no_bottom_phytoplankton = "V6 without near-bottom phytoplankton",
         V7_bottom_layers_smoothed = "V7 bottom layers smoothed (deep)", V8_iteration1_no_texture = "V8 iteration 1 (no texture)")
fb <- bud[bud$component == "total", ]
fb$lab <- factor(LAB[fb$variant], levels = rev(LAB))
fb$ref <- ifelse(fb$reference == "V0_full", "change vs full iteration 2", "change vs iteration 1 (texture gain kept?)")
pa <- ggplot(fb, aes(pct_change, lab)) + geom_vline(xintercept = 0, colour = "grey55") +
  geom_pointrange(aes(xmin = 100 * ci_lo / value_reference, xmax = 100 * ci_hi / value_reference, colour = ci_excludes_zero), size = 0.3) +
  facet_wrap(~ref, scales = "free_x") + scale_colour_manual(values = c(`TRUE` = "#C2410C", `FALSE` = "grey35"), name = "95% CI excludes 0") +
  labs(x = "withheld-region mean squared error, % change (95% block-bootstrap CI)", y = NULL, title = "a  Do the conclusions depend on the integrated layers?") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 9))
tmpl <- init(stk[[1]], NA)
zoom <- function(nm) { r <- tmpl; r[cb] <- PRED[[nm]][match(cb, allc)]; r <- crop(r, eB)
  ggplot() + geom_spatraster(data = r) + scale_fill_viridis_c(limits = c(-1.1, 0.3), oob = scales::squish, na.value = "grey85", name = "log10 TOC") +
    coord_sf(expand = FALSE) + labs(title = LAB[[nm]]) + theme_minimal(base_size = 6.5) +
    theme(axis.text = element_blank(), plot.title = element_text(size = 7, face = "bold"), legend.key.width = unit(2, "mm")) }
pb <- wrap_plots(lapply(c("V0_full", "V3_no_lithology", "V6_no_bottom_phytoplankton", "V7_bottom_layers_smoothed"), zoom), nrow = 1, guides = "collect") +
  plot_annotation(title = "b  SE Pacific: predictions without lithology, without phytoplankton, with smoothed bottom layers")
ggsave(file.path(DIRS$compare, "fig_layer_sensitivity.png"), wrap_elements(pa) / wrap_elements(pb) + plot_layout(heights = c(1, 0.8)),
       width = 230, height = 210, units = "mm", dpi = 250, bg = "white")
msg("layer sensitivity done")
