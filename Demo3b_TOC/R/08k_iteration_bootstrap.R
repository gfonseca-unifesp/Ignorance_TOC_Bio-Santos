# =============================================================================
# 08k_iteration_bootstrap.R — are the differences between iterations larger than their uncertainty?
# =============================================================================
# Same withheld points (nested leave-region-out, 07_calibration.R) in two iterations, paired by location.
# Stratified block bootstrap: within each of the 8 regions, 300-km blocks are resampled with replacement
# (1000 resamples), keeping spatial autocorrelation inside blocks. For each resample the ignorance budget
# (total MSE, regional offset, bias by sediment type, unstructured error above noise) is computed for both
# iterations on the same rows; the 95% CI of the difference (second - first) is reported, overall and by region.
source("R/00_config.R")
set.seed(CFG$seed)
PAIRS <- list(c("iter1_baseline", "iter1b_supply"), c("iter1b_supply", "iter2b_supply_texture"), c("iter1_baseline", "iter2_texture"),
              c("iter1c_clean", "iter2c_clean"), c("iter1b_supply", "iter1c_clean"), c("iter2b_supply_texture", "iter2c_clean"))   # clean chain
key <- function(lon, lat) paste(round(lon, 3), round(lat, 3))
rd  <- function(it) { f <- file.path(ROOT, "outputs", it, "tables", "calibration_outer_points_toc.csv"); if (file.exists(f)) read.csv(f) else NULL }
noise <- median(read.csv(file.path(DIRS$data, "toc_cells.csv"))$y_sd, na.rm = TRUE)
tex <- rast(path_texture_stack())[[paste0("litho_t", 1:6)]]
budget <- function(err, region, lith) {
  reg <- ave(err, region); regl <- ave(err, region, lith); within <- mean((err - regl)^2); nf <- min(noise^2, within)
  c(total = mean(err^2), regional_offset = mean(reg^2), sediment_bias = mean((regl - reg)^2), unstructured_above_noise = within - nf)
}
tot <- list(); regs <- list()
for (pr in PAIRS) {
  a <- rd(pr[1]); b <- rd(pr[2]); if (is.null(a) || is.null(b)) { msg("skipping %s vs %s (missing outputs)", pr[1], pr[2]); next }
  p <- merge(data.frame(k = key(a$lon, a$lat), region = a$region, lon = a$lon, lat = a$lat, e1 = a$pred - a$obs),
             data.frame(k = key(b$lon, b$lat), e2 = b$pred - b$obs), by = "k")
  v <- as.matrix(extract(tex, cbind(p$lon, p$lat))); s <- rowSums(v)
  p$lith <- ifelse(is.na(s) | s < 0.5, "unclassified", paste0("type ", max.col(v, ties.method = "first")))
  exy <- st_coordinates(st_transform(st_as_sf(p, coords = c("lon", "lat"), crs = 4326), "EPSG:8857"))
  p$block <- paste(p$region, floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
  idx_by_block <- split(seq_len(nrow(p)), p$block)
  blocks_by_region <- split(names(idx_by_block), sub(" .*", "", names(idx_by_block)))
  est <- budget(p$e2, p$region, p$lith) - budget(p$e1, p$region, p$lith)
  base <- budget(p$e1, p$region, p$lith)
  B <- 1000
  bs <- t(replicate(B, {
    i <- unlist(lapply(blocks_by_region, function(bl) unlist(idx_by_block[sample(bl, length(bl), replace = TRUE)], use.names = FALSE)), use.names = FALSE)
    budget(p$e2[i], p$region[i], p$lith[i]) - budget(p$e1[i], p$region[i], p$lith[i])
  }))
  ci <- apply(bs, 2, quantile, c(0.025, 0.975))
  tot[[length(tot) + 1]] <- data.frame(first = pr[1], second = pr[2], component = names(est), value_first = base,
    difference = est, ci_lo = ci[1, ], ci_hi = ci[2, ], pct_change = 100 * est / base, ci_excludes_zero = ci[1, ] > 0 | ci[2, ] < 0)
  for (r in sort(unique(p$region))) {
    bl <- blocks_by_region[[as.character(r)]]; j <- which(p$region == r)
    d0 <- mean(p$e2[j]^2) - mean(p$e1[j]^2)
    br <- replicate(B, { i <- unlist(idx_by_block[sample(bl, length(bl), replace = TRUE)], use.names = FALSE); mean(p$e2[i]^2) - mean(p$e1[i]^2) })
    q <- quantile(br, c(0.025, 0.975))
    regs[[length(regs) + 1]] <- data.frame(first = pr[1], second = pr[2], region = r, n = length(j),
      rmse_first = sqrt(mean(p$e1[j]^2)), rmse_second = sqrt(mean(p$e2[j]^2)), mse_difference = d0, ci_lo = q[[1]], ci_hi = q[[2]],
      ci_excludes_zero = q[[1]] > 0 | q[[2]] < 0)
  }
  msg("%s vs %s done", pr[1], pr[2])
}
tot <- do.call(rbind, tot); regs <- do.call(rbind, regs)
write.csv(tot, file.path(DIRS$compare, "iteration_bootstrap_budget.csv"), row.names = FALSE)
write.csv(regs, file.path(DIRS$compare, "iteration_bootstrap_regions.csv"), row.names = FALSE)
print(tot, digits = 3, row.names = FALSE)
print(regs, digits = 3, row.names = FALSE)
