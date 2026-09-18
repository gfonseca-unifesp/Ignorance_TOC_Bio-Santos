# =============================================================================
# 08o_prediction_artefact.R — does iteration 2 import the geometry of texture layers into the TOC predictions?
# =============================================================================
# The iteration-2 map (09b) shows a rectangular block of low TOC in the SE Pacific and a line along the East
# Pacific Rise. Where no observation constrains the model, a categorical or gridded layer can impose its own
# geometry on the prediction. Tests:
#   (1) step test: mean |difference in predicted log10 TOC| between neighbouring 0.1-deg cells in different
#       dominant lithology types vs in the same type, for iteration 1 (no texture) and iteration 2, by distance
#       to the nearest observation cell (3 million random neighbour pairs in each direction). A step ratio
#       (across / within) that rises from iteration 1 to 2 means the predictions inherited polygon edges.
#   (2) where the large changes are: share of the area with |change| > 0.1 log10 units, by distance to data.
#   (3) zoom of the SE Pacific: TOC in both iterations, change with observation cells, texture layers.
IT1 <- Sys.getenv("CMP_FROM", "iter1b_supply"); IT2 <- Sys.getenv("CMP_TO", "iter2b_supply_texture")
Sys.setenv(DEMO3B_ITER = IT2)
source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(tidyterra); library(patchwork) })

rp  <- function(it, f) rast(file.path(ROOT, "outputs", it, "rasters", f))
y1  <- rp(IT1, "toc_prediction_log10.tif"); y2 <- rp(IT2, "toc_prediction_log10.tif")
dkm <- rp(IT2, "toc_distance_nearest_observation_km.tif")
stk <- rast(path_toc_stack())
dom <- which.max(stk[[paste0("litho_t", 1:6)]])
dep <- stk[["depth"]]
DCL <- c("< 100 km from data", "100-500 km", "> 500 km")

# --- (1) step test --------------------------------------------------------------------------------------------------
M <- function(r) as.matrix(r, wide = TRUE)
Y1 <- M(y1); Y2 <- M(y2); D <- M(dom); K <- M(dkm); Z <- M(dep)
nr <- nrow(D); nc <- ncol(D); set.seed(CFG$seed)
kh <- sample.int(nr * (nc - 1), 3e6)                        # right-hand neighbour = k + nr (column-major)
kv <- sample.int(nr * nc, 3e6); kv <- kv[kv %% nr != 0]       # lower neighbour = k + 1
mk <- function(k, off) data.frame(cross = D[k] != D[k + off], d1 = abs(Y1[k] - Y1[k + off]), d2 = abs(Y2[k] - Y2[k + off]),
                                  km = pmax(K[k], K[k + off]), z = pmin(Z[k], Z[k + off]))
P <- rbind(mk(kh, nr), mk(kv, 1)); P <- P[complete.cases(P), ]
rm(Y1, Y2, K, Z); invisible(gc())
P$dist  <- cut(P$km, c(-Inf, 100, 500, Inf), labels = DCL)
P$water <- ifelse(P$z > 1500, "deep (> 1500 m)", "shelf & slope (<= 1500 m)")
step_row <- function(x, dl, wl) data.frame(distance_to_data = dl, water = wl, pairs = nrow(x), pct_pairs_across_boundary = 100 * mean(x$cross),
  step_iter1_across = mean(x$d1[x$cross]), step_iter1_within = mean(x$d1[!x$cross]), step_ratio_iter1 = mean(x$d1[x$cross]) / mean(x$d1[!x$cross]),
  step_iter2_across = mean(x$d2[x$cross]), step_iter2_within = mean(x$d2[!x$cross]), step_ratio_iter2 = mean(x$d2[x$cross]) / mean(x$d2[!x$cross]))
st <- rbind(step_row(P, "all", "all"),
            do.call(rbind, lapply(DCL, function(dl) step_row(P[P$dist == dl, ], dl, "all"))),
            do.call(rbind, lapply(DCL, function(dl) step_row(P[P$dist == dl & P$z > 1500, ], dl, "deep (> 1500 m)"))))
write.csv(st, file.path(DIRS$compare, "prediction_artefact_step_test.csv"), row.names = FALSE)
print(st, digits = 3, row.names = FALSE)

# --- (2) where the large changes are ----------------------------------------------------------------------------------
v <- data.frame(ch = values(y2 - y1)[, 1], km = values(dkm)[, 1], a = values(cellSize(y2, unit = "km"))[, 1]); v <- v[complete.cases(v), ]
v$dist <- cut(v$km, c(-Inf, 100, 500, Inf), labels = DCL); big <- abs(v$ch) > 0.1
t2 <- do.call(rbind, lapply(DCL, function(dl) { s <- v$dist == dl
  data.frame(distance_to_data = dl, pct_ocean_area = 100 * sum(v$a[s]) / sum(v$a), pct_of_large_change_area = 100 * sum(v$a[s & big]) / sum(v$a[big]),
             pct_of_class_with_large_change = 100 * sum(v$a[s & big]) / sum(v$a[s]), mean_abs_change = weighted.mean(abs(v$ch[s]), v$a[s])) }))
t2$enrichment <- t2$pct_of_large_change_area / t2$pct_ocean_area
write.csv(t2, file.path(DIRS$compare, "prediction_artefact_by_distance.csv"), row.names = FALSE)
print(t2, digits = 3, row.names = FALSE)

# --- (3) SE Pacific zoom ------------------------------------------------------------------------------------------------
e <- ext(-140, -90, -50, 0)
obs <- read.csv(file.path(DIRS$data, "toc_cells.csv"))
ob <- obs[obs$x_lon >= xmin(e) & obs$x_lon <= xmax(e) & obs$y_lat >= ymin(e) & obs$y_lat <= ymax(e), ]
msg("SE Pacific box: %d observation cells", nrow(ob))
tocp <- function(r) clamp(10^crop(r, e) - CFG$toc_offset, lower = 0.01)
pm <- function(r, title, sc, pts = FALSE) {
  p <- ggplot() + geom_spatraster(data = r) + sc + coord_sf(expand = FALSE) + labs(title = title, fill = NULL, x = NULL, y = NULL) +
    theme_minimal(base_size = 7) + theme(legend.key.width = unit(2, "mm"), axis.text = element_text(size = 5), plot.title = element_text(face = "bold", size = 7.5))
  if (pts) p <- p + geom_point(data = ob, aes(x_lon, y_lat), inherit.aes = FALSE, size = 0.6, colour = "black")
  p
}
sc_toc <- scale_fill_viridis_c(trans = "log10", limits = c(0.05, 5), oob = scales::squish, na.value = "grey85", breaks = c(0.1, 0.5, 2))
fig <- wrap_plots(
  pm(tocp(y1), "TOC (%), iteration 1 (no texture)", sc_toc),
  pm(tocp(y2), "TOC (%), iteration 2 (+ texture)", sc_toc),
  pm(crop(y2 - y1, e), sprintf("change in log10 TOC; dots = observation cells (%d)", nrow(ob)),
     scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-0.3, 0.3), oob = scales::squish, na.value = "grey85"), pts = TRUE),
  pm(crop(dom, e), "dominant lithology type", scale_fill_viridis_c(na.value = "grey85")),
  pm(crop(stk[["gs_d50"]], e), "grain size D50 (log10 mm)", scale_fill_viridis_c(na.value = "grey85")),
  pm(crop(stk[["porosity"]], e), "porosity (%)", scale_fill_viridis_c(na.value = "grey85")), ncol = 3)
ggsave(file.path(DIRS$compare, "fig_prediction_artefact_zoom.png"), fig, width = 230, height = 175, units = "mm", dpi = 250, bg = "white")
msg("prediction artefact check done")
