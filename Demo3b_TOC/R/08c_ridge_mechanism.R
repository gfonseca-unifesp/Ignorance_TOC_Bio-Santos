# =============================================================================
# 08c_ridge_mechanism.R — narrow linear increases in expected error: rough terrain
# (ridge axes, trench walls) and texture values outside the sampled texture space?
# =============================================================================
# 08b used a ~5 deg relief proxy, which blurs narrow ridge axes and trenches. Here, at 0.1 deg:
#   (1) enrichment of increased expected error (delta > 0.01) by deciles of terrain roughness (TRI)
#       and of |1-deg local relief| in deep water (> 1500 m);
#   (2) for each texture layer, the share of cells outside the range sampled by observations
#       (1st-99th percentile of training cells), in cells with vs without increased error;
#   (3) zoom maps of delta and the texture layers over a ridge (South Atlantic) and a trench (Peru-Chile).
Sys.setenv(DEMO3B_ITER = "iter2_texture")
source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(tidyterra); library(patchwork) })

dlt <- rast(file.path(DIRS$compare, "delta_expected_RMSE_iter2_minus_iter1.tif"))
stk <- rast(path_toc_stack())
dep <- stk[["depth"]]; if (mean(values(dep), na.rm = TRUE) < 0) dep <- -dep
rel1 <- focal(dep, w = 11, fun = "mean", na.rm = TRUE) - dep            # ~1 deg window; > 0 = shallower than surroundings
tex <- c("gs_d50", "gs_d16", "porosity", "litho_gs", paste0("litho_t", 1:6))
X <- c(dlt, dep, stk[["tri"]], rel1, stk[[tex]]); names(X) <- c("delta", "depth", "tri", "relief1", tex)
v <- as.data.frame(X, na.rm = FALSE); v$area <- values(cellSize(dlt, unit = "km"))[, 1]
v <- v[!is.na(v$delta) & !is.na(v$depth), ]
v$red <- v$delta > 0.01
wm <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)])

# (1) roughness and narrow relief, deep water
deep <- v[v$depth > 1500, ]
dec <- function(x) cut(x, unique(quantile(x, 0:10 / 10, na.rm = TRUE)), include.lowest = TRUE, labels = FALSE)
deep$tri_dec <- dec(deep$tri); deep$absrel_dec <- dec(abs(deep$relief1))
t1 <- rbind(
  do.call(rbind, lapply(split(deep, deep$tri_dec), function(x) data.frame(variable = "TRI (roughness)", decile = x$tri_dec[1],
    pct_area_increased_error = 100 * wm(x$red, x$area), mean_delta = wm(x$delta, x$area)))),
  do.call(rbind, lapply(split(deep, deep$absrel_dec), function(x) data.frame(variable = "|1-deg local relief|", decile = x$absrel_dec[1],
    pct_area_increased_error = 100 * wm(x$red, x$area), mean_delta = wm(x$delta, x$area)))))
write.csv(t1, file.path(DIRS$compare, "ridge_mechanism_roughness_deciles.csv"), row.names = FALSE)
print(t1, digits = 3)

# (2) texture values outside the sampled texture space
trn <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
t2 <- do.call(rbind, lapply(tex, function(k) {
  q <- quantile(trn[[k]], c(0.01, 0.99), na.rm = TRUE); out <- v[[k]] < q[1] | v[[k]] > q[2]
  data.frame(variable = k, pct_outside_sampled_range_increased_error = 100 * wm(out[v$red], v$area[v$red]),
             pct_outside_sampled_range_other = 100 * wm(out[!v$red], v$area[!v$red]),
             pct_outside_sampled_range_increased_error_deep_rough = { s <- v$red & v$depth > 1500 & v$tri > quantile(deep$tri, 0.9, na.rm = TRUE); 100 * wm(out[s], v$area[s]) })
}))
write.csv(t2, file.path(DIRS$compare, "ridge_mechanism_texture_outside_sampled.csv"), row.names = FALSE)
print(t2, digits = 3)

# (3) zooms
boxes <- list(`South Atlantic (Mid-Atlantic Ridge)` = ext(-40, 0, -50, 5), `Peru-Chile margin (trench)` = ext(-95, -65, -45, 0))
panel <- function(r, title, pal = "viridis", div = FALSE) {
  p <- ggplot() + geom_spatraster(data = r) + coord_sf(expand = FALSE) + labs(title = title, fill = NULL) +
    theme_minimal(base_size = 6) + theme(axis.text = element_blank(), legend.key.width = unit(2, "mm"))
  if (div) p + scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-0.03, 0.03), oob = scales::squish, na.value = "grey80")
  else p + scale_fill_viridis_c(na.value = "grey80")
}
rows <- lapply(names(boxes), function(b) { e <- boxes[[b]]
  wrap_plots(panel(crop(dlt, e), paste0(b, "\nchange in expected error"), div = TRUE), panel(crop(stk[["tri"]], e), "TRI (roughness)"),
             panel(crop(stk[["gs_d50"]], e), "gs_d50 (log10)"), panel(crop(stk[["litho_gs"]], e), "litho_gs"),
             panel(crop(stk[["porosity"]], e), "porosity"), panel(crop(which.max(stk[[paste0("litho_t", 1:6)]]), e), "dominant lithology type"), nrow = 1) })
ggsave(file.path(DIRS$compare, "fig_ridge_mechanism_zooms.png"), wrap_plots(rows, ncol = 1), width = 260, height = 150, units = "mm", dpi = 250, bg = "white")
msg("ridge mechanism done")
