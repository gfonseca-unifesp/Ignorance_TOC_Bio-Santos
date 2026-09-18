# =============================================================================
# 08d_lithology_edges.R — do the linear increases in expected error follow the lithology map
# (ridge-axis class, polygon boundaries) rather than the seafloor itself?
# =============================================================================
# The zooms (08c) show a thin ridge-axis lithology class (type 6, very high litho_gs) along the
# Mid-Atlantic Ridge and red patches that follow lithology polygon edges. Quantify:
#   (a) area enrichment of increased expected error (delta > 0.01) where type 6 dominates;
#   (b) enrichment within 0.2 deg of a boundary of the dominant lithology type;
#   (c) how often observation cells sample these features;
#   (d) metadata of the lithology files (class meaning).
Sys.setenv(DEMO3B_ITER = "iter2_texture")
source("R/00_config.R")

dlt <- rast(file.path(DIRS$compare, "delta_expected_RMSE_iter2_minus_iter1.tif"))
stk <- rast(path_toc_stack())
dep <- stk[["depth"]]; if (mean(values(dep), na.rm = TRUE) < 0) dep <- -dep
L <- stk[[paste0("litho_t", 1:6)]]
dom <- which.max(L)

# boundaries of the dominant lithology type (4-neighbour change), dilated by 2 cells (0.2 deg)
m <- as.matrix(dom, wide = TRUE)
nb <- function(a, dr, dc) { out <- matrix(NA, nrow(a), ncol(a))
  r <- max(1, 1 + dr):min(nrow(a), nrow(a) + dr); cc <- max(1, 1 + dc):min(ncol(a), ncol(a) + dc)
  out[r, cc] <- a[r - dr, cc - dc]; out }
S4 <- list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
edge <- matrix(FALSE, nrow(m), ncol(m))
for (s in S4) { n <- nb(m, s[1], s[2]); edge <- edge | (!is.na(n) & !is.na(m) & n != m) }
dil <- function(e) { o <- e; for (s in S4) { n <- nb(e, s[1], s[2]); o <- o | (!is.na(n) & n) }; o }
edge2 <- dil(dil(edge))
er <- rast(dom); values(er) <- as.vector(t(edge2))

v <- data.frame(delta = values(dlt)[, 1], depth = values(dep)[, 1], t6 = values(L[[6]])[, 1],
                edge = values(er)[, 1] == 1, area = values(cellSize(dlt, unit = "km"))[, 1])
v <- v[!is.na(v$delta) & !is.na(v$depth), ]
v$red <- v$delta > 0.01
obs <- read.csv(file.path(DIRS$data, "toc_cells.csv")); oxy <- cbind(obs$x_lon, obs$y_lat)
o_t6 <- extract(L[[6]], oxy)[, 1]; o_edge <- extract(er, oxy)[, 1] == 1; o_dep <- extract(dep, oxy)[, 1]

enrich <- function(sel, label, osel) { sel[is.na(sel)] <- FALSE
  data.frame(feature = label,
             pct_ocean_area = 100 * sum(v$area[sel]) / sum(v$area),
             pct_increased_error_area = 100 * sum(v$area[sel & v$red]) / sum(v$area[v$red]),
             enrichment = (sum(v$area[sel & v$red]) / sum(v$area[v$red])) / (sum(v$area[sel]) / sum(v$area)),
             pct_of_feature_with_increased_error = 100 * sum(v$area[sel & v$red]) / sum(v$area[sel]),
             mean_delta = sum(v$delta[sel] * v$area[sel]) / sum(v$area[sel]),
             pct_observation_cells = 100 * mean(osel, na.rm = TRUE)) }
tab <- rbind(
  enrich(v$t6 > 0.5, "lithology type 6 dominant (>50% of cell)", o_t6 > 0.5),
  enrich(v$edge, "within 0.2 deg of a lithology polygon boundary", o_edge),
  enrich(v$edge & v$depth > 1500, "within 0.2 deg of a boundary, deep water (>1500 m)", o_edge & o_dep > 1500),
  enrich(!v$edge & v$depth > 1500, "away from boundaries, deep water (>1500 m)", !o_edge & o_dep > 1500),
  enrich(v$edge & v$depth <= 1500, "within 0.2 deg of a boundary, shelf & slope", o_edge & o_dep <= 1500),
  enrich(!v$edge & v$depth <= 1500, "away from boundaries, shelf & slope", !o_edge & o_dep <= 1500))
write.csv(tab, file.path(DIRS$compare, "lithology_edges_check.csv"), row.names = FALSE)
print(tab, digits = 3)

# metadata of the lithology source files
cache <- DIRS$toc  # F0.2: era "C:/Users/fonse/Demo3a_cache/toc" (= file.path(CACHE, "toc"))
fl <- list.files(cache, pattern = "(litho_maps_type6_|lithology_grain_size_global_8)[.]nc$", recursive = TRUE, full.names = TRUE)
for (f in fl) {
  cat("\n---", basename(f), "---\n")
  d <- describe(f)
  cat(head(grep("long_name|description|title|comment|units|class|litho|source|reference|history", d, value = TRUE, ignore.case = TRUE), 25), sep = "\n")
}
msg("lithology edge check done")
