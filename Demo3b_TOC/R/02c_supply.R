# =============================================================================
# 02c_supply.R — iteration 1b: add carbon supply to the seafloor to the predictor stack
# =============================================================================
# Layers (5 arcmin) from the NN-TOC feature set (Zenodo 10.5281/zenodo.11186224, CC-BY 4.0;
# data/raw/features/FeaturesPhrampusLee_TOCSedRate_updated), extracted by HTTP range requests
# into <cache>/toc/supply (extract_supply.R):
#   TOU_Jorgenson2022.nc                          sediment total oxygen uptake (Jorgensen et al. 2022); -1 = no data
#   SS_POC_LOG_MOL_M3-1_MODIS_Aqua_MISSION_MEANx  surface POC (MODIS Aqua mission mean; values 19-7214, not log)
#   GL_RIVERMOUTH_POC_TGCYR-1_ORNL                river POC input field (ORNL)
#   GL_RIVERMOUTH_TSS_TGYR-1_ORNL                 river suspended-sediment input field (ORNL)
#   GL_TOT_SED_THICK_M_GLOBSED_Straume            total sediment thickness, m (GlobSed; Straume et al. 2019)
# Derived: poc_flux = Martin-curve attenuation of surface POC to the seafloor,
#   log10 F(z) = log10 POC_surf - 0.858 * log10(max(z, 100) / 100)   (Martin et al. 1987)
# Aggregated to 0.1 deg by averaging (flags set to NA first), log10-transformed; small ocean gaps
# filled by focal means. Provenance and possible circularity are recorded in ANALYTIC_LOG.md.
Sys.setenv(DEMO3B_ITER = "iter1b_supply")
source("R/00_config.R")
SUP  <- file.path(DIRS$toc, "supply")
base <- rast(path_base_toc_stack())
tmpl <- base[["sst_mean"]]
ocean <- !is.na(tmpl)
dep <- base[["depth"]]; if (global(dep, "mean", na.rm = TRUE)[1, 1] < 0) dep <- -dep

rd <- function(f) { r <- rast(file.path(SUP, f))[[1]]; crs(r) <- "EPSG:4326"; r }
up <- function(r, name) { r <- resample(r, tmpl, method = "average"); names(r) <- name; r }
fill <- function(r, passes = 6) {
  for (i in seq_len(passes)) {
    miss <- global(ocean & is.na(r), "sum", na.rm = TRUE)[1, 1]
    if (miss == 0) break
    r <- focal(r, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE)
  }
  mask(r, ocean, maskvalues = FALSE)
}
tou  <- rd("TOU_Jorgenson2022.nc"); tou <- ifel(tou < 0, NA, tou)
pocs <- rd("SS_POC_LOG_MOL_M3-1_MODIS_Aqua_MISSION_MEANx.5m.grd")
rpoc <- rd("GL_RIVERMOUTH_POC_TGCYR-1_ORNL.5m.grd")
rtss <- rd("GL_RIVERMOUTH_TSS_TGYR-1_ORNL.5m.grd")
sedt <- rd("GL_TOT_SED_THICK_M_GLOBSED_Straume.5m.nc")

sup <- c(
  fill(log10(clamp(up(tou, "tou"), lower = 1e-3, values = TRUE))),
  fill(log10(up(pocs, "poc_surf"))),
  fill(log10(up(rpoc, "river_poc") + 1e-6)),
  fill(log10(up(rtss, "river_tss") + 1e-5)),
  fill(log10(clamp(up(sedt, "sed_thick"), lower = 1, values = TRUE)))
)
names(sup) <- c("tou", "poc_surf", "river_poc", "river_tss", "sed_thick")
poc_flux <- sup[["poc_surf"]] - 0.858 * log10(max(dep, 100) / 100); names(poc_flux) <- "poc_flux"
sup <- c(sup, poc_flux)             # all six layers in the stack; the iteration block is SUPPLY_VARS (tou excluded, see 00_config.R)
LAYERS <- names(sup)

path_out <- file.path(CACHE, ITERS[["iter1b_supply"]]$stack)
stk <- c(base, sup)
tmp <- sub("[.]tif$", "_writing.tif", path_out)
writeRaster(stk, tmp, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
rm(stk); invisible(gc()); unlink(path_out); stopifnot(file.rename(tmp, path_out))
stk <- rast(path_out)
s <- global(stk[[LAYERS]], c("min", "mean", "max", "notNA"), na.rm = TRUE); s$layer <- rownames(s)
s$ocean_cells_missing <- sapply(LAYERS, function(v) global(ocean & is.na(stk[[v]]), "sum", na.rm = TRUE)[1, 1])
write.csv(s, file.path(DIRS$tables, "supply_ranges.csv"), row.names = FALSE)
print(s)

# descriptive only (not used for selection): rank correlation of each candidate with TOC at observation cells
obs <- read.csv(file.path(DIRS$data, "toc_cells.csv"))
X <- stk[obs$cell][, c(CANDIDATES_BASE, LAYERS)]
rho <- sapply(names(X), function(v) cor(X[[v]], obs$y, method = "spearman", use = "complete.obs"))
rt <- data.frame(variable = names(rho), block = ifelse(names(rho) %in% LAYERS, "supply", "baseline"), rho_with_log10TOC = round(rho, 3))
rt <- rt[order(-abs(rt$rho_with_log10TOC)), ]
write.csv(rt, file.path(DIRS$tables, "supply_univariate_rho.csv"), row.names = FALSE)
print(rt, row.names = FALSE)
msg("supply stack: %s (%d layers)", path_out, nlyr(stk))
