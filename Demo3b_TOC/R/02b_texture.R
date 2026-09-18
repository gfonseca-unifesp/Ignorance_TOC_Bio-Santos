# =============================================================================
# 02b_texture.R — iteration 2: add sediment texture to the predictor stack
# =============================================================================
# Layers (5 arcmin) from the NN-TOC feature set (Zenodo 10.5281/zenodo.11186224, CC-BY 4.0;
# folder data/raw/features/FeaturesPhrampusLee_TOCSedRate_updated), extracted by HTTP range
# requests into <cache>/toc/texture:
#   SF_GRAINSIZE_D50_MM_NGDC.5m.nc, SF_GRAINSIZE_D16_MM_NGDC.5m.nc  grain size (mm)
#   POROSITY_global_prediction.grd                                 porosity (%)
#   lithology_grain_size_global_8.nc                                lithology-derived grain size
#   litho_maps_type1_ ... type6_.nc                                 binary lithology masks
# Aggregated to 0.1 deg by area-weighted averaging (lithology masks become fractions);
# grain sizes are log10-transformed. Provenance and possible circularity (whether any layer
# was predicted using TOC) are recorded in ANALYTIC_LOG.md.

source("R/00_config.R")
TEX  <- file.path(DIRS$toc, "texture")
base <- rast(path_base_toc_stack())
tmpl <- base[["sst_mean"]]
ocean <- !is.na(tmpl)

rd <- function(f) { r <- rast(file.path(TEX, f))[[1]]; crs(r) <- "EPSG:4326"; r }
up <- function(r, name, log = FALSE) {
  r <- resample(r, tmpl, method = "average")
  if (log) r <- log10(clamp(r, lower = 1e-5, values = TRUE))
  names(r) <- name
  mask(r, ocean, maskvalues = FALSE)
}
tex <- c(
  up(rd("SF_GRAINSIZE_D50_MM_NGDC.5m.nc"), "gs_d50", log = TRUE),
  up(rd("SF_GRAINSIZE_D16_MM_NGDC.5m.nc"), "gs_d16", log = TRUE),
  up(rd("POROSITY_global_prediction.grd"), "porosity"),
  up(rd("lithology_grain_size_global_8.nc"), "litho_gs", log = TRUE),
  rast(lapply(1:6, function(k) up(rd(sprintf("litho_maps_type%d_.nc", k)), paste0("litho_t", k))))
)
stk <- c(base, tex)
tmp <- sub("\\.tif$", "_writing.tif", path_texture_stack())
writeRaster(stk, tmp, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
rm(stk); invisible(gc()); unlink(path_texture_stack()); stopifnot(file.rename(tmp, path_texture_stack()))
stk <- rast(path_texture_stack())
s <- global(stk[[TEXTURE_VARS]], c("min", "mean", "max", "notNA"), na.rm = TRUE); s$layer <- rownames(s)
write.csv(s, file.path(DIRS$tables, "texture_ranges.csv"), row.names = FALSE)
print(s)
lf <- global(sum(stk[[paste0("litho_t", 1:6)]]), c("min", "mean", "max"), na.rm = TRUE)
msg("sum of lithology fractions per cell: min %.2f mean %.2f max %.2f", lf$min, lf$mean, lf$max)
msg("texture stack: %s (%d layers)", path_texture_stack(), nlyr(stk))
