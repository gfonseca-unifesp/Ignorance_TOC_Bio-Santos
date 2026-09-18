# =============================================================================
# 02_predictors.R — predictor stack for seafloor TOC
# =============================================================================
# Base layers (temperature, salinity, O2, bathymetry) come from the Demo 3a cache.
# Extra Bio-ORACLE v3 layers: surface chlorophyll, near-bottom phytoplankton carbon,
# bottom current speed, slope, terrain ruggedness (TRI) and topographic position (TPI).
# Distance to coast is computed from the ocean mask.
# Skewed positive predictors are log10-transformed so that the dissimilarity index
# (scaled Euclidean distance) is not dominated by their tails; RF is unaffected.

source("R/00_config.R")
options(timeout = 3600)

erddap_url <- function(dataset, var, t_index) {
  sprintf(paste0("https://erddap.bio-oracle.org/erddap/griddap/%s.nc?",
                 "%s%%5B%d%%5D%%5B0:1:last%%5D%%5B0:1:last%%5D"), dataset, var, t_index)
}
valid_layer <- function(r) {
  s <- global(r, c("sd", "notNA"), na.rm = TRUE)
  z <- global(r == 0, "mean", na.rm = TRUE)[[1]]
  isTRUE(s$notNA > 1e5) && isTRUE(s$sd > 1e-9) && isTRUE(z < 0.9)
}
fetch_layer <- function(dataset, var, t_index, name, agg = "mean") {
  out <- file.path(DIRS$raw, sprintf("%s_t%d_%sdeg.tif", name, t_index, CFG$res_deg))
  if (file.exists(out)) { if (valid_layer(rast(out))) return(out); unlink(out) }
  nc <- file.path(DIRS$raw, sprintf("%s_t%d.nc", name, t_index))
  ok <- FALSE
  for (attempt in 1:3) {
    ok <- tryCatch({
      utils::download.file(erddap_url(dataset, var, t_index), nc, mode = "wb", quiet = TRUE, method = "libcurl")
      file.size(nc) > 1e6 && valid_layer(rast(nc)[[1]])
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok) break
    Sys.sleep(15 * attempt)
  }
  if (!ok) { unlink(nc); msg("  %s t%d unavailable or degenerate", name, t_index); return(NA_character_) }
  r <- rast(nc); if (nlyr(r) > 1) r <- r[[1]]
  fact <- round(CFG$res_deg / 0.05)
  if (fact > 1) r <- aggregate(r, fact = fact, fun = agg, na.rm = TRUE)
  names(r) <- name
  writeRaster(r, out, datatype = "FLT4S", overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  unlink(nc)
  out
}

base <- rast(path_base_stack())
ocean <- !is.na(base[["sst_mean"]])
layers <- list()
for (i in seq_len(nrow(BO_EXTRA))) {
  L <- BO_EXTRA[i, ]
  msg("%s (%s / %s)", L$name, L$dataset, L$var)
  steps <- if (L$decadal) 0:1 else 0
  files <- vapply(steps, function(s) fetch_layer(L$dataset, L$var, s, L$name, L$agg), "")
  files <- files[!is.na(files)]
  if (!length(files)) stop("no valid download for ", L$name)
  r <- if (length(files) > 1) mean(rast(files)) else rast(files)
  names(r) <- L$name
  layers[[L$name]] <- r
}
extra <- rast(layers)

# distance to coast (km), computed at 0.2 deg for speed and resampled
land  <- ifel(ocean, NA, 1)
land2 <- aggregate(land, fact = 2, fun = "max", na.rm = TRUE)
dist2 <- distance(land2) / 1000
dist  <- resample(dist2, base, method = "bilinear"); names(dist) <- "dist_coast_km"

stk <- c(base, extra, dist)
stk <- mask(stk, ocean, maskvalues = FALSE)
# log10 transforms of skewed positive predictors
for (v in c("chl_mean", "phyc_bot", "sws_bot", "slope", "tri", "dist_coast_km")) {
  lo <- max(global(stk[[v]], function(z) quantile(z[z > 0], 0.001, na.rm = TRUE))[[1]], 1e-6)
  stk[[v]] <- log10(stk[[v]] + lo)
}
tmp <- sub("\\.tif$", "_writing.tif", path_base_toc_stack())
writeRaster(stk, tmp, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
rm(stk); invisible(gc()); unlink(path_base_toc_stack()); stopifnot(file.rename(tmp, path_base_toc_stack()))
stk <- rast(path_base_toc_stack())
s <- global(stk[[CANDIDATES]], c("min", "mean", "max", "notNA"), na.rm = TRUE); s$layer <- rownames(s)
write.csv(s, file.path(DIRS$tables, "predictor_ranges_toc.csv"), row.names = FALSE)
print(s)
msg("TOC predictor stack: %s (%d layers)", path_base_toc_stack(), nlyr(stk))
