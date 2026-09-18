# =============================================================================
# 03_predictors.R — Bio-ORACLE v3 predictors: download, aggregate, stack, mask
# =============================================================================
# Layers come straight from the Bio-ORACLE v3 ERDDAP server (the same service
# biooracler::download_layers() wraps), requested at native 0.05 deg and
# aggregated here by the mean to CFG$res_deg. Both decadal steps of the
# 2000-2019 baseline are averaged. Output: one aligned, ocean-masked SpatRaster.

source("R/00_config.R")
options(timeout = 3600)

erddap_url <- function(dataset, var, t_index) {
  sprintf(paste0("https://erddap.bio-oracle.org/erddap/griddap/%s.nc?",
                 "%s%%5B%d%%5D%%5B0:1:last%%5D%%5B0:1:last%%5D"),
          dataset, var, t_index)
}

# A download can arrive with the right size but degenerate content (observed:
# one decadal step of bottom thetao_max served as all zeros). Reject such layers.
valid_layer <- function(r) {
  s <- global(r, c("sd", "notNA"), na.rm = TRUE)
  z <- global(r == 0, "mean", na.rm = TRUE)[[1]]
  isTRUE(s$notNA > 1e5) && isTRUE(s$sd > 1e-3) && isTRUE(z < 0.5)
}

fetch_layer <- function(dataset, var, t_index, name, agg = "mean") {
  out <- file.path(DIRS$raw, sprintf("%s_t%d_%sdeg.tif", name, t_index, CFG$res_deg))
  if (file.exists(out)) {
    if (valid_layer(rast(out))) return(out)
    msg("  cached %s is degenerate — re-downloading", basename(out)); unlink(out)
  }
  nc <- file.path(DIRS$raw, sprintf("%s_t%d.nc", name, t_index))
  for (attempt in 1:4) {
    ok <- tryCatch({
      utils::download.file(erddap_url(dataset, var, t_index), nc,
                           mode = "wb", quiet = TRUE, method = "libcurl")
      file.size(nc) > 1e6 && valid_layer(rast(nc)[[1]])
    }, error = function(e) { msg("  retry %d (%s)", attempt, conditionMessage(e)); FALSE })
    if (ok) break
    msg("  attempt %d failed or degenerate — retrying", attempt)
    Sys.sleep(20 * attempt)
  }
  if (!ok) stop(sprintf("could not obtain a valid layer for %s t%d", name, t_index))
  r <- rast(nc)
  if (nlyr(r) > 1) r <- r[[1]]
  fact <- round(CFG$res_deg / 0.05)
  if (fact > 1) r <- aggregate(r, fact = fact, fun = agg, na.rm = TRUE)
  names(r) <- name
  writeRaster(r, out, datatype = "FLT4S", overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
  unlink(nc)
  out
}

layers <- list()
for (i in seq_len(nrow(BO_LAYERS))) {
  L <- BO_LAYERS[i, ]
  steps <- if (L$decadal) 0:1 else 0
  msg("%s (%s / %s)", L$name, L$dataset, L$var)
  files <- vapply(steps, function(s) fetch_layer(L$dataset, L$var, s, L$name, L$agg), "")
  r <- if (length(files) > 1) mean(rast(files)) else rast(files)
  names(r) <- L$name
  layers[[L$name]] <- r
}

stk <- rast(layers)
for (d in names(DERIVED)) stk[[d]] <- stk[[DERIVED[[d]][1]]] - stk[[DERIVED[[d]][2]]]

# Bathymetry sign: store as positive depth (m) below sea level.
for (b in c("bathy", "bathy_shallow")) {
  if (global(stk[[b]], "mean", na.rm = TRUE)[[1]] < 0) stk[[b]] <- -stk[[b]]
}
names(stk)[names(stk) == "bathy"] <- "depth"                  # mean seabed depth of the cell
names(stk)[names(stk) == "bathy_shallow"] <- "depth_shallow"  # shallowest seabed of the cell

# Ocean mask: cells with a value in every ocean-state layer.
ocean <- !is.na(sum(stk[[setdiff(names(stk), c("depth", "depth_shallow"))]]))
stk   <- mask(stk, ocean, maskvalues = FALSE)
for (b in c("depth", "depth_shallow")) stk[[b]] <- ifel(stk[[b]] < 0, 0, stk[[b]])

# Physical consistency checks before the stack is released to later steps.
chk <- global(c(stk[["sst_max"]] - stk[["sst_mean"]], stk[["sst_mean"]] - stk[["sst_min"]],
                stk[["sbt_max"]] - stk[["sbt_mean"]], stk[["sbt_mean"]] - stk[["sbt_min"]]),
              "min", na.rm = TRUE)[[1]]
if (any(chk < -0.5)) stop("inconsistent temperature layers (max < mean or mean < min): ",
                          paste(round(chk, 2), collapse = ", "))

# write to a temporary name and rename, so no step can read a half-written stack
tmp <- sub("\\.tif$", "_writing.tif", path_pred_stack())
writeRaster(stk, tmp, overwrite = TRUE, datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
rm(stk); gc()
unlink(path_pred_stack())
stopifnot(file.rename(tmp, path_pred_stack()))
stk <- rast(path_pred_stack())

# sanity report
s <- global(stk, c("min", "mean", "max"), na.rm = TRUE)
s$layer <- rownames(s)
write.csv(s, file.path(DIRS$tables, "predictor_ranges.csv"), row.names = FALSE)
print(s)
msg("Predictor stack written: %s (%d x %d, %d layers, %d ocean cells)",
    path_pred_stack(), nrow(stk), ncol(stk), nlyr(stk),
    global(!is.na(stk[["sst_mean"]]), "sum")[[1]])
