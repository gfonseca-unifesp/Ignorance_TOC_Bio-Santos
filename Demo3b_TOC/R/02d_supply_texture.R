# =============================================================================
# 02d_supply_texture.R — iteration 2b: supply stack (02c) + texture layers (02b)
# =============================================================================
Sys.setenv(DEMO3B_ITER = "iter2b_supply_texture")
source("R/00_config.R")
sup <- rast(file.path(CACHE, ITERS[["iter1b_supply"]]$stack))
tex <- rast(path_texture_stack())[[TEXTURE_VARS]]
stopifnot(compareGeom(sup, tex, stopOnError = FALSE))
stk <- c(sup, tex)
path_out <- file.path(CACHE, ITERS[["iter2b_supply_texture"]]$stack)
tmp <- sub("[.]tif$", "_writing.tif", path_out)
writeRaster(stk, tmp, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
rm(stk); invisible(gc()); unlink(path_out); stopifnot(file.rename(tmp, path_out))
msg("supply + texture stack: %s (%d layers)", path_out, nlyr(rast(path_out)))
