## ================================================================
## santos_export_data.R — F4 of ROADMAP_MSv6: the Case 2 part of the data package
## ================================================================
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_export_data.R
## The Santos models and raw biological data belong to Fonseca et al. (2026) and are distributed
## through iMESC; what is released here are the derived ignorance and applicability layers.
suppressPackageStartupMessages({ library(terra) })
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs")
REL  <- file.path(ROOT, "data_release")
dir.create(REL, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")

## ---- 1. the 198 station-campaign samples --------------------------------------------
st <- read.csv(file.path(OUT, "santos_station_ignorance.csv"))
keep <- c("station", "Long", "Lat", "Depth", "Camp", "survey", "set", "ignorance",
          "bio_DInorm_mean", "bio_inAOA_count")
out <- st[, keep]
names(out) <- c("station", "lon", "lat", "depth_m", "campaign", "survey", "set",
                "ignorance", "bio_DI_over_threshold", "bio_models_inside_AOA")
write.csv(out, file.path(REL, "case2_station_ignorance.csv"), row.names = FALSE)
msg("station samples: %d rows (%d stations x 2 surveys; %d training, %d test)", nrow(out),
    length(unique(out$station)) / 2, sum(out$set == "train"), sum(out$set == "test"))

## ---- 2. gridded layers of the v2 map, both campaigns ---------------------------------
r <- c(rast(file.path(OUT, "santos_ignorance_aoa_layers_v2_camp1.tif")),
       rast(file.path(OUT, "santos_ignorance_aoa_layers_v2_camp2.tif")))
names(r) <- c(paste0(names(r)[1:4], "_2019"), paste0(names(r)[5:8], "_2021"))
writeRaster(r, file.path(REL, "case2_grid_layers.tif"), overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"))
msg("grid raster: %d layers (%s)", nlyr(r), paste(names(r), collapse = ", "))

## ---- 3. the sampling priority of F3 ---------------------------------------------------
for (f in c("santos_priority_A_stations", "santos_priority_B_sites", "santos_priority_curve", "santos_priority_summary"))
  file.copy(file.path(OUT, paste0(f, ".csv")),
            file.path(REL, sub("^santos_priority", "case2_sampling_priority", paste0(f, ".csv"))), overwrite = TRUE)
msg("case 2 export done -> %s", REL)
