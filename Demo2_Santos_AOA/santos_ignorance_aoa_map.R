## ================================================================
## Demo 2b — Santos Basin: ignorance field (S2) combined with the AOA
## ================================================================
## Run after santos_aoa.R:
##   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_ignorance_aoa_map.R
##
## Layers
##  * Station ignorance (Supplementary Fig. S2), rebuilt from the savepoint models:
##    |obs - pred| / SD(obs) per indicator, averaged over the 12 indicators.
##    Out-of-fold predictions (mean over the 5 x 5 CV) for the 160 training stations,
##    predict() for the 38 held-out test stations. Reproduces S2 (basin means 0.39 / 0.46).
##  * Biological-tier dissimilarity: mean over the 12 models of DI / AOA threshold,
##    computed on the MODELLED environmental grid (Env_predictions_grid) against the
##    observed environment of the training stations (santos_aoa.R, station-grouped CV).
##    Continuous; every grid cell is below 1 (inside all 12 AOAs).
##  * Environmental-tier extrapolation: cells where the mean DI / threshold of the 44
##    environmental models (predictors Long, Lat, Survey, Depth) exceeds 1, i.e. where the
##    environmental predictions that feed the biological models are themselves extrapolated.
##    Patches smaller than SIEVE_CELLS raster cells are removed (cartographic generalisation
##    of the thin bands between sampled isobaths); the unsieved layer is exported too.
## ================================================================

.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({
  library(caret); library(randomForest); library(terra); library(sf)
  library(ggplot2); library(tidyterra); library(rnaturalearth)
})
set.seed(20260914)

ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")  # F0.2: ver .Renviron.example
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs")
RES  <- 0.02           # raster resolution (deg); grid spacing is ~2 km
SIEVE_CELLS <- 40      # minimum patch size kept in the extrapolation mask (~40 x 4.9 km2)
BIO <- c("Meio_N", "Meio_Nema", "Meio_Cop", "Meio_Kino", "Meio_Poly", "Meio_S",
         "Macro_N", "Macro_Annel", "Macro_Ploy", "Macro_Crust", "Macro_Moll", "Macro_S")
SURVEY <- c(`1` = "2019 (Survey 1)", `2` = "2021 (Survey 2)")

## ---- 1. station ignorance (S2) ---------------------------------------------------
sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
W <- sd1$Wat_Sed_Org_Geo_mean2; G <- sd1$Geo_mean; B <- sd1[["Biod_mean_Numeric_2024-05-20_log10"]]
rfW <- attr(W, "rf"); rm(sd1); invisible(gc())
stopifnot(identical(rownames(B), rownames(W)), identical(rownames(G), rownames(W)))

O <- as.matrix(B[, BIO]); P <- O * NA
for (v in BIO) {
  m  <- rfW[[paste0(v, "~Wat_Sed_Org_Geo_mean2")]]$m
  pr <- m$pred
  for (tn in names(m$bestTune)) pr <- pr[pr[[tn]] == m$bestTune[[tn]], ]
  oof <- tapply(pr$pred, pr$rowIndex, mean)
  P[rownames(m$trainingData)[as.integer(names(oof))], v] <- oof
  te <- setdiff(rownames(W), rownames(m$trainingData))
  P[te, v] <- predict(m, W[te, predictors(m), drop = FALSE])
}
stopifnot(!anyNA(P))
std_err <- sweep(abs(O - P), 2, apply(O, 2, sd), "/")

st <- read.csv(file.path(OUT, "santos_aoa_stations.csv"))
stopifnot(identical(st$station, rownames(W)))
st$ignorance <- rowMeans(std_err)
st$survey <- SURVEY[as.character(st$Camp)]
write.csv(cbind(st[, c("station", "Long", "Lat", "Camp", "survey", "Depth", "set", "ignorance",
                       "bio_DInorm_mean", "bio_inAOA_count")], std_err),
          file.path(OUT, "santos_station_ignorance.csv"), row.names = FALSE)

## ---- 2. grid layers -> rasters ---------------------------------------------------
x <- read.csv(file.path(OUT, "santos_aoa_consensus.csv"))
tmpl <- rast(ext(min(x$Long) - RES, max(x$Long) + RES, min(x$Lat) - RES, max(x$Lat) + RES),
             resolution = RES, crs = "EPSG:4326")
to_rast <- function(d, field) {
  r <- rasterize(as.matrix(d[, c("Long", "Lat")]), tmpl, values = d[[field]], fun = mean)
  focal(r, 3, "mean", na.policy = "only", na.rm = TRUE)          # fill single empty cells
}
layers <- list()
for (cp in c(1, 2)) {
  d <- x[x$Camp == cp, ]
  bio <- to_rast(d, "DInorm_mean"); env <- to_rast(d, "env_DInorm_mean")
  env_out <- ifel(env > 1, 1, NA)
  env_out_s <- sieve(ifel(is.na(env_out), 0, 1), threshold = SIEVE_CELLS, directions = 8)
  env_out_s <- mask(ifel(env_out_s == 1, 1, NA), bio)
  names(bio) <- "bio_DInorm"; names(env) <- "env_DInorm"; names(env_out) <- "env_outside_AOA"; names(env_out_s) <- "env_outside_AOA_sieved"
  s <- c(bio, env, env_out, env_out_s)
  writeRaster(s, file.path(OUT, sprintf("santos_ignorance_aoa_layers_camp%d.tif", cp)), overwrite = TRUE,
              gdal = "COMPRESS=DEFLATE")
  layers[[as.character(cp)]] <- s
}

## ---- 3. statistics --------------------------------------------------------------------
sp <- function(a, b) { ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)); c(unname(ct$estimate), ct$p.value) }
r_all <- sp(st$ignorance, st$bio_DInorm_mean)
r_tr  <- sp(st$ignorance[st$set == "train"], st$bio_DInorm_mean[st$set == "train"])
r_te  <- sp(st$ignorance[st$set == "test"],  st$bio_DInorm_mean[st$set == "test"])
cell_pct <- function(cp, lyr) { s <- layers[[as.character(cp)]]
  100 * global(!is.na(s[[lyr]]), "sum")[[1]] / global(!is.na(s[["bio_DInorm"]]), "sum")[[1]] }
S <- function(metric, value, description) data.frame(metric, value = round(value, 3), description)
summ <- rbind(
  S("ignorance_mean_2019", mean(st$ignorance[st$Camp == 1]), "station mean, reproduces S2 (0.39)"),
  S("ignorance_mean_2021", mean(st$ignorance[st$Camp == 2]), "station mean, reproduces S2 (0.46)"),
  S("spearman_ignorance_DI_all", r_all[1], "198 stations: ignorance vs biological-tier DI / threshold"),
  S("spearman_ignorance_DI_all_p", r_all[2], ""),
  S("spearman_ignorance_DI_train", r_tr[1], "160 training stations (CV-based DI)"),
  S("spearman_ignorance_DI_test", r_te[1], "38 held-out stations (aoa() DI)"),
  S("spearman_ignorance_DI_test_p", r_te[2], ""),
  S("grid_bio_DInorm_median_2019", median(x$DInorm_mean[x$Camp == 1]), "all cells < 1: inside every biological AOA"),
  S("grid_bio_DInorm_median_2021", median(x$DInorm_mean[x$Camp == 2]), ""),
  S("grid_bio_DInorm_max", max(x$DInorm_mean), ""),
  S("pct_grid_env_outside_AOA", 100 * mean(x$env_DInorm_mean > 1), "grid rows, environmental tier mean DI / threshold > 1"),
  S("pct_grid_env_outside_AOA_deeper2400", 100 * mean(x$env_DInorm_mean[x$Depth > 2400] > 1), "cells deeper than the deepest station"),
  S("pct_raster_env_outside_sieved_2019", cell_pct(1, "env_outside_AOA_sieved"), sprintf("after removing patches < %d cells", SIEVE_CELLS)),
  S("pct_raster_env_outside_sieved_2021", cell_pct(2, "env_outside_AOA_sieved"), "")
)
write.csv(summ, file.path(OUT, "santos_ignorance_aoa_summary.csv"), row.names = FALSE)
print(summ, right = FALSE)

## ---- 4. combined map -----------------------------------------------------------------
land <- ne_countries(scale = 50, returnclass = "sf")
hatch_for <- function(poly, spacing = 0.06) {           # diagonal hatch lines clipped to polygon
  bb <- st_bbox(poly); w <- bb["xmax"] - bb["xmin"]; h <- bb["ymax"] - bb["ymin"]
  offs <- seq(-h, w, by = spacing)
  lines <- st_sfc(lapply(offs, function(o) st_linestring(rbind(c(bb["xmin"] + o, bb["ymin"]),
                                                             c(bb["xmin"] + o + h, bb["ymax"])))), crs = 4326)
  suppressWarnings(st_intersection(lines, st_union(poly)))
}
bio_df <- do.call(rbind, lapply(names(layers), function(cp) {
  d <- as.data.frame(layers[[cp]][["bio_DInorm"]], xy = TRUE, na.rm = TRUE); d$survey <- SURVEY[cp]; d }))
outp <- do.call(rbind, lapply(names(layers), function(cp) {
  p <- st_as_sf(as.polygons(layers[[cp]][["env_outside_AOA_sieved"]], dissolve = TRUE))
  if (!nrow(p)) return(NULL)
  h <- st_sf(geometry = hatch_for(p)); h$survey <- SURVEY[cp]
  p$survey <- SURVEY[cp]; list(poly = p[, "survey"], hatch = h)
}))
polys <- do.call(rbind, outp[, "poly"]); hatches <- do.call(rbind, outp[, "hatch"])
lab <- data.frame(survey = SURVEY, txt = sprintf("mean ignorance = %.2f (n = 99)",
                  tapply(st$ignorance, st$Camp, mean)))
ign_cap <- as.numeric(quantile(st$ignorance, 0.975))

g <- ggplot() +
  geom_sf(data = land, fill = "#E6DFCF", colour = "grey45", linewidth = 0.2) +
  geom_raster(data = bio_df, aes(x, y, fill = bio_DInorm)) +
  scale_fill_gradient(low = "#F4F6F8", high = "#3E5C76", limits = c(0, 1),
                      name = "Dissimilarity to sampled conditions\n(biological tier, DI / AOA threshold; 1 = edge)") +
  geom_sf(data = polys, fill = "grey30", alpha = 0.18, colour = "grey25", linewidth = 0.25) +
  geom_sf(data = hatches, colour = "grey25", linewidth = 0.15) +
  geom_point(data = st, aes(Long, Lat, colour = pmin(ignorance, ign_cap)), size = 2.3) +
  geom_point(data = st, aes(Long, Lat), shape = 21, size = 2.3, colour = "grey15", stroke = 0.25) +
  scale_colour_distiller(palette = "YlOrRd", direction = 1,
                         name = "Station ignorance\n(mean standardised |obs - pred|, out-of-sample)") +
  geom_text(data = lab, aes(x = -48.85, y = -22.72, label = txt), hjust = 0, vjust = 1, size = 2.7) +
  facet_wrap(~survey) +
  coord_sf(xlim = c(-49, -40.6), ylim = c(-28.1, -22.6), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude",
       title = "Santos Basin: where the benthic models are ignorant, and where they extrapolate",
       subtitle = sprintf(paste0("Points: out-of-sample ignorance (Suppl. Fig. S2). Background: dissimilarity of the modelled environment to the training stations.\n",
                                 "Hatched: environmental predictions extrapolated (environmental-tier DI > AOA threshold). Stations: ignorance vs dissimilarity, Spearman rho = %.2f (test stations %.2f)."),
                          r_all[1], r_te[1])) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", legend.key.width = unit(11, "mm"), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"), strip.background = element_rect(fill = "grey92"))
ggsave(file.path(OUT, "santos_ignorance_aoa_map.png"), g, width = 250, height = 140, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "santos_ignorance_aoa_map.pdf"), g, width = 250, height = 140, units = "mm", bg = "white")
writeLines(capture.output(sessioninfo::session_info()), file.path(OUT, "sessionInfo_ignorance_map.txt"))
cat("done\n")
