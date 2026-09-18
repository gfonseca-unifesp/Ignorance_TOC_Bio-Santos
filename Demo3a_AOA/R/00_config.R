# =============================================================================
# Demo 3a — Global Area of Applicability (AOA) of a marine-ectotherm SDM
# 00_config.R : parameters, paths, packages and helpers shared by all steps
# =============================================================================
# Run every script from the project root (the folder containing run_all.R).

.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
})

CFG <- list(
  seed = 20260914,

  # --- taxon (see outputs/tables/taxon_candidates.csv, step 01) --------------
  species       = "Pocillopora damicornis",
  target_group  = "Scleractinia",   # target-group background (same sampling bias)
  max_depth_m   = 100,              # drop records reported deeper than this
  year_min      = 1970,             # drop older records (georeferencing quality)
  max_coord_unc = 10000,            # m; drop records with larger coordinate uncertainty
  gbif_max      = 50000,            # GBIF records pulled (OBIS is the primary source)

  # --- grid and habitat domain ----------------------------------------------------
  res_deg   = 0.1,                  # analysis resolution (Bio-ORACLE native = 0.05)
  shelf_m   = 200,                  # depth defining the continental-shelf domain
  habitat_depth_m = 100,            # habitat domain: cells whose shallowest seabed is <= this
  tg_enddepth     = 100,            # target-group effort only from records <= this depth

  # --- sampling design ---------------------------------------------------------
  bg_ratio        = 1,              # background per presence; 1:1 for RF (Barbet-Massin et al. 2012).
                                    # No case weights: caret adds a '.weights' column that ffs/trainDI treat as a predictor
  obis_grid_prec  = 4,              # geohash precision of OBIS effort grid (~0.35 deg)
  bg_buffer_km    = 1500,           # accessible area: background within this distance of presences
  exclude_basin   = "atlantic",     # species-specific barrier: Indo-Pacific taxon cannot reach the
                                    # Atlantic/Mediterranean (Isthmus of Panama, Suez); NULL to disable

  # --- cross-validation --------------------------------------------------------
  cv_block_km = 1000,               # spatial block size (Equal Earth projection)
  cv_k        = 5,

  # --- model -------------------------------------------------------------------
  ffs_trees   = 200,
  final_trees = 500,
  cores       = max(1, parallel::detectCores() - 2)
)

# --- predictors: Bio-ORACLE v3 (baseline 2000-2019, two decadal steps averaged)
# name = short layer name; dataset/var = ERDDAP identifiers.
# "bot" layers = benthic layer at the minimum bottom depth of each cell (depthmin).
# agg = function used to aggregate 0.05 deg -> CFG$res_deg. bathymetry_min is the
# shallowest seabed of a 0.05 deg cell (negative metres), so "max" keeps the
# shallowest seabed of the coarser cell (used for the habitat domain, not as a predictor).
BO_LAYERS <- tibble::tribble(
  ~name,           ~dataset,                                ~var,               ~decadal, ~agg,
  "sst_mean",      "thetao_baseline_2000_2019_depthsurf",   "thetao_mean",      TRUE,     "mean",
  "sst_min",       "thetao_baseline_2000_2019_depthsurf",   "thetao_min",       TRUE,     "mean",
  "sst_max",       "thetao_baseline_2000_2019_depthsurf",   "thetao_max",       TRUE,     "mean",
  "sbt_mean",      "thetao_baseline_2000_2019_depthmin",    "thetao_mean",      TRUE,     "mean",
  "sbt_min",       "thetao_baseline_2000_2019_depthmin",    "thetao_min",       TRUE,     "mean",
  "sbt_max",       "thetao_baseline_2000_2019_depthmin",    "thetao_max",       TRUE,     "mean",
  "sss_mean",      "so_baseline_2000_2019_depthsurf",       "so_mean",          TRUE,     "mean",
  "sbs_mean",      "so_baseline_2000_2019_depthmin",        "so_mean",          TRUE,     "mean",
  "o2b_mean",      "o2_baseline_2000_2018_depthmin",        "o2_mean",          TRUE,     "mean",
  "bathy",         "terrain_characteristics",               "bathymetry_mean",  FALSE,    "mean",
  "bathy_shallow", "terrain_characteristics",               "bathymetry_min",   FALSE,    "max"
)
# Derived locally (max - min) to avoid two extra downloads:
DERIVED <- list(sst_range = c("sst_max", "sst_min"),
                sbt_range = c("sbt_max", "sbt_min"))

# Atlantic + Caribbean + Gulf of Mexico + Mediterranean + Black Sea basin (lon, lat).
# Edges run over land (Central American isthmus, Andes, African interior, Sinai,
# Caucasus) and over open water at the basin limits (Drake Passage 67W, Agulhas 20E).
ATLANTIC_MED_POLY <- matrix(c(
  -100, 31,  -94.5, 17.5,  -91, 16,  -87.5, 14.3,  -85.2, 12,  -84.2, 10.2,  -81.5, 8.7,
  -79.5, 9.15,  -77.8, 8,  -76.5, 5,  -78, 0,  -78, -8,  -69, -18,  -70, -30,  -71.5, -45,
  -69.5, -53.5,  -67, -56,  -67, -65,  20, -65,  20, -34.8,  26, -15,  30, 10,  32.4, 30.6,
  34.8, 30.8,  36.5, 33.5,  44, 40,  44, 48,  40, 80,  -100, 80,  -100, 31), ncol = 2, byrow = TRUE)

# --- paths ---------------------------------------------------------------------
ROOT  <- normalizePath(getwd(), winslash = "/")
# Raw downloads are large (~4 GB transient): keep them outside OneDrive.
CACHE <- Sys.getenv("DEMO3A_CACHE", "C:/Users/fonse/Demo3a_cache")
DIRS <- list(
  cache   = CACHE,
  raw     = file.path(CACHE, "raw"),
  data    = file.path(ROOT, "data"),
  rasters = file.path(ROOT, "outputs", "rasters"),
  tables  = file.path(ROOT, "outputs", "tables"),
  figs    = file.path(ROOT, "outputs", "figures"),
  models  = file.path(ROOT, "outputs", "models")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

sp_tag <- function() gsub(" ", "_", tolower(CFG$species))
tg_tag <- function() sprintf("%s_le%dm", gsub(" ", "_", tolower(CFG$target_group)), CFG$tg_enddepth)
path_pred_stack <- function() file.path(DIRS$cache, sprintf("predictors_%sdeg.tif", CFG$res_deg))

dir.create(file.path(CACHE, "terra_tmp"), showWarnings = FALSE, recursive = TRUE)
terraOptions(progress = 0, memfrac = 0.5, tempdir = file.path(CACHE, "terra_tmp"))

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
