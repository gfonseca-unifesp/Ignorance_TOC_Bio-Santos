# =============================================================================
# Demo 3b — Global ignorance map of seafloor total organic carbon (TOC)
# 00_config.R : parameters, paths, packages and helpers shared by all steps
# =============================================================================
# Why TOC: every record is a measurement, so prediction error is observed directly
# (no pseudo-absences). The compilation is large (NN-TOC labels, ~110,000 records)
# and strongly biased in space: North Sea/Baltic, US east coast, China shelf, upwelling
# margins dense; Southern Hemisphere, Indo-Pacific and Southern Ocean sparse.
# Run every script from the project root (the folder containing run_all.R).

.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(terra); library(sf); library(dplyr) })

CFG <- list(
  seed = 20260915,
  res_deg = 0.1,
  toc_offset = 0.01,          # response = log10(TOC% + offset)
  toc_max = 50,               # drop implausible TOC values (%)
  cv_block_km = 1000,         # spatial CV blocks (Equal Earth)
  cv_k = 5,
  n_regions = 8,              # leave-region-out transfer test
  ffs_trees = 150,
  final_trees = 500,
  cores = max(1, parallel::detectCores() - 2)
)

# --- analytic iterations ----------------------------------------------------------------
# Each iteration keeps its own outputs (outputs/<iteration>/), so the knowledge gain between
# iterations stays documented (see ANALYTIC_LOG.md). Select with the environment variable
# DEMO3B_ITER; default = latest iteration.
CANDIDATES_BASE <- c("sbt_mean", "sbt_min", "sbt_max", "sbs_mean", "o2b_mean", "sst_mean", "sss_mean",
                     "chl_mean", "phyc_bot", "sws_bot", "depth", "slope", "tri", "tpi", "dist_coast_km")
# sediment texture (NN-TOC feature set, 5 arcmin): grain size D50/D16, porosity,
# lithology-derived grain size, and fractions of six lithology types
TEXTURE_VARS <- c("gs_d50", "gs_d16", "porosity", "litho_gs", paste0("litho_t", 1:6))
# carbon supply to the seafloor (NN-TOC feature set; built by 02c_supply.R): sediment oxygen uptake,
# surface POC, Martin-curve POC flux at the seafloor, river POC and suspended-sediment input, sediment thickness
# tou (sediment oxygen uptake) is in the stack but not in the block: it is undefined over much of the Arctic
# (381 observation cells) and its provenance is unverified; kept for a sensitivity refit
SUPPLY_VARS <- c("poc_surf", "poc_flux", "river_poc", "river_tss", "sed_thick")
LITHO_VARS  <- c("litho_gs", paste0("litho_t", 1:6))   # sample-based lithology (the part of texture that carries the gain; 08q)
# block = pre-declared variables added to the iteration-1 selection (no data-driven selection; 03_model.R)
ITERS <- list(
  iter1_baseline = list(candidates = CANDIDATES_BASE, stack = "predictors_toc_0.1deg.tif", block = NULL,
                        label = "I1 baseline: ocean-state, productivity and terrain predictors"),
  iter1b_supply  = list(candidates = c(CANDIDATES_BASE, SUPPLY_VARS), stack = "predictors_toc_supply_0.1deg.tif", block = SUPPLY_VARS,
                        label = "I1b + carbon supply (oxygen uptake, POC flux, river input, sediment thickness)"),
  iter2_texture  = list(candidates = c(CANDIDATES_BASE, TEXTURE_VARS), stack = "predictors_toc_texture_0.1deg.tif", block = TEXTURE_VARS,
                        label = "I2 + sediment texture (grain size, porosity, lithology)"),
  # presentation sequence for the manuscript and report: base model = iter1b (supply), then texture on top of it
  iter2b_supply_texture = list(candidates = c(CANDIDATES_BASE, SUPPLY_VARS, TEXTURE_VARS), stack = "predictors_toc_supply_texture_0.1deg.tif",
                               block = TEXTURE_VARS, base = "iter1b_supply",
                               label = "I2b + sediment texture on top of carbon supply"),
  # clean chain (provenance check, 08q-08u; ANALYTIC_LOG). Stacks built by 02e_clean_stacks.R keep the layer names but:
  #   dist_coast_km = distance to land masses >= 25,000 km2 (not to any island: tiny islands drew straight-edged artefacts)
  #   sbt_mean, sbs_mean, o2b_mean, phyc_bot, sws_bot = smoothed (0.5-deg mean) where depth > 1500 m (seams in Bio-ORACLE)
  # iter1c re-uses the iter1b selection (no block); iter2c adds sample-based lithology only (porosity and grain size add nothing)
  iter1c_clean = list(candidates = c(CANDIDATES_BASE, SUPPLY_VARS), stack = "predictors_toc_supply_clean_0.1deg.tif",
                      block = NULL, base = "iter1b_supply",
                      label = "I1c base model with carbon supply; corrected distance to land and smoothed bottom layers"),
  iter2c_clean = list(candidates = c(CANDIDATES_BASE, SUPPLY_VARS, TEXTURE_VARS), stack = "predictors_toc_supply_texture_clean_0.1deg.tif",
                      block = LITHO_VARS, base = "iter1c_clean",
                      label = "I2c + sample-based lithology on the clean base model")
)
ITER <- Sys.getenv("DEMO3B_ITER", "iter2_texture")
stopifnot(ITER %in% names(ITERS))
CANDIDATES <- ITERS[[ITER]]$candidates

# Extra Bio-ORACLE v3 layers (the Demo 3a cache already holds temperature, salinity, O2, bathymetry)
BO_EXTRA <- tibble::tribble(
  ~name,      ~dataset,                               ~var,                          ~decadal, ~agg,
  "chl_mean", "chl_baseline_2000_2018_depthsurf",     "chl_mean",                    TRUE,     "mean",
  "phyc_bot", "phyc_baseline_2000_2020_depthmin",     "phyc_mean",                   TRUE,     "mean",
  "sws_bot",  "sws_baseline_2000_2019_depthmin",      "sws_mean",                    TRUE,     "mean",
  "slope",    "terrain_characteristics",              "slope",                       FALSE,    "mean",
  "tri",      "terrain_characteristics",              "terrain_ruggedness_index",    FALSE,    "mean",
  "tpi",      "terrain_characteristics",              "topographic_position_index",  FALSE,    "mean"
)

ROOT  <- normalizePath(getwd(), winslash = "/")
CACHE <- Sys.getenv("DEMO3A_CACHE", "C:/Users/fonse/Demo3a_cache")
OUT_ITER <- file.path(ROOT, "outputs", ITER)
DIRS <- list(
  raw       = file.path(CACHE, "raw"),
  toc       = file.path(CACHE, "toc"),
  data      = file.path(ROOT, "data"),                 # shared across iterations (toc_cells.csv)
  data_iter = file.path(OUT_ITER, "data"),             # iteration-specific training table
  rasters   = file.path(OUT_ITER, "rasters"),
  tables    = file.path(OUT_ITER, "tables"),
  figs      = file.path(OUT_ITER, "figures"),
  models    = file.path(OUT_ITER, "models"),
  logs      = file.path(ROOT, "logs", ITER),
  compare   = file.path(ROOT, "outputs", "comparison")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(file.path(CACHE, "terra_tmp"), showWarnings = FALSE, recursive = TRUE)
terraOptions(progress = 0, memfrac = 0.5, tempdir = file.path(CACHE, "terra_tmp"))

path_base_stack     <- function() file.path(CACHE, "predictors_0.1deg.tif")          # built by Demo3a_AOA/R/02_predictors.R
path_base_toc_stack <- function() file.path(CACHE, "predictors_toc_0.1deg.tif")      # + extra layers (02_predictors.R)
path_texture_stack  <- function() file.path(CACHE, "predictors_toc_texture_0.1deg.tif")  # + texture (02b_texture.R)
path_toc_stack      <- function() file.path(CACHE, ITERS[[ITER]]$stack)               # stack of the active iteration

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
