## ================================================================
## Demo 2 — Santos Basin: consensus Area of Applicability of the
## structured benthic RF models (Meyer & Pebesma 2021, CAST)
## ================================================================
## Run from anywhere:
##   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_aoa.R
##
## Differences from the first recipe (all found by inspecting the savepoint):
##  * Savepoint_part1/part2 are two full iMESC savepoints, not halves of one.
##    Models + stations live in part1; the environmental prediction grid is in part2.
##  * Models are caret::train (method "rf") stored in ATTRIBUTES:
##    attr(saved_data$Wat_Sed_Org_Geo_mean2, "rf")[["<y>~Wat_Sed_Org_Geo_mean2"]]$m (biological tier)
##    attr(saved_data$Geo_mean, "rf")[["<env>~Geo_mean"]]$m (environmental tier: Long, Lat, Camp.1, Depth)
##  * x/y interface: the outcome name is in the list name, not in model$terms.
##  * CV is repeated 5x5: trainDI() stops when a row is tested more than once, so
##    folds of one repetition must be given explicitly (CVtest/CVtrain).
##  * The published CV is random, and each station occurs in both campaigns at the
##    same coordinates, so the twin of a test row sits in the training folds.
##    Primary scheme: station-grouped folds (both campaigns together).
##    Sensitivity: published folds (Rep1).
##  * The biological models see PREDICTED environment on the grid (shrunk towards
##    the training mean), so their AOA alone is falsely reassuring. The structured
##    AOA propagates applicability: a biological model is applicable in a cell only
##    if it is inside its own AOA AND inside the AOA of every environmental model
##    whose variable it uses (non-zero importance weight).
##  * DI is scaled per model; averages use DI / threshold (1 = AOA edge).
##  * Stations: training rows use the CV-based training DI (aoa() on training rows gives DI = 0).
## ================================================================

.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(CAST); library(caret); library(randomForest) })
SEED <- 20260914
set.seed(SEED)

ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")  # F0.2: ver .Renviron.example
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")

BIO <- c("Meio_N", "Meio_Nema", "Meio_Cop", "Meio_Kino", "Meio_Poly", "Meio_S",
         "Macro_N", "Macro_Annel", "Macro_Ploy", "Macro_Crust", "Macro_Moll", "Macro_S")
ENV_PRED <- c("Long", "Lat", "Camp.1", "Depth")

## ---- 1. load -------------------------------------------------------
sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
W   <- sd1$Wat_Sed_Org_Geo_mean2          # 198 stations x 44 observed env
G   <- sd1$Geo_mean                       # 198 stations: Long, Lat, Camp.1, Depth
GB  <- sd1$Grid_bat                       # grid: Long, Lat, Depth, Camp.1 (2 campaigns)
rfW <- attr(W, "rf"); rfG <- attr(G, "rf")
rm(sd1); invisible(gc())
EG  <- readRDS(file.path(ROOT, "Savepoint_part2.rds"))$saved_data$Env_predictions_grid
invisible(gc())
stopifnot(identical(rownames(GB), rownames(EG)), identical(rownames(W), rownames(G)))

bio_models <- setNames(lapply(BIO, function(y) rfW[[paste0(y, "~Wat_Sed_Org_Geo_mean2")]]$m), BIO)
env_models <- setNames(lapply(rfG, `[[`, "m"), sub("~Geo_mean$", "", names(rfG)))
stopifnot(all(vapply(bio_models, inherits, TRUE, "train")),
          all(vapply(env_models, inherits, TRUE, "train")))
bio_vars <- unique(unlist(lapply(bio_models, predictors)))
stopifnot(all(bio_vars %in% names(env_models)), all(bio_vars %in% names(EG)))
msg("bio models: %d | env models: %d | grid rows: %d (%d cells x campaigns)",
    length(bio_models), length(env_models), nrow(GB), nrow(unique(GB[, c("Long", "Lat")])))

## ---- 2. CV fold sets ---------------------------------------------------
rep1 <- function(m) {                                  # positions of repetition 1 of the iMESC 5x5 CV
  k <- grep("Rep1$", names(m$control$index))           # index is named Fold1.Rep1..; indexOut Resample01..
  if (!length(k)) k <- seq_len(m$control$number)
  stopifnot(length(unique(unlist(m$control$indexOut[k]))) == nrow(m$trainingData))
  k
}
cv_sets <- function(m, scheme) {
  if (scheme == "published") {
    k <- rep1(m)
    return(list(CVtest = unname(m$control$indexOut[k]), CVtrain = unname(m$control$index[k])))
  }
  rn <- rownames(m$trainingData)                       # station-grouped: campaign twins together
  stopifnot(all(rn %in% rownames(G)))
  st <- paste(G[rn, "Long"], G[rn, "Lat"])
  f  <- CreateSpacetimeFolds(data.frame(st = st), spacevar = "st", k = 5, seed = SEED)
  list(CVtest = f$indexOut, CVtrain = f$index)
}

run_aoa <- function(m, newdata, scheme) {
  cvs <- cv_sets(m, scheme)
  tdi <- suppressMessages(trainDI(m, CVtest = cvs$CVtest, CVtrain = cvs$CVtrain, verbose = FALSE))
  a   <- suppressMessages(aoa(newdata[, tdi$variables, drop = FALSE], model = m, trainDI = tdi, verbose = FALSE))
  w   <- unlist(tdi$weight)
  list(DI = a$DI, AOA = a$AOA, thr = tdi$threshold, tdi = tdi,
       used = tdi$variables[w > 0])
}

SCHEMES <- c("station_grouped", "published")
res <- list()
for (sch in SCHEMES) {
  msg("scheme: %s — environmental tier (%d models)", sch, length(env_models))
  res[[sch]]$env <- lapply(env_models, run_aoa, newdata = GB[, ENV_PRED], scheme = sch)
  msg("scheme: %s — biological tier (%d models)", sch, length(bio_models))
  res[[sch]]$bio <- lapply(bio_models, run_aoa, newdata = EG, scheme = sch)
}

## ---- 3. grid consensus ----------------------------------------------------
consensus <- function(r) {
  envA  <- sapply(r$env, `[[`, "AOA");  envDn <- sapply(r$env, function(x) x$DI / x$thr)
  bioA  <- sapply(r$bio, `[[`, "AOA");  bioD  <- sapply(r$bio, `[[`, "DI")
  bioDn <- sapply(r$bio, function(x) x$DI / x$thr)
  # structured: bio model inside its AOA AND inside the AOA of every env model it uses
  strA  <- sapply(names(r$bio), function(v) {
    used <- r$bio[[v]]$used
    bioA[, v] * as.integer(rowSums(envA[, used, drop = FALSE] == 0) == 0)
  })
  list(envA = envA, bioA = bioA, strA = strA,
       tab = data.frame(
         inAOA_count        = rowSums(bioA),            # 0..12, biological tier alone (original recipe)
         DI_mean            = rowMeans(bioD),           # raw, as in the original recipe (not comparable across models)
         DInorm_mean        = rowMeans(bioDn),          # mean DI/threshold, biological tier
         env_inAOA_count    = rowSums(envA),            # 0..44
         env_DInorm_mean    = rowMeans(envDn),
         struct_inAOA_count = rowSums(strA)))           # 0..12, propagated through the hierarchy
}
C  <- lapply(res, consensus)
out <- data.frame(Long = GB$Long, Lat = GB$Lat, Camp = GB$Camp.1, Depth = GB$Depth,
                  C$station_grouped$tab,
                  inAOA_count_publishedCV        = C$published$tab$inAOA_count,
                  struct_inAOA_count_publishedCV = C$published$tab$struct_inAOA_count)
write.csv(out, file.path(OUT, "santos_aoa_consensus.csv"), row.names = FALSE)

## ---- 4. per-model table ----------------------------------------------------
per_model <- do.call(rbind, lapply(SCHEMES, function(sch) do.call(rbind, lapply(c("env", "bio"), function(tier) {
  do.call(rbind, lapply(names(res[[sch]][[tier]]), function(v) {
    x <- res[[sch]][[tier]][[v]]
    data.frame(scheme = sch, tier = tier, model = v, DI_threshold = x$thr,
               trainDI_median = median(x$tdi$trainDI, na.rm = TRUE),
               pct_grid_outside_AOA_camp1 = 100 * mean(x$AOA[GB$Camp.1 == 1] == 0),
               pct_grid_outside_AOA_camp2 = 100 * mean(x$AOA[GB$Camp.1 == 2] == 0),
               n_vars_weighted = length(x$used))
  }))
}))))
write.csv(per_model, file.path(OUT, "santos_aoa_per_model.csv"), row.names = FALSE)

## ---- 5. stations: training (CV DI) and held-out test (aoa) -------------------
tr_rows <- rownames(bio_models[[1]]$trainingData)
te_rows <- setdiff(rownames(W), tr_rows)
st_tab <- lapply(names(bio_models), function(v) {
  x  <- res$station_grouped$bio[[v]]
  m  <- bio_models[[v]]
  di <- setNames(rep(NA_real_, nrow(W)), rownames(W))
  di[rownames(m$trainingData)] <- x$tdi$trainDI
  a_te <- suppressMessages(aoa(W[te_rows, x$tdi$variables, drop = FALSE], model = m, trainDI = x$tdi, verbose = FALSE))
  di[te_rows] <- a_te$DI
  di / x$thr
})
st_Dn <- do.call(cbind, st_tab)
stations <- data.frame(station = rownames(W), Long = G$Long, Lat = G$Lat, Camp = G$Camp.1, Depth = G$Depth,
                       set = ifelse(rownames(W) %in% tr_rows, "train", "test"),
                       bio_inAOA_count = rowSums(st_Dn <= 1), bio_DInorm_mean = rowMeans(st_Dn))
write.csv(stations, file.path(OUT, "santos_aoa_stations.csv"), row.names = FALSE)

## ---- 6. DI x CV error (published CV, Rep1 predictions) --------------------------
di_err <- do.call(rbind, lapply(names(bio_models), function(v) {
  m  <- bio_models[[v]]
  x  <- res$published$bio[[v]]
  pr <- m$pred
  for (tn in names(m$bestTune)) pr <- pr[pr[[tn]] == m$bestTune[[tn]], ]
  k  <- rep1(m)
  pr <- pr[pr$Resample %in% c(names(m$control$index)[k], names(m$control$indexOut)[k]), ]
  stopifnot(nrow(pr) == nrow(m$trainingData))
  pr$DI <- x$tdi$trainDI[pr$rowIndex]
  pr$abs_err_std <- abs(pr$obs - pr$pred) / sd(pr$obs)
  ct <- suppressWarnings(cor.test(pr$DI, pr$abs_err_std, method = "spearman", exact = FALSE))
  data.frame(model = v, n = nrow(pr), spearman_rho = unname(ct$estimate), p = ct$p.value)
}))
write.csv(di_err, file.path(OUT, "santos_DI_vs_CVerror.csv"), row.names = FALSE)

## ---- 7. summary -------------------------------------------------------------------
deep <- GB$Depth > max(G$Depth)
S <- function(metric, value, description) data.frame(metric, value = round(value, 2), description)
summ <- rbind(
  S("grid_rows", nrow(out), "cells x 2 campaigns"),
  S("pct_rows_bio_12of12", 100 * mean(out$inAOA_count == 12), "biological tier alone (station-grouped CV)"),
  S("pct_rows_bio_12of12_publishedCV", 100 * mean(out$inAOA_count_publishedCV == 12), "biological tier alone (published random CV)"),
  S("pct_rows_env_44of44", 100 * mean(out$env_inAOA_count == 44), "environmental tier: all 44 env models inside AOA"),
  S("pct_rows_struct_12of12", 100 * mean(out$struct_inAOA_count == 12), "structured AOA: all 12 indicators applicable"),
  S("pct_rows_struct_0of12", 100 * mean(out$struct_inAOA_count == 0), "structured AOA: no indicator applicable"),
  S("pct_rows_struct_12of12_publishedCV", 100 * mean(out$struct_inAOA_count_publishedCV == 12), "structured, published random CV"),
  S("pct_grid_deeper_than_stations", 100 * mean(deep), sprintf("grid cells deeper than %g m", max(G$Depth))),
  S("pct_deep_rows_struct_0of12", 100 * mean(out$struct_inAOA_count[deep] == 0), "of those deep cells, no indicator applicable"),
  S("pct_test_stations_bio_12of12", 100 * mean(stations$bio_inAOA_count[stations$set == "test"] == 12), "held-out stations inside all 12 AOAs"),
  S("median_spearman_DI_CVerror", median(di_err$spearman_rho), "12 biological models, published CV Rep1"),
  S("n_models_positive_DI_error", sum(di_err$spearman_rho > 0), "of 12")
)
write.csv(summ, file.path(OUT, "santos_aoa_summary.csv"), row.names = FALSE)
print(summ, right = FALSE)
print(per_model[per_model$tier == "bio", c("scheme", "model", "DI_threshold", "pct_grid_outside_AOA_camp1")], row.names = FALSE)
print(di_err, row.names = FALSE)

## ---- 8. preview (the manuscript figure is built from the CSV) ----------------------
suppressPackageStartupMessages(library(ggplot2))
pv <- rbind(data.frame(out[out$Camp == 1, c("Long", "Lat")], panel = "a  Biological tier alone (0-12)", n = out$inAOA_count[out$Camp == 1]),
            data.frame(out[out$Camp == 1, c("Long", "Lat")], panel = "b  Structured AOA, propagated (0-12)", n = out$struct_inAOA_count[out$Camp == 1]))
g <- ggplot(pv, aes(Long, Lat, colour = n)) + geom_point(size = 0.05, shape = 15) +
  scale_colour_viridis_c(option = "mako", limits = c(0, 12), name = "Indicators\ninside AOA") +
  facet_wrap(~panel) + coord_quickmap() + theme_bw(base_size = 8) +
  labs(title = "Santos Basin, campaign 1: in how many of the 12 benthic models is each cell inside the AOA")
ggsave(file.path(OUT, "santos_aoa_preview.png"), g, width = 200, height = 95, units = "mm", dpi = 250)

writeLines(capture.output(sessioninfo::session_info()), file.path(OUT, "sessionInfo.txt"))
msg("done: %s", OUT)
