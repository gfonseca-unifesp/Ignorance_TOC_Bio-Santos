## ================================================================
## Demo 2g — Santos Basin: the where-or-when intervals, resampled by STATION
## (MSv10, Fase 1 do ROADMAP_MSv10; ponto N3 da segunda revisão pré-submissão)
## ================================================================
## santos_T4b_intervals.R bootstrapped the 3,158 rows of each arm as if they were independent. They
## are not: the same station appears in many runs and in both test surveys, so those intervals are
## too narrow. The manuscript's own rule (Box 2) asks for blocks or stations as the replicate unit.
##
## This script recomputes the same quantities with the station as the resampling unit, and keeps the
## old script untouched:
##   A  where-or-when (T4b): mean gain of each arm and the paired difference against "new stations,
##      same survey", with a cluster bootstrap over stations (B = 2000, seed 20260915);
##      sensitivity: a two-stage bootstrap (runs, then stations within the run).
##   B  more-data test (T1): one row per station and test survey, so the same station bootstrap
##      applies directly.
##   C  guide test (T3): the pairing is by run (n = 40), which is already the right unit; reported
##      here so that the three tests are read side by side.
## Writes outputs/audit/M5_*_cluster.csv; nothing is overwritten.
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_T4b_intervals_cluster.R
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/Gustavo/OneDrive/Documentos/Ignorance_MS")
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs", "audit")
SEED <- 20260915; B <- 2000
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")

## ---- the two bootstraps -------------------------------------------------------------------------
# cluster bootstrap: resample stations with replacement, keep every row of the station drawn
boot_station <- function(df, value, id = "st", B = 2000, seed = SEED) {
  set.seed(seed)
  ids <- unique(df[[id]]); rows <- split(seq_len(nrow(df)), df[[id]])
  b <- replicate(B, { s <- sample(ids, replace = TRUE); mean(df[[value]][unlist(rows[s], use.names = FALSE)]) })
  c(mean(df[[value]]), quantile(b, c(0.025, 0.975)))
}
# two-stage: resample runs, then stations within each run drawn
boot_two_stage <- function(df, value, run = "rep", id = "st", B = 2000, seed = SEED) {
  set.seed(seed)
  runs <- unique(df[[run]])
  by_run <- split(df, df[[run]])
  b <- replicate(B, {
    rs <- sample(runs, replace = TRUE)
    v <- unlist(lapply(rs, function(r) { x <- by_run[[as.character(r)]]
      ids <- unique(x[[id]]); s <- sample(ids, replace = TRUE)
      rows <- split(seq_len(nrow(x)), x[[id]])
      x[[value]][unlist(rows[s], use.names = FALSE)] }), use.names = FALSE)
    mean(v) })
  c(mean(df[[value]]), quantile(b, c(0.025, 0.975)))
}

## ---- A. where or when ----------------------------------------------------------------------------
rows <- read.csv(file.path(OUT, "T4b_where_when_rows.csv"))
rows$gain <- rows$ign0 - rows$ign
REF <- "space: new stations, same survey"

gains <- do.call(rbind, lapply(split(rows, rows$strategy), function(x) {
  s1 <- boot_station(x, "gain"); s2 <- boot_two_stage(x, "gain")
  data.frame(strategy = x$strategy[1], n_rows = nrow(x), n_stations = length(unique(x$st)),
             gain = s1[1], lo = s1[2], hi = s1[3], lo_two_stage = s2[2], hi_two_stage = s2[3]) }))

P <- merge(rows[rows$strategy != REF, ],
           rows[rows$strategy == REF, c("rep", "base_survey", "test_survey", "st", "gain")],
           by = c("rep", "base_survey", "test_survey", "st"), suffixes = c("", "_ref"))
P$diff <- P$gain - P$gain_ref
paired <- do.call(rbind, lapply(split(P, P$strategy), function(x) {
  s1 <- boot_station(x, "diff"); s2 <- boot_two_stage(x, "diff")
  data.frame(strategy = x$strategy[1], n_pairs = nrow(x),
             minus_new_stations_same_survey = s1[1], d_lo = s1[2], d_hi = s1[3],
             d_lo_two_stage = s2[2], d_hi_two_stage = s2[3],
             ratio_vs_new_stations_same_survey = mean(x$gain) / mean(x$gain_ref)) }))

A <- merge(gains, paired, by = "strategy", all.x = TRUE)
A$resample_unit <- "station (cluster); sensitivity: runs then stations"
A$test_set <- "all test stations, both surveys"
write.csv(A, file.path(OUT, "M5_T4b_intervals_cluster.csv"), row.names = FALSE)
print(A[, c("strategy", "gain", "lo", "hi", "minus_new_stations_same_survey", "d_lo", "d_hi",
            "ratio_vs_new_stations_same_survey")], digits = 3, row.names = FALSE)
print(A[, c("strategy", "lo_two_stage", "hi_two_stage", "d_lo_two_stage", "d_hi_two_stage")], digits = 3, row.names = FALSE)

## ---- B. the more-data test ------------------------------------------------------------------------
t1 <- read.csv(file.path(OUT, "T1_time_rows.csv"))
t1w <- reshape(t1[, c("fold", "st", "test_survey", "train", "ign")], idvar = c("fold", "st", "test_survey"),
               timevar = "train", direction = "wide")
names(t1w) <- make.names(names(t1w))
stopifnot(!any(duplicated(t1w[, c("st", "test_survey")])))     # one row per station and test survey
cols <- setdiff(names(t1w), c("fold", "st", "test_survey"))
both <- grep("both", cols, value = TRUE); only <- setdiff(cols, both)
Bt <- do.call(rbind, lapply(sort(unique(t1w$test_survey)), function(sv) {
  x <- t1w[t1w$test_survey == sv, ]
  o <- only[which.max(vapply(only, function(k) sum(grepl(substr(sv, 1, 4), k)), 1))]
  x$gain <- x[[o]] - x[[both]]
  x <- x[is.finite(x$gain), ]
  s <- boot_station(x, "gain")
  data.frame(test_survey = sv, n_stations = nrow(x), column_own_survey_only = o, column_both = both,
             gain_from_adding_the_other_survey = s[1], lo = s[2], hi = s[3],
             resample_unit = "station (cluster)") }))
write.csv(Bt, file.path(OUT, "M5_T1_intervals_cluster.csv"), row.names = FALSE)
print(Bt, digits = 3, row.names = FALSE)

## ---- C. the guide test, for the record --------------------------------------------------------------
t3 <- read.csv(file.path(OUT, "T3_sampling_runs.csv"))
Ct <- do.call(rbind, lapply(split(t3, t3$pct_pool), function(x) {
  w <- reshape(x[, c("rep", "strategy", "ignorance")], idvar = "rep", timevar = "strategy", direction = "wide")
  names(w) <- sub("ignorance\\.", "", names(w))
  do.call(rbind, lapply(setdiff(names(w), c("rep", "random")), function(g) {
    d <- w$random - w[[g]]; d <- d[is.finite(d)]
    set.seed(SEED); b <- replicate(B, mean(sample(d, replace = TRUE)))
    data.frame(pct_pool = x$pct_pool[1], guide = g, n_runs = length(d), mean_advantage = mean(d),
               lo = quantile(b, 0.025), hi = quantile(b, 0.975),
               detectable_advantage_80pct_power = 2.8 * sd(d) / sqrt(length(d)),
               resample_unit = "run (the unit of pairing)") })) }))
write.csv(Ct, file.path(OUT, "M5_T3_power_cluster.csv"), row.names = FALSE)
print(Ct[Ct$pct_pool != 20, ], digits = 3, row.names = FALSE)
msg("done")
