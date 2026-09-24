## ================================================================
## Demo 2f — Santos Basin: intervals for the where-or-when test, and the power of the guide test
## (MSv9, pre-submission review points M5 and M6)
## ================================================================
## Reads the rows already written by santos_audit_T4.R and santos_audit.R; refits nothing.
##   A  where-or-when (T4b): gain of each way of adding the same number of rows, with a bootstrap
##      interval over stations, and the paired difference against adding new stations in the same
##      survey. Reports the test set of each comparison.
##   B  the more-data test (T1) beside it, to reconcile the two: T1 adds a survey already represented
##      in the training set; T4b adds the first rows of a survey that is not.
##   C  power of the guide test (T3): the smallest advantage over random stations that the design
##      could have detected, from the spread of the runs.
## ================================================================
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/Gustavo/OneDrive/Documentos/Ignorance_MS")
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs", "audit")
SEED <- 20260915
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")

boot_ci <- function(x, B = 2000) { set.seed(SEED); b <- replicate(B, mean(sample(x, replace = TRUE))); c(mean(x), quantile(b, c(0.025, 0.975))) }

## ---- A. where or when, with intervals ---------------------------------------------------------
rows <- read.csv(file.path(OUT, "T4b_where_when_rows.csv"))
rows$gain <- rows$ign0 - rows$ign
A <- do.call(rbind, lapply(split(rows, rows$strategy), function(x) {
  ci <- boot_ci(x$gain)
  data.frame(strategy = x$strategy[1], n_rows = nrow(x), gain = ci[1], lo = ci[2], hi = ci[3]) }))
ref <- "space: new stations, same survey"
P <- merge(rows[rows$strategy != ref, ], rows[rows$strategy == ref, c("rep", "base_survey", "test_survey", "st", "gain")],
           by = c("rep", "base_survey", "test_survey", "st"), suffixes = c("", "_ref"))
P$diff <- P$gain - P$gain_ref
B <- do.call(rbind, lapply(split(P, P$strategy), function(x) {
  ci <- boot_ci(x$diff)
  data.frame(strategy = x$strategy[1], n_pairs = nrow(x), minus_new_stations_same_survey = ci[1], lo = ci[2], hi = ci[3],
             ratio_vs_new_stations_same_survey = mean(x$gain) / mean(x$gain_ref)) }))
A <- merge(A, B, by = "strategy", all.x = TRUE)
A$test_set <- "all test stations, both surveys"
write.csv(A, file.path(OUT, "M5_T4b_intervals.csv"), row.names = FALSE)
print(A, digits = 3, row.names = FALSE)

## ---- B. the more-data test, for the reconciliation --------------------------------------------
t1 <- read.csv(file.path(OUT, "T1_time_rows.csv"))
t1w <- reshape(t1[, c("fold", "st", "test_survey", "train", "ign")], idvar = c("fold", "st", "test_survey"),
               timevar = "train", direction = "wide")
names(t1w) <- make.names(names(t1w))
cols <- setdiff(names(t1w), c("fold", "st", "test_survey"))
both <- grep("both", cols, value = TRUE); only <- setdiff(cols, both)
C <- do.call(rbind, lapply(sort(unique(t1w$test_survey)), function(sv) {
  x <- t1w[t1w$test_survey == sv, ]
  o <- only[which.max(vapply(only, function(k) sum(grepl(substr(sv, 1, 4), k)), 1))]
  g <- x[[o]] - x[[both]]
  ci <- boot_ci(g[is.finite(g)])
  data.frame(test_survey = sv, n = sum(is.finite(g)), column_own_survey_only = o, column_both = both,
             gain_from_adding_the_other_survey = ci[1], lo = ci[2], hi = ci[3]) }))
write.csv(C, file.path(OUT, "M5_T1_intervals.csv"), row.names = FALSE)
print(C, digits = 3, row.names = FALSE)

## ---- C. what the guide test could have detected ------------------------------------------------
t3 <- read.csv(file.path(OUT, "T3_sampling_runs.csv"))
D <- do.call(rbind, lapply(split(t3, t3$pct_pool), function(x) {
  w <- reshape(x[, c("rep", "strategy", "ignorance")], idvar = "rep", timevar = "strategy", direction = "wide")
  names(w) <- sub("ignorance\\.", "", names(w))
  guides <- setdiff(names(w), c("rep", "random"))
  do.call(rbind, lapply(guides, function(g) {
    dif <- w$random - w[[g]]                       # positive = the guide beats random
    sdv <- sd(dif, na.rm = TRUE); n <- sum(is.finite(dif))
    ci <- boot_ci(dif[is.finite(dif)])
    data.frame(pct_pool = x$pct_pool[1], guide = g, n_runs = n, mean_advantage = ci[1], lo = ci[2], hi = ci[3],
               sd_of_paired_difference = sdv,
               detectable_advantage_80pct_power = 2.8 * sdv / sqrt(n),
               mean_ignorance_random = mean(w$random, na.rm = TRUE)) })) }))
write.csv(D, file.path(OUT, "M5_T3_power.csv"), row.names = FALSE)
print(D, digits = 3, row.names = FALSE)
msg("done")
