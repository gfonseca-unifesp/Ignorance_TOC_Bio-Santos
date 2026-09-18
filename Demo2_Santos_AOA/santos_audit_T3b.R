## ================================================================
## Demo 2c-b — Santos: why did directed sampling fail? (exploratory follow-up to T3)
## ================================================================
## T3 (santos_audit.R) found environmental novelty and model uncertainty WORSE than random stations
## at low effort. Two explanations are tested, keeping the T3 result as reported:
##   (a) batch extremes vs coverage: novelty picks the n most dissimilar stations at once, which can
##       cluster in one extreme; "env_coverage" instead fills the importance-weighted environmental
##       space sequentially (maximin), i.e. feature-space coverage;
##   (b) cold start: with ~16 starting stations the model cannot tell informative from idiosyncratic
##       stations; the same strategies are repeated from a warm start (50% of the pool -> 70%, 90%).
## Same metric, data and station codes as santos_audit.R; 40 random test/pool splits per design.
## ================================================================
.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(ranger); library(FNN); library(ggplot2) })
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")  # F0.2: ver .Renviron.example
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs", "audit")
SEED <- 20260915; NT <- 4
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
BIO <- c("Meio_N", "Meio_Nema", "Meio_Cop", "Meio_Kino", "Meio_Poly", "Meio_S",
         "Macro_N", "Macro_Annel", "Macro_Ploy", "Macro_Crust", "Macro_Moll", "Macro_S")
sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
W <- sd1$Wat_Sed_Org_Geo_mean2; G <- sd1$Geo_mean; B <- sd1[["Biod_mean_Numeric_2024-05-20_log10"]]; rm(sd1); invisible(gc())
ENV <- colnames(W)
d <- cbind(data.frame(g_Long = G$Long, g_Lat = G$Lat, survey = G$Camp.1, g_Depth = G$Depth), as.data.frame(W), as.data.frame(B[, BIO]))
names(d)[match(ENV, names(d))] <- make.names(ENV, unique = TRUE); ENV <- make.names(ENV, unique = TRUE)
d$st <- vapply(strsplit(rownames(W), "_"), function(p) paste(tail(p, 2), collapse = "_"), "")
sdY <- apply(d[, BIO], 2, sd); stations <- unique(d$st)
fit_models <- function(tr) lapply(BIO, function(y) ranger(x = d[tr, ENV], y = d[[y]][tr], num.trees = 300, min.node.size = 5,
                                                          quantreg = TRUE, importance = "impurity", num.threads = NT, seed = SEED))
ign_rows <- function(ms, rows) { P <- sapply(ms, function(m) predict(m, d[rows, ENV], num.threads = NT)$predictions)
  rowMeans(abs(sweep(as.matrix(d[rows, BIO]) - P, 2, sdY, "/"))) }
XG <- cbind(d$g_Long * cos(mean(d$g_Lat) * pi / 180), d$g_Lat)
ZE <- scale(as.matrix(d[, ENV]), scale = pmax(apply(d[, ENV], 2, sd), 1e-9))
weights_of <- function(ms) rowMeans(sapply(ms, function(m) { v <- pmax(m$variable.importance[ENV], 0); v / max(sum(v), 1e-12) }))
maximin <- function(X, crow, cst, trow, n_add) {
  dmin <- get.knnx(X[trow, , drop = FALSE], X[crow, , drop = FALSE], k = 1)$nn.dist[, 1]; picked <- character(0)
  for (j in seq_len(n_add)) { bs <- tapply(dmin, cst, mean); bs <- bs[setdiff(names(bs), picked)]; b <- names(bs)[which.max(bs)]; picked <- c(picked, b)
    dmin <- pmin(dmin, get.knnx(X[crow[cst == b], , drop = FALSE], X[crow, , drop = FALSE], k = 1)$nn.dist[, 1]) }
  picked
}
select_st <- function(s, cand, samp, ms, n_add) {
  if (s == "random") return(sample(cand, n_add))
  crow <- which(d$st %in% cand); cst <- d$st[crow]; trow <- which(d$st %in% samp)
  if (s == "geographic") return(maximin(XG, crow, cst, trow, n_add))
  Z <- sweep(ZE, 2, weights_of(ms), "*")
  if (s == "env_coverage") return(maximin(Z, crow, cst, trow, n_add))
  sc <- if (s == "novelty") get.knnx(Z[trow, , drop = FALSE], Z[crow, , drop = FALSE], k = 1)$nn.dist[, 1] else
    rowMeans(sapply(seq_along(ms), function(j) { q <- predict(ms[[j]], d[crow, ENV], type = "quantiles", quantiles = c(0.1, 0.9), num.threads = NT)$predictions
      (q[, 2] - q[, 1]) / sdY[j] }))
  bs <- tapply(sc, cst, mean); names(sort(bs, decreasing = TRUE))[seq_len(n_add)]
}
STRATS <- c("random", "geographic", "env_coverage", "novelty", "uncertainty")
DESIGNS <- list(cold = c(0.2, 0.4, 0.6), warm = c(0.5, 0.7, 0.9))
NREP <- 40; runs <- list()
for (des in names(DESIGNS)) {
  fr <- DESIGNS[[des]]
  for (rep in seq_len(NREP)) {
    set.seed(SEED + rep)
    test_st <- sample(stations, round(0.2 * length(stations))); pool <- setdiff(stations, test_st)
    te <- which(d$st %in% test_st); start <- sample(pool, round(fr[1] * length(pool)))
    ms0 <- fit_models(which(d$st %in% start)); e0 <- mean(ign_rows(ms0, te))
    for (s in STRATS) {
      set.seed(SEED + 1000 * rep + match(s, STRATS)); samp <- start; ms <- ms0
      runs[[length(runs) + 1]] <- data.frame(design = des, rep = rep, strategy = s, step = 0, pct_pool = 100 * fr[1], ignorance = e0)
      for (k in 2:3) {
        n_add <- round(fr[k] * length(pool)) - length(samp)
        samp <- c(samp, select_st(s, setdiff(pool, samp), samp, ms, n_add)); ms <- fit_models(which(d$st %in% samp))
        runs[[length(runs) + 1]] <- data.frame(design = des, rep = rep, strategy = s, step = k - 1, pct_pool = 100 * fr[k], ignorance = mean(ign_rows(ms, te)))
      }
    }
    if (rep %% 10 == 0) msg("design %s: repetition %d of %d", des, rep, NREP)
  }
}
R <- do.call(rbind, runs)
write.csv(R, file.path(OUT, "T3b_sampling_runs.csv"), row.names = FALSE)
s0 <- R[R$step == 0, c("design", "rep", "strategy", "ignorance")]; names(s0)[4] <- "start"
g <- merge(R[R$step > 0, ], s0, by = c("design", "rep", "strategy")); g$gain <- g$start - g$ignorance
rnd <- g[g$strategy == "random", c("design", "rep", "step", "gain")]; names(rnd)[4] <- "gain_random"
g <- merge(g, rnd, by = c("design", "rep", "step")); g$adv <- g$gain - g$gain_random
set.seed(SEED)
S <- do.call(rbind, lapply(split(g, list(g$design, g$strategy, g$step), drop = TRUE), function(x) {
  b <- replicate(2000, mean(sample(x$adv, replace = TRUE)))
  data.frame(design = x$design[1], strategy = x$strategy[1], step = x$step[1], pct_pool = x$pct_pool[1], mean_gain = mean(x$gain),
             mean_gain_random = mean(x$gain_random), advantage_vs_random = mean(x$adv), ci_lo = quantile(b, 0.025), ci_hi = quantile(b, 0.975),
             reps_better = sum(x$adv > 0)) }))
S <- S[order(S$design, S$step, S$strategy), ]
write.csv(S, file.path(OUT, "T3b_sampling_summary.csv"), row.names = FALSE)
print(S, digits = 3, row.names = FALSE)
msg("santos T3b done")
