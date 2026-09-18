## ================================================================
## Demo 2c — Santos Basin under the updated ignorance audit
## ================================================================
## Does the approach developed on the global TOC case hold for a structured, regional,
## biological model? Three tests on the 12 benthic indicators (99 stations x 2 surveys),
## refitted with ranger on the 44 observed environmental variables (the biological tier):
##   T1 more data vs more knowledge (time): does adding the other survey's data (twice the rows)
##      reduce error on a survey, and does knowledge of one survey transfer to the other?
##   T2 decision rule by depth zone (shelf <= 200 m, upper slope 200-1000 m, lower slope > 1000 m):
##      transfer (zone withheld) vs interpolation (global model with the zone's other stations)
##      vs zone-only model, plus a learning curve (half of the zone's stations).
##   T3 directed vs random addition of stations (H7): random; geographic gap filling; driver space
##      (standardised Long, Lat, Depth = the predictors of the environmental tier); environmental
##      novelty (importance-weighted 44-variable space, DI logic); model uncertainty (QRF interval).
## Error metric = station ignorance of Supplementary Fig. S2: |obs - pred| / SD(obs), mean over indicators.
## Station-grouped designs throughout (both surveys of a station stay together unless the test is temporal).
## ================================================================
.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(ranger); library(FNN); library(ggplot2); library(patchwork) })
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")  # F0.2: ver .Renviron.example
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs", "audit"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- 20260915; set.seed(SEED); NT <- 4
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
BIO <- c("Meio_N", "Meio_Nema", "Meio_Cop", "Meio_Kino", "Meio_Poly", "Meio_S",
         "Macro_N", "Macro_Annel", "Macro_Ploy", "Macro_Crust", "Macro_Moll", "Macro_S")

sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
W <- sd1$Wat_Sed_Org_Geo_mean2; G <- sd1$Geo_mean; B <- sd1[["Biod_mean_Numeric_2024-05-20_log10"]]
rm(sd1); invisible(gc())
stopifnot(identical(rownames(B), rownames(W)), identical(rownames(G), rownames(W)))
ENV <- colnames(W)
d <- data.frame(g_Long = G$Long, g_Lat = G$Lat, survey = G$Camp.1, g_Depth = G$Depth, check.names = FALSE)
d <- cbind(d, as.data.frame(W), as.data.frame(B[, BIO]))
stopifnot(all(vapply(d[, ENV], is.numeric, TRUE)), !anyNA(d[, c(ENV, BIO)]))
names(d)[match(ENV, names(d))] <- make.names(ENV, unique = TRUE); ENV <- make.names(ENV, unique = TRUE)
# station identity: the station code ends each row name (e.g. "..._Outono_A_A01" and "..._Primavera_A_A01" = station A_A01);
# coordinates of the same station differ slightly between surveys for some stations
code <- vapply(strsplit(rownames(W), "_"), function(p) paste(tail(p, 2), collapse = "_"), "")
d$st <- code   # 100 codes: 98 stations sampled in both surveys + 2 sampled once (G_G11, P_P01; 35 km apart)
msg("station identity: station code (%d stations; %d sampled in both surveys)", length(unique(d$st)), sum(table(d$st) == 2))
sdY <- apply(d[, BIO], 2, sd)
msg("%d rows, %d stations, %d env predictors, surveys: %s", nrow(d), length(unique(d$st)), length(ENV), paste(table(d$survey), collapse = "/"))

fit_models <- function(tr) lapply(BIO, function(y) ranger(x = d[tr, ENV], y = d[[y]][tr], num.trees = 300, min.node.size = 5,
                                                          quantreg = TRUE, importance = "impurity", num.threads = NT, seed = SEED))
predict_models <- function(ms, rows) sapply(ms, function(m) predict(m, d[rows, ENV], type = "response", num.threads = NT)$predictions)
ign_rows <- function(ms, rows) rowMeans(abs(sweep(as.matrix(d[rows, BIO]) - predict_models(ms, rows), 2, sdY, "/")))
boot_ci <- function(x, B = 2000) { set.seed(SEED); b <- replicate(B, mean(sample(x, replace = TRUE))); c(mean(x), quantile(b, c(0.025, 0.975))) }

## ---- T1 more data vs more knowledge (time) -------------------------------------------------------
stations <- unique(d$st); set.seed(SEED)
fold_of <- setNames(sample(rep_len(1:5, length(stations))), stations); d$fold <- fold_of[d$st]
t1 <- list()
for (k in 1:5) {
  te <- which(d$fold == k)
  for (trset in c("2019 only", "2021 only", "both surveys")) {
    tr <- which(d$fold != k & (trset == "both surveys" | d$survey == ifelse(trset == "2019 only", 1, 2)))
    ms <- fit_models(tr)
    t1[[length(t1) + 1]] <- data.frame(fold = k, train = trset, n_train_rows = length(tr), st = d$st[te], test_survey = ifelse(d$survey[te] == 1, "2019", "2021"),
                                       ign = ign_rows(ms, te))
  }
}
T1 <- do.call(rbind, t1)
write.csv(T1, file.path(OUT, "T1_time_rows.csv"), row.names = FALSE)
t1s <- do.call(rbind, lapply(split(T1, list(T1$train, T1$test_survey)), function(x) { ci <- boot_ci(x$ign)
  data.frame(train = x$train[1], test_survey = x$test_survey[1], n_train_rows = round(mean(x$n_train_rows)), mean_ignorance = ci[1], ci_lo = ci[2], ci_hi = ci[3]) }))
pair_diff <- function(a, b, s) { x <- T1[T1$train == a & T1$test_survey == s, ]; y <- T1[T1$train == b & T1$test_survey == s, ]
  z <- merge(x[, c("st", "ign")], y[, c("st", "ign")], by = "st"); ci <- boot_ci(z$ign.y - z$ign.x)
  data.frame(comparison = sprintf("%s -> %s, tested on %s", a, b, s), difference = ci[1], ci_lo = ci[2], ci_hi = ci[3]) }
t1d <- rbind(pair_diff("2019 only", "both surveys", "2019"), pair_diff("2021 only", "both surveys", "2021"),
             pair_diff("2021 only", "2019 only", "2021"), pair_diff("2019 only", "2021 only", "2019"))
write.csv(t1s, file.path(OUT, "T1_time_summary.csv"), row.names = FALSE); write.csv(t1d, file.path(OUT, "T1_time_differences.csv"), row.names = FALSE)
print(t1s, digits = 3, row.names = FALSE); print(t1d, digits = 3, row.names = FALSE)

## ---- T2 decision rule by depth zone -------------------------------------------------------------------
d$zone <- cut(d$g_Depth, c(-Inf, 200, 1000, Inf), labels = c("shelf (<=200 m)", "upper slope (200-1000 m)", "lower slope (>1000 m)"))
t2 <- list()
for (z in levels(d$zone)) {
  inz <- which(d$zone == z); outz <- which(d$zone != z)
  zst <- unique(d$st[inz]); kz <- min(5, length(zst)); set.seed(SEED)
  zf <- setNames(sample(rep_len(seq_len(kz), length(zst))), zst)
  msT <- fit_models(outz)
  for (k in seq_len(kz)) {
    te <- inz[zf[d$st[inz]] == k]; trz <- inz[zf[d$st[inz]] != k]
    trst <- unique(d$st[trz])
    half <- lapply(1:2, function(r) { set.seed(SEED + 10 * k + r); keep <- sample(trst, ceiling(length(trst) / 2)); c(outz, trz[d$st[trz] %in% keep]) })
    g50 <- (ign_rows(fit_models(half[[1]]), te) + ign_rows(fit_models(half[[2]]), te)) / 2
    t2[[length(t2) + 1]] <- data.frame(zone = z, fold = k, st = d$st[te], T = ign_rows(msT, te), G50 = g50,
                                       G = ign_rows(fit_models(c(outz, trz)), te), R = ign_rows(fit_models(trz), te))
  }
  msg("T2 zone done: %s (%d stations)", z, length(zst))
}
T2 <- do.call(rbind, t2)
write.csv(T2, file.path(OUT, "T2_zones_rows.csv"), row.names = FALSE)
t2s <- do.call(rbind, lapply(split(T2, T2$zone), function(x) {
  st_ids <- unique(x$st); set.seed(SEED)
  bs <- t(replicate(2000, { s <- sample(st_ids, replace = TRUE); y <- do.call(rbind, lapply(s, function(q) x[x$st == q, ]))
    c(TG = mean(y$T - y$G), GR = mean(y$G - y$R), G50G = mean(y$G50 - y$G)) }))
  ci <- apply(bs, 2, quantile, c(0.025, 0.975))
  dec <- if (ci[1, "GR"] > 0) "2 regionalize the model" else if (ci[1, "TG"] > 0 && ci[1, "G50G"] > 0) "1 collect data and load into the model" else
         if (ci[1, "TG"] > 0) "1 collect data (learning curve flat)" else "3 new drivers or finer resolution"
  data.frame(zone = x$zone[1], n_stations = length(st_ids), T_transfer = mean(x$T), G50 = mean(x$G50), G_interpolation = mean(x$G), R_zone_model = mean(x$R),
             T_minus_G = mean(x$T - x$G), TG_lo = ci[1, "TG"], TG_hi = ci[2, "TG"], G_minus_R = mean(x$G - x$R), GR_lo = ci[1, "GR"], GR_hi = ci[2, "GR"],
             G50_minus_G = mean(x$G50 - x$G), G50G_lo = ci[1, "G50G"], G50G_hi = ci[2, "G50G"], path = dec)
}))
write.csv(t2s, file.path(OUT, "T2_zones_decision.csv"), row.names = FALSE)
print(t2s, digits = 3, row.names = FALSE)

## ---- T3 directed vs random addition of stations (H7) ---------------------------------------------------------
XG <- cbind(d$g_Long * cos(mean(d$g_Lat) * pi / 180), d$g_Lat)
XD <- scale(cbind(d$g_Long, d$g_Lat, d$g_Depth))
ZE <- scale(as.matrix(d[, ENV]), scale = pmax(apply(d[, ENV], 2, sd), 1e-9))
select_st <- function(s, cand, samp, ms, n_add) {
  if (s == "random") return(sample(cand, n_add))
  crow <- which(d$st %in% cand); cst <- d$st[crow]; trow <- which(d$st %in% samp)
  if (s %in% c("geographic", "drivers")) {
    X <- if (s == "geographic") XG else XD
    dmin <- get.knnx(X[trow, , drop = FALSE], X[crow, , drop = FALSE], k = 1)$nn.dist[, 1]; picked <- character(0)
    for (j in seq_len(n_add)) {
      bs <- tapply(dmin, cst, mean); bs <- bs[setdiff(names(bs), picked)]; b <- names(bs)[which.max(bs)]; picked <- c(picked, b)
      dmin <- pmin(dmin, get.knnx(X[crow[cst == b], , drop = FALSE], X[crow, , drop = FALSE], k = 1)$nn.dist[, 1])
    }
    return(picked)
  }
  sc <- if (s == "novelty") {
    imp <- rowMeans(sapply(ms, function(m) { v <- pmax(m$variable.importance[ENV], 0); v / max(sum(v), 1e-12) }))
    Z <- sweep(ZE, 2, imp, "*"); get.knnx(Z[trow, , drop = FALSE], Z[crow, , drop = FALSE], k = 1)$nn.dist[, 1]
  } else {
    rowMeans(sapply(seq_along(ms), function(j) { q <- predict(ms[[j]], d[crow, ENV], type = "quantiles", quantiles = c(0.1, 0.9), num.threads = NT)$predictions
      (q[, 2] - q[, 1]) / sdY[j] }))
  }
  bs <- tapply(sc, cst, mean); names(sort(bs, decreasing = TRUE))[seq_len(n_add)]
}
STRATS <- c("random", "geographic", "drivers", "novelty", "uncertainty")
NREP <- 40; runs <- list()
for (rep in seq_len(NREP)) {
  set.seed(SEED + rep)
  test_st <- sample(stations, round(0.2 * length(stations))); pool <- setdiff(stations, test_st)
  te <- which(d$st %in% test_st); start <- sample(pool, round(0.2 * length(pool)))
  ms0 <- fit_models(which(d$st %in% start)); e0 <- mean(ign_rows(ms0, te))
  for (s in STRATS) {
    set.seed(SEED + 1000 * rep + match(s, STRATS)); samp <- start; ms <- ms0
    runs[[length(runs) + 1]] <- data.frame(rep = rep, strategy = s, pct_pool = 20, n_stations = length(samp), ignorance = e0)
    for (target in c(0.4, 0.6)) {
      n_add <- round(target * length(pool)) - length(samp)
      samp <- c(samp, select_st(s, setdiff(pool, samp), samp, ms, n_add)); ms <- fit_models(which(d$st %in% samp))
      runs[[length(runs) + 1]] <- data.frame(rep = rep, strategy = s, pct_pool = 100 * target, n_stations = length(samp), ignorance = mean(ign_rows(ms, te)))
    }
  }
  if (rep %% 5 == 0) msg("T3 repetition %d of %d", rep, NREP)
}
T3 <- do.call(rbind, runs)
write.csv(T3, file.path(OUT, "T3_sampling_runs.csv"), row.names = FALSE)
st0 <- T3[T3$pct_pool == 20, c("rep", "strategy", "ignorance")]; names(st0)[3] <- "start"
g <- merge(T3[T3$pct_pool > 20, ], st0, by = c("rep", "strategy")); g$gain <- g$start - g$ignorance
rnd <- g[g$strategy == "random", c("rep", "pct_pool", "gain")]; names(rnd)[3] <- "gain_random"
g <- merge(g, rnd, by = c("rep", "pct_pool")); g$adv <- g$gain - g$gain_random
t3s <- do.call(rbind, lapply(split(g, list(g$strategy, g$pct_pool), drop = TRUE), function(x) { ci <- boot_ci(x$adv)
  data.frame(strategy = x$strategy[1], pct_pool = x$pct_pool[1], mean_gain = mean(x$gain), mean_gain_random = mean(x$gain_random),
             advantage_vs_random = ci[1], ci_lo = ci[2], ci_hi = ci[3], reps_better = sum(x$adv > 0), reps = nrow(x)) }))
write.csv(t3s, file.path(OUT, "T3_sampling_summary.csv"), row.names = FALSE)
print(t3s, digits = 3, row.names = FALSE)

## ---- figure ------------------------------------------------------------------------------------------------------
t1s$train <- factor(t1s$train, levels = c("2019 only", "2021 only", "both surveys"))
pa <- ggplot(t1s, aes(test_survey, mean_ignorance, colour = train)) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = c("#1D6A73", "#C2410C", "#111111"), name = "training data") +
  labs(x = "survey predicted (station-grouped CV)", y = "station ignorance (mean standardised |error|)",
       title = "a  More data vs more knowledge", subtitle = "doubling the rows with the other survey vs transfer between surveys") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
lz <- do.call(rbind, lapply(seq_len(nrow(t2s)), function(i) data.frame(zone = t2s$zone[i], path = t2s$path[i],
  model = factor(c("transfer (zone withheld)", "global + zone data", "zone-only model"), levels = c("transfer (zone withheld)", "global + zone data", "zone-only model")),
  value = c(t2s$T_transfer[i], t2s$G_interpolation[i], t2s$R_zone_model[i]))))
lz$lab <- paste0(lz$zone, "\n", lz$path)
pb <- ggplot(lz, aes(value, lab, colour = model)) + geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("#C2410C", "#1D6A73", "#7F77DD"), name = NULL) +
  labs(x = "station ignorance on the same zone stations", y = NULL, title = "b  Two errors per depth zone -> decision") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
lab3 <- c(geographic = "geographic gap filling", drivers = "driver space (Long, Lat, Depth)", novelty = "environmental novelty (DI)", uncertainty = "model uncertainty (QRF)")
h <- t3s[t3s$strategy != "random", ]; h$lab <- factor(lab3[h$strategy], levels = rev(lab3)); h$effort <- factor(paste0("to ", h$pct_pool, "% of pool"))
pc <- ggplot(h, aes(advantage_vs_random, lab, colour = effort)) + geom_vline(xintercept = 0, colour = "grey55") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = c("#1D6A73", "#C2410C"), name = NULL) +
  labs(x = "extra reduction of station ignorance vs random stations (95% CI)", y = NULL,
       title = "c  Which space should guide new stations?", subtitle = sprintf("%d random test/pool splits", NREP)) +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
ggsave(file.path(OUT, "santos_audit.png"), (pa | pb) / pc + plot_layout(heights = c(1, 0.9)), width = 230, height = 200, units = "mm", dpi = 250, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
msg("santos audit done")
