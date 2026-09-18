## ================================================================
## Demo 2d — Santos Basin: directing collection in time, and the structure of the environment in space and time
## ================================================================
## Follows santos_audit.R (same data, station identity and error metric: station ignorance =
## mean standardised |obs - pred| over the 12 indicators).
##   T4a  decision rule in time. For each survey: transfer (model of the other survey only) vs interpolation
##        (both surveys) vs survey-only model, plus a learning curve (other survey + half of the survey's
##        stations). Same rule as T2.
##   T4b  where vs when to add data, at equal effort. From half of the pool stations sampled in one survey,
##        add the same number of rows as (i) new stations in the same survey (space), (ii) the same stations in
##        the other survey (time), (iii) new stations in the other survey (space + time), (iv) random rows among
##        all unsampled station-survey rows. Skill = station ignorance on test stations, in both surveys.
##   T4c  structure of the environment. Within a survey: Spearman correlation between geographic (or driver-space)
##        and environmental distance among stations. In time: environmental change at the same station between
##        surveys, against the change met by moving between stations. Overlap between the stations picked by
##        geographic gap filling, driver space and environmental novelty at the cold start of T3.
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
names(d)[match(ENV, names(d))] <- make.names(ENV, unique = TRUE); ENV <- make.names(ENV, unique = TRUE)
d$st <- vapply(strsplit(rownames(W), "_"), function(p) paste(tail(p, 2), collapse = "_"), "")
d$survey <- as.integer(as.character(d$survey)); stopifnot(all(d$survey %in% 1:2))
SV <- c("2019", "2021")
sdY <- apply(d[, BIO], 2, sd)
stations <- unique(d$st)
msg("%d rows, %d stations, surveys %s", nrow(d), length(stations), paste(table(d$survey), collapse = "/"))

fit_models <- function(tr) lapply(BIO, function(y) ranger(x = d[tr, ENV], y = d[[y]][tr], num.trees = 300, min.node.size = 5,
                                                          importance = "impurity", num.threads = NT, seed = SEED))
predict_models <- function(ms, rows) sapply(ms, function(m) predict(m, d[rows, ENV], num.threads = NT)$predictions)
ign_rows <- function(ms, rows) rowMeans(abs(sweep(as.matrix(d[rows, BIO]) - predict_models(ms, rows), 2, sdY, "/")))
boot_ci <- function(x, B = 2000) { set.seed(SEED); b <- replicate(B, mean(sample(x, replace = TRUE))); c(mean(x), quantile(b, c(0.025, 0.975))) }

## ---- T4a decision rule in time ------------------------------------------------------------------------------------
set.seed(SEED); fold_of <- setNames(sample(rep_len(1:5, length(stations))), stations); d$fold <- fold_of[d$st]   # as T1
t4a <- list()
for (k in 1:5) for (s in 1:2) {
  o <- 3 - s
  te <- which(d$fold == k & d$survey == s)
  tr_o <- which(d$fold != k & d$survey == o); tr_s <- which(d$fold != k & d$survey == s)
  sst <- unique(d$st[tr_s])
  g50 <- rowMeans(sapply(1:2, function(r) { set.seed(SEED + 10 * k + s + r); keep <- sample(sst, ceiling(length(sst) / 2))
    ign_rows(fit_models(c(tr_o, tr_s[d$st[tr_s] %in% keep])), te) }))
  t4a[[length(t4a) + 1]] <- data.frame(survey = SV[s], fold = k, st = d$st[te], T = ign_rows(fit_models(tr_o), te), G50 = g50,
                                       G = ign_rows(fit_models(c(tr_o, tr_s)), te), R = ign_rows(fit_models(tr_s), te))
}
T4a <- do.call(rbind, t4a)
write.csv(T4a, file.path(OUT, "T4a_time_rows.csv"), row.names = FALSE)
t4as <- do.call(rbind, lapply(split(T4a, T4a$survey), function(x) {
  set.seed(SEED)
  bs <- t(replicate(2000, { i <- sample(nrow(x), replace = TRUE); c(TG = mean(x$T[i] - x$G[i]), GR = mean(x$G[i] - x$R[i]), G50G = mean(x$G50[i] - x$G[i])) }))
  ci <- apply(bs, 2, quantile, c(0.025, 0.975))
  dec <- if (ci[1, "GR"] > 0) "2 separate model per survey" else if (ci[1, "TG"] > 0 && ci[1, "G50G"] > 0) "1 collect in this time interval and load into the model" else
         if (ci[1, "TG"] > 0) "1 collect in this time interval (learning curve flat)" else "3 new drivers or finer resolution"
  data.frame(survey_predicted = x$survey[1], n_stations = nrow(x), T_transfer = mean(x$T), G50 = mean(x$G50), G_both = mean(x$G), R_survey_only = mean(x$R),
             T_minus_G = mean(x$T - x$G), TG_lo = ci[1, "TG"], TG_hi = ci[2, "TG"], G_minus_R = mean(x$G - x$R), GR_lo = ci[1, "GR"], GR_hi = ci[2, "GR"],
             G50_minus_G = mean(x$G50 - x$G), G50G_lo = ci[1, "G50G"], G50G_hi = ci[2, "G50G"], path = dec)
}))
write.csv(t4as, file.path(OUT, "T4a_time_decision.csv"), row.names = FALSE)
print(t4as, digits = 3, row.names = FALSE)
msg("T4a done")

## ---- T4b where vs when to add data ---------------------------------------------------------------------------------
STR <- c("space: new stations, same survey", "time: same stations, other survey", "space + time: new stations, other survey",
         "random rows (any station, either survey)")
NREP <- 40; runs <- list()
for (rep in seq_len(NREP)) {
  set.seed(SEED + rep)
  test_st <- sample(stations, round(0.2 * length(stations))); pool <- setdiff(stations, test_st)
  te <- which(d$st %in% test_st)
  for (b in 1:2) {
    o <- 3 - b
    base_st <- sample(pool, round(length(pool) / 2)); new_st <- setdiff(pool, base_st)
    base <- which(d$st %in% base_st & d$survey == b)
    add <- list(which(d$st %in% new_st & d$survey == b), which(d$st %in% base_st & d$survey == o), which(d$st %in% new_st & d$survey == o))
    n_add <- min(lengths(add))
    add <- lapply(add, function(x) if (length(x) > n_add) sample(x, n_add) else x)
    add[[4]] <- sample(setdiff(which(d$st %in% pool), base), n_add)
    e0 <- ign_rows(fit_models(base), te)
    for (j in 1:4) runs[[length(runs) + 1]] <- data.frame(rep = rep, base_survey = SV[b], strategy = STR[j], n_base = length(base), n_add = n_add,
                                                         test_survey = SV[d$survey[te]], st = d$st[te], ign0 = e0, ign = ign_rows(fit_models(c(base, add[[j]])), te))
  }
  if (rep %% 10 == 0) msg("T4b repetition %d of %d", rep, NREP)
}
T4b <- do.call(rbind, runs)
write.csv(T4b, file.path(OUT, "T4b_where_when_rows.csv"), row.names = FALSE)
T4b$gain <- T4b$ign0 - T4b$ign
T4b$test_set <- ifelse(T4b$test_survey == T4b$base_survey, "survey already sampled", "other survey")
unit <- rbind(transform(aggregate(gain ~ rep + base_survey + strategy, T4b, mean), test_set = "both surveys"),
              aggregate(gain ~ rep + base_survey + strategy + test_set, T4b, mean))
sp <- unit[unit$strategy == STR[1], c("rep", "base_survey", "test_set", "gain")]; names(sp)[4] <- "gain_space"
u <- merge(unit, sp, by = c("rep", "base_survey", "test_set")); u$adv <- u$gain - u$gain_space
t4bs <- do.call(rbind, lapply(split(u, list(u$strategy, u$test_set), drop = TRUE), function(x) { g <- boot_ci(x$gain); a <- boot_ci(x$adv)
  data.frame(strategy = x$strategy[1], test_set = x$test_set[1], units = nrow(x), mean_gain = g[1], gain_lo = g[2], gain_hi = g[3],
             advantage_vs_space = a[1], adv_lo = a[2], adv_hi = a[3], units_better_than_space = sum(x$adv > 0)) }))
t4bs$mean_rows_added <- mean(T4b$n_add); t4bs$mean_base_rows <- mean(T4b$n_base)
write.csv(t4bs, file.path(OUT, "T4b_where_when_summary.csv"), row.names = FALSE)
print(t4bs, digits = 3, row.names = FALSE)
msg("T4b done")

## ---- T4c structure of the environment in space and time -------------------------------------------------------------
Zs  <- scale(as.matrix(d[, ENV]), scale = pmax(apply(d[, ENV], 2, sd), 1e-9))
imp <- rowMeans(sapply(fit_models(seq_len(nrow(d))), function(m) { v <- pmax(m$variable.importance[ENV], 0); v / max(sum(v), 1e-12) }))
Zw  <- sweep(Zs, 2, imp, "*")
km  <- function(i) cbind(d$g_Long[i] * 111.32 * cos(mean(d$g_Lat) * pi / 180), d$g_Lat[i] * 110.57)
XD  <- scale(cbind(d$g_Long, d$g_Lat, d$g_Depth))
space <- do.call(rbind, lapply(1:2, function(s) { i <- which(d$survey == s); ed <- as.vector(dist(Zw[i, ]))
  data.frame(survey = SV[s], n_stations = length(i), rho_geographic_env = cor(as.vector(dist(km(i))), ed, method = "spearman"),
             rho_depth_env = cor(as.vector(dist(d$g_Depth[i])), ed, method = "spearman"),
             rho_drivers_env = cor(as.vector(dist(XD[i, ])), ed, method = "spearman"), median_env_distance_between_stations = median(ed)) }))
both <- names(which(table(d$st) == 2))
r1 <- vapply(both, function(q) which(d$st == q & d$survey == 1), 1L); r2 <- vapply(both, function(q) which(d$st == q & d$survey == 2), 1L)
e_time <- sqrt(rowSums((Zw[r1, ] - Zw[r2, ])^2))
nn_env <- unlist(lapply(1:2, function(s) { i <- which(d$survey == s); nn <- get.knnx(km(i), km(i), k = 2)$nn.index[, 2]; sqrt(rowSums((Zw[i, ] - Zw[i[nn], ])^2)) }))
nn_km  <- unlist(lapply(1:2, function(s) { i <- which(d$survey == s); get.knnx(km(i), km(i), k = 2)$nn.dist[, 2] }))
bins <- do.call(rbind, lapply(1:2, function(s) { i <- which(d$survey == s); gd <- as.vector(dist(km(i))); ed <- as.vector(dist(Zw[i, ]))
  br <- seq(0, ceiling(max(gd) / 25) * 25, by = 25); bb <- cut(gd, br, include.lowest = TRUE)
  data.frame(survey = SV[s], km_mid = head(br, -1) + 12.5, median_env = as.vector(tapply(ed, bb, median)), pairs = as.vector(table(bb))) }))
bins <- bins[bins$pairs >= 20, ]
equiv_km <- vapply(SV, function(s) { x <- bins[bins$survey == s, ]; k <- which(x$median_env >= median(e_time))[1]; if (is.na(k)) NA_real_ else x$km_mid[k] }, 1)
time_tab <- data.frame(stations_in_both_surveys = length(both), median_env_change_same_station = median(e_time),
                       median_env_distance_nearest_station = median(nn_env), median_km_nearest_station = median(nn_km),
                       pct_stations_change_gt_nearest_station = 100 * mean(e_time > median(nn_env)),
                       km_between_stations_with_same_env_distance_2019 = equiv_km[1], km_between_stations_with_same_env_distance_2021 = equiv_km[2])
write.csv(space, file.path(OUT, "T4c_structure_space.csv"), row.names = FALSE)
write.csv(time_tab, file.path(OUT, "T4c_structure_time.csv"), row.names = FALSE)
print(space, digits = 3, row.names = FALSE); print(time_tab, digits = 3, row.names = FALSE)

XG <- cbind(d$g_Long * cos(mean(d$g_Lat) * pi / 180), d$g_Lat)
select_st <- function(s, cand, samp, ms, n_add) {   # as santos_audit.R (T3)
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
  im <- rowMeans(sapply(ms, function(m) { v <- pmax(m$variable.importance[ENV], 0); v / max(sum(v), 1e-12) }))
  Z <- sweep(Zs, 2, im, "*"); sc <- get.knnx(Z[trow, , drop = FALSE], Z[crow, , drop = FALSE], k = 1)$nn.dist[, 1]
  bs <- tapply(sc, cst, mean); names(sort(bs, decreasing = TRUE))[seq_len(n_add)]
}
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
ov <- do.call(rbind, lapply(seq_len(NREP), function(rep) {
  set.seed(SEED + rep)
  test_st <- sample(stations, round(0.2 * length(stations))); pool <- setdiff(stations, test_st)
  start <- sample(pool, round(0.2 * length(pool))); cand <- setdiff(pool, start); n_add <- round(0.4 * length(pool)) - length(start)
  ms0 <- fit_models(which(d$st %in% start))
  sel <- lapply(c(geographic = "geographic", drivers = "drivers", novelty = "novelty"), function(s) select_st(s, cand, start, ms0, n_add))
  e <- n_add^2 / length(cand)
  data.frame(rep = rep, n_add = n_add, n_cand = length(cand), geographic_novelty = jac(sel$geographic, sel$novelty),
             drivers_novelty = jac(sel$drivers, sel$novelty), geographic_drivers = jac(sel$geographic, sel$drivers), random_expected = e / (2 * n_add - e))
}))
ovs <- data.frame(pair = c("geographic_novelty", "drivers_novelty", "geographic_drivers", "random_expected"),
                  mean_jaccard = colMeans(ov[, c("geographic_novelty", "drivers_novelty", "geographic_drivers", "random_expected")]))
write.csv(ovs, file.path(OUT, "T4c_selection_overlap.csv"), row.names = FALSE)
print(ovs, digits = 3, row.names = FALSE)

## ---- figure ----------------------------------------------------------------------------------------------------------
lt <- do.call(rbind, lapply(seq_len(nrow(t4as)), function(i) data.frame(lab = paste0("predict ", t4as$survey_predicted[i], "\npath: collect in this interval"),
  model = factor(c("transfer (other survey only)", "both surveys", "this survey only"), levels = c("transfer (other survey only)", "both surveys", "this survey only")),
  value = c(t4as$T_transfer[i], t4as$G_both[i], t4as$R_survey_only[i]))))
pa <- ggplot(lt, aes(value, lab, colour = model)) + geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("#C2410C", "#1D6A73", "#7F77DD"), name = NULL) +
  labs(x = "station ignorance (same test stations)", y = NULL, title = "a  Two errors in time -> decision") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
h <- t4bs; h$strategy <- factor(h$strategy, levels = rev(STR)); h$test_set <- factor(h$test_set, levels = c("both surveys", "survey already sampled", "other survey"))
pb <- ggplot(h, aes(mean_gain, strategy, colour = test_set)) + geom_vline(xintercept = 0, colour = "grey55") +
  geom_pointrange(aes(xmin = gain_lo, xmax = gain_hi), position = position_dodge(width = 0.55), size = 0.3) +
  scale_colour_manual(values = c("#111111", "#1D6A73", "#C2410C"), name = "ignorance measured in") +
  labs(x = "reduction of station ignorance from adding the same number of rows (95% CI)", y = NULL,
       title = "b  Where or when to add data?", subtitle = sprintf("%d splits x 2 starting surveys; ~%d rows added to ~%d", NREP, round(mean(T4b$n_add)), round(mean(T4b$n_base)))) +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
pc <- ggplot(bins, aes(km_mid, median_env, colour = survey)) + geom_line() + geom_point(size = 1) +
  geom_hline(yintercept = median(e_time), linetype = "dashed") +
  annotate("text", x = max(bins$km_mid), y = median(e_time), label = "same station, 2019 vs 2021", hjust = 1, vjust = -0.5, size = 2.5) +
  scale_colour_manual(values = c("#1D6A73", "#C2410C"), name = NULL) +
  labs(x = "distance between stations (km)", y = "environmental distance (weighted, standardised)",
       title = "c  The environment changes in space and in time",
       subtitle = sprintf("Spearman rho (2019 / 2021)\nwith distance %.2f / %.2f; with depth %.2f / %.2f",
                          space$rho_geographic_env[1], space$rho_geographic_env[2], space$rho_depth_env[1], space$rho_depth_env[2])) +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
ggsave(file.path(OUT, "santos_T4.png"), ((pa | pc) + plot_layout(widths = c(0.9, 1.1))) / pb + plot_layout(heights = c(1, 0.9)), width = 230, height = 200, units = "mm", dpi = 250, bg = "white")
msg("santos T4 done")
