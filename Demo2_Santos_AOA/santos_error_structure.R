## ================================================================
## Santos Basin — where does the out-of-sample error have structure, and at which level of ignorance?
## ================================================================
## Protocol in Ignorance_MS/error_structure_reasoning.md. Biological tier refitted as in santos_audit.R (ranger,
## 44 observed environmental variables, 12 indicators), both surveys, station-grouped 5-fold CV (folds as T1).
## Residual per row: bias = mean standardised signed error over indicators; magnitude = station ignorance
## (mean standardised absolute error). Declared partitions (none is a predictor of the biological tier):
##   survey (time)        units = rows; permutation swaps survey labels within stations (paired)
##   depth zone           units = stations (mean of their rows); shelf <= 200 m, upper slope 200-1000 m, lower slope > 1000 m
##   transect (space)     units = stations; transects A-H (south to north) and P
## Effect size eta2; p from 2000 permutations; Benjamini-Hochberg; structured = eta2 >= 0.01 and q < 0.05.
## Out-of-sample gain: refit with the partition as a predictor on the same folds; paired station bootstrap.
## ================================================================
.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(ranger); library(ggplot2) })
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")  # F0.2: ver .Renviron.example
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs", "audit")
SEED <- 20260915; NT <- 4; NPERM <- 2000
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
BIO <- c("Meio_N", "Meio_Nema", "Meio_Cop", "Meio_Kino", "Meio_Poly", "Meio_S", "Macro_N", "Macro_Annel", "Macro_Ploy", "Macro_Crust", "Macro_Moll", "Macro_S")
sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
W <- sd1$Wat_Sed_Org_Geo_mean2; G <- sd1$Geo_mean; B <- sd1[["Biod_mean_Numeric_2024-05-20_log10"]]; rm(sd1); invisible(gc())
ENV <- colnames(W)
d <- data.frame(g_Long = G$Long, g_Lat = G$Lat, survey = G$Camp.1, g_Depth = G$Depth, check.names = FALSE)
d <- cbind(d, as.data.frame(W), as.data.frame(B[, BIO]))
names(d)[match(ENV, names(d))] <- make.names(ENV, unique = TRUE); ENV <- make.names(ENV, unique = TRUE)
d$st <- vapply(strsplit(rownames(W), "_"), function(p) paste(tail(p, 2), collapse = "_"), "")
d$survey <- as.integer(as.character(d$survey))
sdY <- apply(d[, BIO], 2, sd); stations <- unique(d$st)
set.seed(SEED); fold_of <- setNames(sample(rep_len(1:5, length(stations))), stations); d$fold <- fold_of[d$st]   # as T1
d$survey_lab <- c("2019", "2021")[d$survey]
d$depth_zone <- as.character(cut(d$g_Depth, c(-Inf, 200, 1000, Inf), labels = c("shelf", "upper slope", "lower slope")))
d$transect <- substr(d$st, 1, 1)

cv_err <- function(extra = NULL) {
  E <- matrix(NA_real_, nrow(d), length(BIO))
  for (k in 1:5) { tr <- d$fold != k; X <- d[, ENV]; if (!is.null(extra)) X$partition <- factor(d[[extra]])
    for (j in seq_along(BIO)) { m <- ranger(x = X[tr, , drop = FALSE], y = d[[BIO[j]]][tr], num.trees = 300, min.node.size = 5,
                                            respect.unordered.factors = "order", num.threads = NT, seed = SEED)
      E[!tr, j] <- (predict(m, X[!tr, , drop = FALSE], num.threads = NT)$predictions - d[[BIO[j]]][!tr]) / sdY[j] } }
  E }
E0 <- cv_err(); d$bias <- rowMeans(E0); d$magnitude <- rowMeans(abs(E0))
eta2 <- function(y, g) { m <- tapply(y, g, mean); n <- tapply(y, g, length); sum(n * (m - mean(y))^2) / sum((y - mean(y))^2) }

rows <- list(); lvl <- list()
for (pn in c("survey_lab", "depth_zone", "transect")) {
  for (resp in c("bias", "magnitude")) {
    if (pn == "survey_lab") {
      U <- d[, c("st", "survey_lab", resp)]; names(U) <- c("st", "level", "y")
      obs <- eta2(U$y, U$level); set.seed(SEED)
      perm <- replicate(NPERM, { flip <- sample(c(TRUE, FALSE), length(stations), replace = TRUE); f <- setNames(flip, stations)[U$st]
        lv <- ifelse(f, ifelse(U$level == "2019", "2021", "2019"), U$level); eta2(U$y, lv) })
    } else {
      U <- aggregate(d[[resp]], list(st = d$st, level = d[[pn]]), mean); names(U)[3] <- "y"
      obs <- eta2(U$y, U$level); set.seed(SEED); perm <- replicate(NPERM, eta2(U$y, sample(U$level)))
    }
    rows[[length(rows) + 1]] <- data.frame(partition = sub("_lab", "", pn), type = if (pn == "survey_lab") "time" else "space",
      in_model = FALSE, response = resp, n_units = nrow(U), n_levels = length(unique(U$level)), eta2 = obs, f2 = obs / (1 - obs),
      p_perm = (1 + sum(perm >= obs)) / (NPERM + 1))
    lm_ <- tapply(U$y, U$level, mean)
    lvl[[length(lvl) + 1]] <- data.frame(partition = sub("_lab", "", pn), response = resp, level = names(lm_), mean = as.numeric(lm_))
  }
}
T <- do.call(rbind, rows); T$q_BH <- p.adjust(T$p_perm, "BH")
T$structured <- T$eta2 >= 0.01 & T$q_BH < 0.05
T$level_of_ignorance <- ifelse(T$structured, paste0("L2 unnamed (", T$type, ")"), "no structure")
L <- do.call(rbind, lvl)

gain <- do.call(rbind, lapply(c("survey_lab", "depth_zone", "transect"), function(pn) {
  E1 <- cv_err(pn); g0 <- rowMeans(abs(E0)); g1 <- rowMeans(abs(E1))
  dd <- aggregate(g1 - g0, list(st = d$st), mean)$x; set.seed(SEED); bs <- replicate(1000, mean(sample(dd, replace = TRUE)))
  data.frame(partition = sub("_lab", "", pn), ignorance_without = mean(g0), ignorance_with = mean(g1), difference = mean(g1 - g0),
             ci_lo = quantile(bs, 0.025), ci_hi = quantile(bs, 0.975), pct_change = 100 * (mean(g1) - mean(g0)) / mean(g0))
}))
write.csv(T, file.path(OUT, "S5_error_structure_tests.csv"), row.names = FALSE)
write.csv(L, file.path(OUT, "S5_error_structure_level_means.csv"), row.names = FALSE)
write.csv(gain, file.path(OUT, "S5_error_structure_gain.csv"), row.names = FALSE)
print(T, digits = 3, row.names = FALSE); print(L, digits = 3, row.names = FALSE); print(gain, digits = 3, row.names = FALSE)

T$lab <- factor(T$partition, levels = rev(c("survey", "depth_zone", "transect")))
fig <- ggplot(T, aes(eta2, lab, colour = level_of_ignorance, shape = response)) + geom_vline(xintercept = 0.01, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.8, position = position_dodge(width = 0.5)) + scale_x_sqrt() +
  scale_colour_manual(values = c(`no structure` = "grey60", `L2 unnamed (time)` = "#C2410C", `L2 unnamed (space)` = "#1D6A73"), name = NULL) +
  labs(x = "eta2 (sqrt scale; dashed = 0.01)", y = NULL, title = "Santos Basin: structure of the out-of-sample error",
       subtitle = "bias = signed error; magnitude = station ignorance\nstructured = eta2 >= 0.01 and FDR-corrected permutation p < 0.05") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "santos_error_structure.png"), fig, width = 170, height = 90, units = "mm", dpi = 250, bg = "white")
msg("error structure (Santos) done")
