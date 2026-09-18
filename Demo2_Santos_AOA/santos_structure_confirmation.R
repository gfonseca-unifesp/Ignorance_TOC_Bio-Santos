## ================================================================
## Santos Basin — are the structured partitions confirmed in stations not used to discover them?
## ================================================================
## Same residuals as santos_error_structure.R (station-grouped 5-fold CV, both surveys; bias and magnitude).
## 100 random splits of the 100 stations into halves (both surveys of a station stay in the same half).
##   discovery half: survey (paired permutation within stations), depth zone and transect (station means);
##                   structured = eta2 >= 0.01 and BH q < 0.05
##   confirmation half: only the discoveries; confirmed = eta2 >= 0.01 and p < 0.05
## ================================================================
.libPaths(c("C:/Users/fonse/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(ranger) })
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")  # F0.2: ver .Renviron.example
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs", "audit")
SEED <- 20260915; NT <- 4; NSPLIT <- 100; NPERM <- 499
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n", sep = "")
f_res <- file.path(OUT, "S5_error_structure_residuals.csv")
if (!file.exists(f_res)) {
  BIO <- c("Meio_N", "Meio_Nema", "Meio_Cop", "Meio_Kino", "Meio_Poly", "Meio_S", "Macro_N", "Macro_Annel", "Macro_Ploy", "Macro_Crust", "Macro_Moll", "Macro_S")
  sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
  W <- sd1$Wat_Sed_Org_Geo_mean2; G <- sd1$Geo_mean; B <- sd1[["Biod_mean_Numeric_2024-05-20_log10"]]; rm(sd1); invisible(gc())
  ENV <- colnames(W)
  d <- data.frame(g_Lat = G$Lat, survey = G$Camp.1, g_Depth = G$Depth, check.names = FALSE)
  d <- cbind(d, as.data.frame(W), as.data.frame(B[, BIO]))
  names(d)[match(ENV, names(d))] <- make.names(ENV, unique = TRUE); ENV <- make.names(ENV, unique = TRUE)
  d$st <- vapply(strsplit(rownames(W), "_"), function(p) paste(tail(p, 2), collapse = "_"), "")
  d$survey <- as.integer(as.character(d$survey)); sdY <- apply(d[, BIO], 2, sd); stations <- unique(d$st)
  set.seed(SEED); fold_of <- setNames(sample(rep_len(1:5, length(stations))), stations); d$fold <- fold_of[d$st]   # as santos_error_structure.R
  E <- matrix(NA_real_, nrow(d), length(BIO))
  for (k in 1:5) { tr <- d$fold != k
    for (j in seq_along(BIO)) { m <- ranger(x = d[tr, ENV], y = d[[BIO[j]]][tr], num.trees = 300, min.node.size = 5, respect.unordered.factors = "order", num.threads = NT, seed = SEED)
      E[!tr, j] <- (predict(m, d[!tr, ENV], num.threads = NT)$predictions - d[[BIO[j]]][!tr]) / sdY[j] } }
  res <- data.frame(st = d$st, survey = c("2019", "2021")[d$survey], depth_zone = as.character(cut(d$g_Depth, c(-Inf, 200, 1000, Inf), labels = c("shelf", "upper slope", "lower slope"))),
                    transect = substr(d$st, 1, 1), bias = rowMeans(E), magnitude = rowMeans(abs(E)))
  write.csv(res, f_res, row.names = FALSE)
}
res <- read.csv(f_res); stations <- unique(res$st)
eta2 <- function(y, g) { m <- tapply(y, g, mean); n <- tapply(y, g, length); sum(n * (m - mean(y))^2) / sum((y - mean(y))^2) }
test_one <- function(r, pn, resp) {
  if (pn == "survey") {
    y <- r[[resp]]; g <- r$survey; st <- unique(r$st); o <- eta2(y, g)
    perm <- replicate(NPERM, { f <- setNames(sample(c(TRUE, FALSE), length(st), replace = TRUE), st)[r$st]; eta2(y, ifelse(f, ifelse(g == "2019", "2021", "2019"), g)) })
  } else {
    u <- aggregate(r[[resp]], list(st = r$st, level = r[[pn]]), mean); o <- eta2(u$x, u$level); perm <- replicate(NPERM, eta2(u$x, sample(u$level)))
  }
  c(eta2 = o, p = (1 + sum(perm >= o)) / (NPERM + 1))
}
out <- list(); set.seed(SEED)
for (s in seq_len(NSPLIT)) {
  disc <- sample(stations, length(stations) / 2); rows <- list()
  for (pn in c("survey", "depth_zone", "transect")) for (resp in c("bias", "magnitude")) {
    a <- test_one(res[res$st %in% disc, ], pn, resp)
    rows[[length(rows) + 1]] <- data.frame(split = s, partition = pn, response = resp, eta2_discovery = a[["eta2"]], p_discovery = a[["p"]])
  }
  R <- do.call(rbind, rows); R$q_discovery <- p.adjust(R$p_discovery, "BH"); R$discovered <- R$eta2_discovery >= 0.01 & R$q_discovery < 0.05
  R$eta2_confirmation <- NA_real_; R$p_confirmation <- NA_real_
  for (i in which(R$discovered)) { b <- test_one(res[!res$st %in% disc, ], R$partition[i], R$response[i]); R$eta2_confirmation[i] <- b[["eta2"]]; R$p_confirmation[i] <- b[["p"]] }
  R$confirmed <- R$discovered & R$eta2_confirmation >= 0.01 & R$p_confirmation < 0.05
  out[[s]] <- R
  if (s %% 20 == 0) msg("split %d of %d", s, NSPLIT)
}
A <- do.call(rbind, out)
S <- do.call(rbind, lapply(split(A, list(A$partition, A$response), drop = TRUE), function(x) data.frame(partition = x$partition[1], response = x$response[1],
  pct_splits_discovered = 100 * mean(x$discovered), pct_discoveries_confirmed = if (any(x$discovered)) 100 * mean(x$confirmed[x$discovered]) else NA,
  median_eta2_discovery = median(x$eta2_discovery), median_eta2_confirmation = if (any(x$discovered)) median(x$eta2_confirmation[x$discovered]) else NA)))
write.csv(A, file.path(OUT, "S5_error_structure_confirmation_splits.csv"), row.names = FALSE)
write.csv(S, file.path(OUT, "S5_error_structure_confirmation_summary.csv"), row.names = FALSE)
print(S, digits = 3, row.names = FALSE)
cat("\n--- illustrative split (split 1) ---\n"); print(A[A$split == 1, ], digits = 3, row.names = FALSE)
msg("Santos structure confirmation done")
