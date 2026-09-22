## ================================================================
## Demo 2e — Santos Basin: when to sample (MSv7, comentário C126 do autor)
## ================================================================
## Follows santos_audit_T4.R (same data, station identity, importance-weighted environmental space and
## station ignorance). Reads T4a_time_rows.csv; refits nothing except the importance weights (as T4c).
##   W1  transfer vs interpolation in time (T4a) split by depth zone: where does time cost most knowledge?
##   W2  environmental change at the same station between surveys, by depth zone, against the spatial difference
##       to the nearest station in the same survey.
##   W3  which environmental variables carry the change between surveys (share of the squared weighted change).
##   W4  what "survey" means: the season and cruise leg of every sample, by survey and depth zone.
## ================================================================
suppressPackageStartupMessages({ library(ranger); library(FNN) })
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/Gustavo/OneDrive/Documentos/Ignorance_MS")  # ver .Renviron.example
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs", "audit")
SEED <- 20260915; set.seed(SEED); NT <- 4
BIO <- c("Meio_N", "Meio_Nema", "Meio_Cop", "Meio_Kino", "Meio_Poly", "Meio_S",
         "Macro_N", "Macro_Annel", "Macro_Ploy", "Macro_Crust", "Macro_Moll", "Macro_S")

sd1 <- readRDS(file.path(ROOT, "Savepoint_part1.rds"))$saved_data
W <- sd1$Wat_Sed_Org_Geo_mean2; G <- sd1$Geo_mean; B <- sd1[["Biod_mean_Numeric_2024-05-20_log10"]]
rm(sd1); invisible(gc())
stopifnot(identical(rownames(B), rownames(W)), identical(rownames(G), rownames(W)))
ENV0 <- colnames(W)
d <- data.frame(g_Long = G$Long, g_Lat = G$Lat, survey = as.integer(as.character(G$Camp.1)), g_Depth = G$Depth, check.names = FALSE)
d <- cbind(d, as.data.frame(W), as.data.frame(B[, BIO]))
names(d)[match(ENV0, names(d))] <- make.names(ENV0, unique = TRUE); ENV <- make.names(ENV0, unique = TRUE)
p <- strsplit(rownames(W), "_")
d$st     <- vapply(p, function(x) paste(tail(x, 2), collapse = "_"), "")
d$leg    <- vapply(p, `[`, "", 2)
d$season <- vapply(p, `[`, "", 3)
d$depth_zone <- as.character(cut(d$g_Depth, c(-Inf, 200, 1000, Inf), labels = c("shelf", "upper slope", "lower slope")))
SV <- c("2019", "2021"); ZN <- c("shelf", "upper slope", "lower slope")
zone_of <- tapply(d$depth_zone, d$st, `[`, 1)

## ---- W1 transfer vs interpolation in time, by depth zone ---------------------------------------------------------
T4a <- read.csv(file.path(OUT, "T4a_time_rows.csv"))
T4a$depth_zone <- zone_of[T4a$st]
boot <- function(x, B = 2000) { set.seed(SEED); b <- replicate(B, mean(sample(x, replace = TRUE))); c(mean(x), quantile(b, c(0.025, 0.975))) }
w1 <- do.call(rbind, lapply(ZN, function(z) { x <- T4a[T4a$depth_zone == z, ]
  tg <- boot(x$T - x$G); g5 <- boot(x$G50 - x$G)
  data.frame(depth_zone = z, n_station_surveys = nrow(x), G_interpolation = mean(x$G), T_transfer = mean(x$T),
             T_minus_G = tg[1], TG_lo = tg[2], TG_hi = tg[3], T_minus_G_pct_of_G = 100 * tg[1] / mean(x$G),
             G50_minus_G = g5[1], G50G_lo = g5[2], G50G_hi = g5[3]) }))
write.csv(w1, file.path(OUT, "W1_time_transfer_by_zone.csv"), row.names = FALSE)
print(w1, digits = 3, row.names = FALSE)

## ---- W2 environmental change in time vs in space, by depth zone --------------------------------------------------
Zs  <- scale(as.matrix(d[, ENV]), scale = pmax(apply(d[, ENV], 2, sd), 1e-9))
ms  <- lapply(BIO, function(y) ranger(x = d[, ENV], y = d[[y]], num.trees = 300, min.node.size = 5,
                                      importance = "impurity", num.threads = NT, seed = SEED))
imp <- rowMeans(sapply(ms, function(m) { v <- pmax(m$variable.importance[ENV], 0); v / max(sum(v), 1e-12) }))
Zw  <- sweep(Zs, 2, imp, "*")
km  <- function(i) cbind(d$g_Long[i] * 111.32 * cos(mean(d$g_Lat) * pi / 180), d$g_Lat[i] * 110.57)
both <- names(which(table(d$st) == 2))
r1 <- vapply(both, function(q) which(d$st == q & d$survey == 1), 1L); r2 <- vapply(both, function(q) which(d$st == q & d$survey == 2), 1L)
e_time <- sqrt(rowSums((Zw[r1, ] - Zw[r2, ])^2))
nn_env <- rep(NA_real_, nrow(d))
for (s in 1:2) { i <- which(d$survey == s); nn <- get.knnx(km(i), km(i), k = 2)$nn.index[, 2]; nn_env[i] <- sqrt(rowSums((Zw[i, ] - Zw[i[nn], ])^2)) }
st_tab <- data.frame(st = both, depth_zone = zone_of[both], depth_m = d$g_Depth[r1], env_change_between_surveys = e_time,
                     env_distance_nearest_station = (nn_env[r1] + nn_env[r2]) / 2)
st_tab$ratio_time_to_space <- st_tab$env_change_between_surveys / st_tab$env_distance_nearest_station
w2 <- do.call(rbind, lapply(c(ZN, "all"), function(z) { x <- if (z == "all") st_tab else st_tab[st_tab$depth_zone == z, ]
  data.frame(depth_zone = z, n_stations = nrow(x), median_env_change_between_surveys = median(x$env_change_between_surveys),
             median_env_distance_nearest_station = median(x$env_distance_nearest_station),
             median_ratio_time_to_space = median(x$ratio_time_to_space),
             pct_stations_time_gt_space = 100 * mean(x$ratio_time_to_space > 1)) }))
kw <- kruskal.test(env_change_between_surveys ~ depth_zone, data = st_tab)
w2$kruskal_p_change_by_zone <- kw$p.value
write.csv(st_tab, file.path(OUT, "W2_time_change_stations.csv"), row.names = FALSE)
write.csv(w2, file.path(OUT, "W2_time_change_by_zone.csv"), row.names = FALSE)
print(w2, digits = 3, row.names = FALSE)

## ---- W3 which variables change between surveys -------------------------------------------------------------------
sq <- (Zw[r1, ] - Zw[r2, ])^2
w3 <- data.frame(variable = ENV0, importance_weight = imp, share_of_change_pct = 100 * colSums(sq) / sum(sq),
                 median_abs_change_sd = apply(abs(Zs[r1, ] - Zs[r2, ]), 2, median))
w3 <- w3[order(-w3$share_of_change_pct), ]
w3$cumulative_pct <- cumsum(w3$share_of_change_pct)
write.csv(w3, file.path(OUT, "W3_time_change_variables.csv"), row.names = FALSE)
print(head(w3, 12), digits = 3, row.names = FALSE)

## ---- W4 season and cruise leg of every sample -----------------------------------------------------------------------
w4 <- as.data.frame(table(survey = SV[d$survey], leg = d$leg, season = d$season, depth_zone = d$depth_zone))
w4 <- w4[w4$Freq > 0, ]; names(w4)[names(w4) == "Freq"] <- "samples"
write.csv(w4, file.path(OUT, "W4_season_leg_by_survey.csv"), row.names = FALSE)
print(w4, row.names = FALSE)
sessioninfo::session_info(to_file = file.path(OUT, "sessionInfo_when.txt"))
