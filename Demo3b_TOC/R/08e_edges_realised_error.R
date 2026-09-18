# =============================================================================
# 08e_edges_realised_error.R — did REALISED error also rise where the ignorance map now
# expects more error (lithology boundaries, ridge-axis class, rough deep seafloor)?
# =============================================================================
# The red features in the delta map are changes in EXPECTED error. Here the same withheld points
# (nested leave-region-out, 07_calibration.R) are compared between iterations: paired change in
# absolute error, and change in expected error, by feature. Points are spatially autocorrelated,
# so the normal-approximation CIs are optimistic.
Sys.setenv(DEMO3B_ITER = "iter2_texture")
source("R/00_config.R")

stk <- rast(path_toc_stack())
dep <- stk[["depth"]]; if (mean(values(dep), na.rm = TRUE) < 0) dep <- -dep
L <- stk[[paste0("litho_t", 1:6)]]; dom <- which.max(L)
m <- as.matrix(dom, wide = TRUE)
nb <- function(a, dr, dc) { out <- matrix(NA, nrow(a), ncol(a))
  r <- max(1, 1 + dr):min(nrow(a), nrow(a) + dr); cc <- max(1, 1 + dc):min(ncol(a), ncol(a) + dc)
  out[r, cc] <- a[r - dr, cc - dc]; out }
S4 <- list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
edge <- matrix(FALSE, nrow(m), ncol(m))
for (s in S4) { n <- nb(m, s[1], s[2]); edge <- edge | (!is.na(n) & !is.na(m) & n != m) }
dil <- function(e) { o <- e; for (s in S4) { n <- nb(e, s[1], s[2]); o <- o | (!is.na(n) & n) }; o }
er <- rast(dom); values(er) <- as.vector(t(dil(dil(edge))))

rd <- function(it) read.csv(file.path(ROOT, "outputs", it, "tables", "calibration_outer_points_toc.csv"))
p1 <- rd(Sys.getenv("CMP_FROM", "iter1b_supply")); p2 <- rd(Sys.getenv("CMP_TO", "iter2b_supply_texture"))
key <- function(d) paste(round(d$lon, 3), round(d$lat, 3))
a <- data.frame(k = key(p1), lon = p1$lon, lat = p1$lat, ae1 = abs(p1$pred - p1$obs), ex1 = p1$exp_C2)
b <- data.frame(k = key(p2), ae2 = abs(p2$pred - p2$obs), ex2 = p2$exp_C2)
p <- merge(a, b, by = "k")
msg("paired withheld points: %d", nrow(p))
xy <- cbind(p$lon, p$lat)
p$edge <- extract(er, xy)[, 1] == 1; p$t6 <- extract(L[[6]], xy)[, 1] > 0.5
p$depth <- extract(dep, xy)[, 1]; p$tri <- extract(stk[["tri"]], xy)[, 1]
dv <- values(dep)[, 1]; tri90 <- quantile(values(stk[["tri"]])[which(dv > 1500), 1], 0.9, na.rm = TRUE)

grp <- list(`all withheld points` = rep(TRUE, nrow(p)),
            `within 0.2 deg of a lithology boundary` = p$edge, `away from lithology boundaries` = !p$edge,
            `boundary, deep water (>1500 m)` = p$edge & p$depth > 1500, `boundary, shelf & slope` = p$edge & p$depth <= 1500,
            `lithology type 6 dominant` = p$t6, `roughest 10% of deep seafloor` = p$depth > 1500 & p$tri > tri90)
res <- do.call(rbind, lapply(names(grp), function(g) { s <- which(grp[[g]] %in% TRUE); x <- p$ae2[s] - p$ae1[s]
  data.frame(group = g, n = length(s),
             RMSE_iter1 = sqrt(mean(p$ae1[s]^2)), RMSE_iter2 = sqrt(mean(p$ae2[s]^2)),
             mean_change_abs_error = mean(x), ci95_lo = mean(x) - 1.96 * sd(x) / sqrt(length(s)), ci95_hi = mean(x) + 1.96 * sd(x) / sqrt(length(s)),
             mean_change_expected_error = mean(p$ex2[s] - p$ex1[s]),
             wilcoxon_p = if (length(s) > 5) wilcox.test(p$ae2[s], p$ae1[s], paired = TRUE)$p.value else NA_real_) }))
write.csv(res, file.path(DIRS$compare, "lithology_edges_realised_error.csv"), row.names = FALSE)
print(res, digits = 3)
msg("realised-error check done")
