# =============================================================================
# 08x_structure_confirmation.R — are the structured partitions confirmed in blocks not used to discover them?
# =============================================================================
# Discovery/confirmation split of the error-structure test (08w). The 300-km blocks are split at random into two halves
# (all units of a block stay in the same half, so autocorrelated cells never cross halves).
#   discovery half:    test all declared partitions x responses; structured = eta2 >= 0.01 and BH q < 0.05
#   confirmation half: re-test only what the discovery half flagged; confirmed = eta2 >= 0.01 and p < 0.05
# Repeated over 100 random splits. Reported per partition: how often it is discovered, how often a discovery is confirmed,
# and the median eta2 in each half; plus one illustrative split (the first).
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1c_clean"))
source("R/00_config.R")
ITS <- c("iter1c_clean", "iter2c_clean"); NSPLIT <- 100; NPERM <- 499
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
to_xyz <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
tex <- rast(path_texture_stack())[[paste0("litho_t", 1:6)]]
eta2 <- function(y, g) { m <- tapply(y, g, mean); n <- tapply(y, g, length); sum(n * (m - mean(y))^2) / sum((y - mean(y))^2) }
ptest <- function(y, g) { o <- eta2(y, g); c(eta2 = o, p = (1 + sum(replicate(NPERM, eta2(y, sample(g))) >= o)) / (NPERM + 1)) }

out <- list(); ill <- list()
for (it in ITS) {
  m <- readRDS(file.path(ROOT, "outputs", it, "models", "rf_spatialCV_toc.rds")); vars <- setdiff(names(m$trainingData), c(".outcome", ".weights"))
  d <- read.csv(file.path(ROOT, "outputs", it, "data", "toc_training_table.csv"))
  d$pred_cv <- NA_real_; d$pred_cv[m$pred$rowIndex] <- m$pred$pred; d$e <- d$pred_cv - d$y; d$ae <- abs(d$e)
  set.seed(CFG$seed); d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster
  exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
  d$blk <- paste(floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
  lv <- as.matrix(extract(tex, cbind(d$x_lon, d$y_lat))); ls <- rowSums(lv)
  P <- list(region = RN[d$region],
            latitude_band = as.character(cut(d$y_lat, seq(-90, 90, 30), include.lowest = TRUE)),
            longitude_sector = as.character(cut(d$x_lon, seq(-180, 180, 60), include.lowest = TRUE)),
            depth_zone = as.character(cut(d$depth, c(-Inf, 200, 1000, 3000, Inf))),
            lithology_class = ifelse(is.na(ls) | ls < 0.5, "unclassified", paste0("type ", max.col(lv, ties.method = "first"))))
  U <- lapply(P, function(g) { u <- aggregate(cbind(bias = d$e, magnitude = d$ae) ~ blk + level, data = data.frame(e = d$e, ae = d$ae, blk = d$blk, level = g), FUN = mean)
    cnt <- table(paste(d$blk, g, sep = "|")); u$n <- as.numeric(cnt[paste(u$blk, u$level, sep = "|")]); u[u$n >= 5, ] })
  blocks <- unique(d$blk)
  set.seed(CFG$seed)
  for (s in seq_len(NSPLIT)) {
    disc_blk <- sample(blocks, floor(length(blocks) / 2))
    rows <- list()
    for (pn in names(P)) for (resp in c("bias", "magnitude")) {
      u <- U[[pn]]; dh <- u[u$blk %in% disc_blk, ]; ch <- u[!u$blk %in% disc_blk, ]
      a <- ptest(dh[[resp]], dh$level)
      rows[[length(rows) + 1]] <- data.frame(iteration = it, split = s, partition = pn, response = resp, eta2_discovery = a[["eta2"]], p_discovery = a[["p"]],
                                             ch = I(list(ch[, c("level", resp)])))
    }
    R <- do.call(rbind, rows); R$q_discovery <- p.adjust(R$p_discovery, "BH")
    R$discovered <- R$eta2_discovery >= 0.01 & R$q_discovery < 0.05
    R$eta2_confirmation <- NA_real_; R$p_confirmation <- NA_real_
    for (i in which(R$discovered)) { ch <- R$ch[[i]]; b <- ptest(ch[[2]], ch$level); R$eta2_confirmation[i] <- b[["eta2"]]; R$p_confirmation[i] <- b[["p"]] }
    R$confirmed <- R$discovered & R$eta2_confirmation >= 0.01 & R$p_confirmation < 0.05
    R$ch <- NULL; out[[length(out) + 1]] <- R
  }
  msg("%s: %d splits done", it, NSPLIT)
}
A <- do.call(rbind, out)
S <- do.call(rbind, lapply(split(A, list(A$iteration, A$partition, A$response), drop = TRUE), function(x) data.frame(
  iteration = x$iteration[1], partition = x$partition[1], response = x$response[1],
  pct_splits_discovered = 100 * mean(x$discovered), pct_discoveries_confirmed = if (any(x$discovered)) 100 * mean(x$confirmed[x$discovered]) else NA,
  median_eta2_discovery = median(x$eta2_discovery), median_eta2_confirmation = if (any(x$discovered)) median(x$eta2_confirmation[x$discovered]) else NA)))
S <- S[order(S$iteration, S$partition, S$response), ]
write.csv(A, file.path(DIRS$compare, "error_structure_confirmation_splits.csv"), row.names = FALSE)
write.csv(S, file.path(DIRS$compare, "error_structure_confirmation_summary.csv"), row.names = FALSE)
print(S, digits = 3, row.names = FALSE)
cat("\n--- illustrative split (split 1, iteration 1) ---\n")
print(A[A$split == 1 & A$iteration == "iter1c_clean", c("partition", "response", "eta2_discovery", "q_discovery", "discovered", "eta2_confirmation", "p_confirmation", "confirmed")], digits = 3, row.names = FALSE)
msg("structure confirmation done")
