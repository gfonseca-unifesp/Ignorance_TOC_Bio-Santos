# =============================================================================
# 08n_space_env_structure.R — is the environment structured in geographic space?
# =============================================================================
# If conditions change smoothly with distance, filling geographic gaps and filling environmental gaps pick
# similar sites, and their gains coincide. Within each region of the H7 experiment (08h, same splits):
#   rho_geo_env   Spearman correlation between pairwise great-circle distance and pairwise distance in the
#                 importance-weighted standardised predictor space (the DI space); <= 1,500 cells per region
#   jaccard       overlap between the blocks picked by geographic gap filling ("space") and by environmental
#                 novelty at the first step of 08h (same test blocks and starts), mean of the 3 starts,
#                 and the overlap expected between two random picks of the same size
#   gains         extra gain over random of each strategy at 40% of the pool (targeted_sampling_by_region_toc.csv)
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1b_supply"))
source("R/00_config.R")
suppressPackageStartupMessages({ library(FNN) })

model <- readRDS(file.path(DIRS$models, "rf_spatialCV_toc.rds"))
vars  <- setdiff(names(model$trainingData), c(".outcome", ".weights"))
d     <- read.csv(file.path(DIRS$data_iter, "toc_training_table.csv"))
to_xyz   <- function(lon, lat) { lo <- lon * pi / 180; la <- lat * pi / 180; cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) }
chord_km <- function(z) 2 * 6371 * asin(pmin(1, z / 2))
set.seed(CFG$seed)
d$region <- kmeans(to_xyz(d$x_lon, d$y_lat), centers = CFG$n_regions, nstart = 25, iter.max = 100)$cluster   # same as 08h
exy <- st_coordinates(st_transform(st_as_sf(d, coords = c("x_lon", "y_lat"), crs = 4326), "EPSG:8857"))
d$rblock <- paste(d$region, floor(exy[, 1] / 3e5), floor(exy[, 2] / 3e5))
RN <- c("1 S Atlantic-Scotia", "2 NE Atlantic", "3 Arabian Sea", "4 Angola-Benguela", "5 Siberian Arctic", "6 NW Atlantic", "7 NE Pacific", "8 W Pacific")
imp <- pmax(model$finalModel$variable.importance[vars], 0); imp <- imp / sum(imp)
Zw  <- sweep(scale(as.matrix(d[, vars])), 2, imp, "*")
X3  <- to_xyz(d$x_lon, d$y_lat)

res <- do.call(rbind, lapply(sort(unique(d$region)), function(r) {
  inr <- which(d$region == r); outr <- which(d$region != r)
  blk <- unique(d$rblock[inr]); set.seed(CFG$seed + r)
  test_b <- sample(blk, round(0.2 * length(blk))); pool <- setdiff(blk, test_b)
  ov <- sapply(1:3, function(rep) {
    set.seed(CFG$seed + 10 * r + rep)
    start_b <- sample(pool, max(1, round(0.2 * length(pool))))
    cand <- setdiff(pool, start_b); n_add <- round(0.4 * length(pool)) - length(start_b)
    sampled_idx <- c(outr, inr[d$rblock[inr] %in% start_b])
    ci <- inr[d$rblock[inr] %in% cand]; cb <- d$rblock[ci]
    dmin <- chord_km(get.knnx(X3[sampled_idx, , drop = FALSE], X3[ci, , drop = FALSE], k = 1)$nn.dist[, 1]); sp <- character(0)
    for (j in seq_len(n_add)) {
      bs <- tapply(dmin, cb, mean); bs <- bs[setdiff(names(bs), sp)]; b <- names(bs)[which.max(bs)]; sp <- c(sp, b)
      dmin <- pmin(dmin, chord_km(get.knnx(X3[ci[cb == b], , drop = FALSE], X3[ci, , drop = FALSE], k = 1)$nn.dist[, 1]))
    }
    nv <- names(sort(tapply(get.knnx(Zw[sampled_idx, , drop = FALSE], Zw[ci, , drop = FALSE], k = 1)$nn.dist[, 1], cb, mean), decreasing = TRUE))[seq_len(n_add)]
    e <- n_add^2 / length(cand)
    c(jaccard = length(intersect(sp, nv)) / length(union(sp, nv)), jaccard_random = e / (2 * n_add - e), n_add = n_add, n_cand = length(cand))
  })
  set.seed(CFG$seed + 1000 + r); s <- if (length(inr) > 1500) sample(inr, 1500) else inr
  data.frame(region = RN[r], n_cells = length(inr), rho_geo_env = cor(chord_km(as.vector(dist(X3[s, ]))), as.vector(dist(Zw[s, ])), method = "spearman"),
             jaccard_space_novelty = mean(ov["jaccard", ]), jaccard_random = mean(ov["jaccard_random", ]), blocks_added = mean(ov["n_add", ]))
}))
h7 <- read.csv(file.path(DIRS$tables, "targeted_sampling_by_region_toc.csv"))
h7 <- h7[h7$pct_pool == 40, ]
res$gain_vs_random_space   <- h7$gain_minus_random[h7$strategy == "space"][match(res$region, h7$region[h7$strategy == "space"])]
res$gain_vs_random_novelty <- h7$gain_minus_random[h7$strategy == "novelty"][match(res$region, h7$region[h7$strategy == "novelty"])]
res$abs_diff_space_novelty <- abs(res$gain_vs_random_space - res$gain_vs_random_novelty)
write.csv(res, file.path(DIRS$tables, "space_env_structure_toc.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
msg("rho(structure, |space - novelty|) across regions = %.2f; rho(structure, overlap) = %.2f",
    cor(res$rho_geo_env, res$abs_diff_space_novelty, method = "spearman"), cor(res$rho_geo_env, res$jaccard_space_novelty, method = "spearman"))
