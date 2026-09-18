# =============================================================================
# 08l_collection_by_region.R — per region: which path, and how to collect (directed or random; which space)
# =============================================================================
# Joins the decision rule (08g) with the directed-sampling test (08h) for the same iteration.
# For each region and strategy, the advantage over random is averaged over the two effort levels (40%, 60%)
# and the three random starts. A strategy is recommended only if its mean advantage is positive at both effort
# levels; among those, the one with the largest mean advantage. Otherwise random collection suffices.
#   novelty     -> directed in environmental space (DI logic)
#   uncertainty -> directed by model uncertainty (expected error)
#   space       -> directed in geographic space
Sys.setenv(DEMO3B_ITER = Sys.getenv("DEMO3B_ITER", "iter1b_supply"))
source("R/00_config.R")
rule <- read.csv(file.path(DIRS$tables, "regionalization_by_region_toc.csv"))
h7   <- read.csv(file.path(DIRS$tables, "targeted_sampling_by_region_toc.csv"))
LAB  <- c(novelty = "directed: environmental space", uncertainty = "directed: model uncertainty", space = "directed: geographic space")
rec <- do.call(rbind, lapply(split(h7, h7$region), function(x) {
  s <- do.call(rbind, lapply(names(LAB), function(k) { y <- x[x$strategy == k, ]
    data.frame(strategy = k, adv40 = y$gain_minus_random[y$pct_pool == 40], adv60 = y$gain_minus_random[y$pct_pool == 60],
               gain_random = mean(x$gain[x$strategy == "random"])) }))
  s$mean_adv <- (s$adv40 + s$adv60) / 2
  ok <- s[s$adv40 > 0 & s$adv60 > 0, ]
  best <- if (nrow(ok)) ok[which.max(ok$mean_adv), ] else NULL
  data.frame(region = x$region[1],
             collection = if (is.null(best)) "random suffices" else LAB[[best$strategy]],
             advantage_vs_random = if (is.null(best)) NA else best$mean_adv,
             advantage_pct_of_random_gain = if (is.null(best)) NA else 100 * best$mean_adv / s$gain_random[1],
             adv_environmental = s$mean_adv[s$strategy == "novelty"], adv_uncertainty = s$mean_adv[s$strategy == "uncertainty"],
             adv_geographic = s$mean_adv[s$strategy == "space"])
}))
out <- merge(rule[, c("region", "RMSE_transfer_T", "RMSE_global_G", "RMSE_regional_R", "T_minus_G", "G_minus_R", "G50_minus_G", "decision")], rec, by = "region")
out$path <- ifelse(grepl("^collect", out$decision), "1 collect data and load into the global model",
            ifelse(grepl("^regionalize", out$decision), "2 regionalize the model",
            ifelse(grepl("^new drivers", out$decision), "3 new drivers or finer resolution", "stop: noise floor")))
write.csv(out, file.path(DIRS$tables, "collection_by_region_toc.csv"), row.names = FALSE)
print(out[, c("region", "path", "collection", "advantage_vs_random", "advantage_pct_of_random_gain", "adv_environmental", "adv_uncertainty", "adv_geographic")], digits = 3, row.names = FALSE)
