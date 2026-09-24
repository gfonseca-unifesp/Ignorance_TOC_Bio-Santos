# =============================================================================
# 09_audit_figure.R — conceptual figure: (a) ignorance audit loop, (b) ignorance budget
# =============================================================================
# (a) The loop: measure the ignorance map -> test it against withheld regions -> diagnose the
#     type of mismatch -> act by ignorance level -> re-measure the gain. Annotated with the
#     seafloor-TOC case (numbers read from the outputs). Presentation sequence: base model with carbon
#     supply (iter1b_supply) -> + sediment texture (iter2b_supply_texture).
# (b) The ignorance budget: mean squared error of withheld-region predictions (nested
#     leave-region-out, 07_calibration.R) decomposed additively:
#       MSE = mean(regional mean error^2)                                  regional offset (transfer)
#           + mean((mean error by region x sediment type - regional mean)^2)  bias by sediment type (L2)
#           + mean((error - mean error by region x sediment type)^2)          within-group error,
#     split into the noise floor (L5, within-cell SD of observations) and the remainder
#     (L4 / not yet named). One bar per analytic iteration.

source("R/00_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
dir.create(DIRS$compare, recursive = TRUE, showWarnings = FALSE)

IT1 <- Sys.getenv("CMP_FROM", "iter1b_supply"); IT2 <- Sys.getenv("CMP_TO", "iter2b_supply_texture")
ITS <- setNames(c("Iteration 1\nbase model (with carbon supply)", if (grepl("clean", IT2)) "Iteration 2\n+ lithology" else "Iteration 2\n+ sediment texture"), c(IT1, IT2))
tabp <- function(it, f) { p <- file.path(ROOT, "outputs", it, "tables", f); if (file.exists(p)) read.csv(p) else NULL }
getv <- function(s, m) { v <- s$value[s$metric == m]; if (length(v)) v else NA }
noise <- median(read.csv(file.path(DIRS$data, "toc_cells.csv"))$y_sd, na.rm = TRUE)
tex <- rast(path_texture_stack())[[paste0("litho_t", 1:6)]]
lith_class <- function(lon, lat) {
  v <- as.matrix(extract(tex, cbind(lon, lat))); s <- rowSums(v)
  k <- max.col(v, ties.method = "first"); ifelse(is.na(s) | s < 0.5, "unclassified", paste0("type ", k))
}
COMP <- c("Noise floor (L5)", "Unstructured error above noise (L4 / not yet named)",
          "Bias by sediment type within regions (L2)", "Regional offset (transfer to unsampled regions)")

# --- (b) budget ----------------------------------------------------------------------------
avail <- names(ITS)[vapply(names(ITS), function(it) !is.null(tabp(it, "calibration_outer_points_toc.csv")), TRUE)]
budget <- do.call(rbind, lapply(avail, function(it) {
  x <- tabp(it, "calibration_outer_points_toc.csv")
  x$err <- x$pred - x$obs; x$lith <- lith_class(x$lon, x$lat)
  reg  <- ave(x$err, x$region)
  regl <- ave(x$err, x$region, x$lith)
  mse <- mean(x$err^2); within <- mean((x$err - regl)^2); nf <- min(noise^2, within)
  data.frame(iteration = it, component = COMP,
             mse = c(nf, within - nf, mean((regl - reg)^2), mean(reg^2)),
             total_mse = mse, rmse = sqrt(mse), n = nrow(x))
}))
budget$share <- 100 * budget$mse / budget$total_mse
write.csv(budget, file.path(DIRS$compare, "ignorance_budget.csv"), row.names = FALSE)
print(budget, digits = 3)
stopifnot(all(abs(tapply(budget$mse, budget$iteration, sum) - tapply(budget$total_mse, budget$iteration, `[`, 1)) < 1e-9))

cal <- do.call(rbind, lapply(avail, function(it) { s <- tabp(it, "toc_summary.csv"); m <- tabp(it, "calibration_metrics_toc.csv")
  b <- m[m$method == getv(s, "calibration_best_method") & m$subset == "all withheld points", ]
  data.frame(iteration = it, R2_cv = as.numeric(getv(s, "R2_log_spatialCV")), coverage = b$coverage90, rho = b$spearman_points,
             rho_bins = b$spearman_bins, slope = b$slope) }))
tot <- tapply(budget$total_mse, budget$iteration, `[`, 1)
comp_mse <- function(it, pat) budget$mse[budget$iteration == it & grepl(pat, budget$component)]
pct_change <- function(pat) 100 * (comp_mse(IT2, pat) / comp_mse(IT1, pat) - 1)

budget$iteration_lab <- factor(ITS[budget$iteration], levels = rev(ITS))
budget$component <- factor(budget$component, levels = COMP)
pal <- setNames(c("#B4B2A9", "#F0997B", "#7F77DD", "#5DCAA5"), COMP)
lab_end <- merge(unique(budget[, c("iteration", "iteration_lab", "total_mse", "rmse")]), cal, by = "iteration")
lab_end$txt <- sprintf("withheld RMSE %.3f | spatial-CV R2 %.2f\n90%% interval coverage %.2f\ncalibration across dissimilarity bins rho %.2f;\nranking of individual points rho %.2f",
                       lab_end$rmse, lab_end$R2_cv, lab_end$coverage, lab_end$rho_bins, lab_end$rho)
pb <- ggplot(budget, aes(mse, iteration_lab, fill = component)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(share >= 4, sprintf("%.0f%%", share), "")), position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.6) +
  geom_text(data = lab_end, aes(x = total_mse, y = iteration_lab, label = txt), inherit.aes = FALSE, hjust = -0.04, size = 2.4, lineheight = 0.95) +
  # MSv10 (F3.3): os rótulos da legenda usam os termos do texto (Table 3, Box 3); os níveis do fator,
  # que vêm da tabela do orçamento, não mudam
  scale_fill_manual(values = pal, name = NULL,
                    labels = c("Replicate floor (L5)", "Unstructured remainder (L4; candidate L2)",
                               "Bias by sediment type within regions (L2)",
                               "Regional offset (transfer to unsampled regions)")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.7))) +
  labs(x = "mean squared error of withheld-region predictions (log10 TOC)^2", y = NULL,
       title = "b  The ignorance budget: what each iteration removed, and what remains") +
  guides(fill = guide_legend(ncol = 2)) +
  theme_classic(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 9))

# --- (a) loop ---------------------------------------------------------------------------------
c1 <- cal[cal$iteration == IT1, ]; c2 <- if (IT2 %in% cal$iteration) cal[cal$iteration == IT2, ] else NULL
XC <- c(14.5, 38.5, 62.5, 86.5)                     # column centres; x = 1.5 is kept free for the loop path
steps <- data.frame(x = XC, y = 52, w = 20, h = 10,
  title = c("1  Measure", "2  Test", "3  Diagnose", "4  Act by level"),
  sub = c("ignorance map\n(DI, distance, error model)", "against withheld\nregions", "type the mismatch\n(hypotheses by level)", "the gap's level sets\nthe response"))
levels_df <- data.frame(x = XC, y = 30, w = 20, h = 12,
  title = c("L2  structure", "L3  sampling", "L4  model", "L5  noise"),
  sub = c("integrate data,\nrevise conceptual model", "collect where flagged", "better representation", "report the floor"),
  hl = c(TRUE, TRUE, FALSE, FALSE))
remeasure <- data.frame(x = 50.5, y = 11, w = 23, h = 9, title = "5  Re-measure", sub = "quantify the knowledge gain")
case <- data.frame(
  x = c(XC[2], XC[3], XC[1], XC[2], 81.5),
  y = c(43.5, 43.5, 25.8, 25.8, 11),
  txt = c(sprintf("TOC: error magnitude calibrated\n(90%% interval coverage %.2f)", c1$coverage),
          "H1 dissimilarity: weak\nH2 noise: rejected\nH3 missing drivers",
          if (grepl("clean", IT2)) "TOC: + carbon supply,\n+ sample-based lithology" else "TOC: + carbon supply,\n+ sediment texture",
          { rg <- tabp(IT1, "regionalization_by_region_toc.csv")
            if (is.null(rg)) "TOC: collect data,\nstrategy chosen by test" else
              sprintf("TOC: collect in %d of %d regions,\nstrategy chosen by test", sum(!grepl("^new drivers|^stop", rg$decision)), nrow(rg)) },   # marginal regionalize calls read as collect (08v, Text S2.7)
          if (is.null(c2)) "TOC iteration 2: pending" else
            sprintf("TOC iteration 2: withheld MSE %+.1f%%\n(regional offset %+.0f%%, sediment bias %+.0f%%)\nremaining: unstructured error above noise",
                    100 * (tot[[IT2]] / tot[[IT1]] - 1), pct_change("Regional offset"), pct_change("Bias by sediment"))))
box <- function(df, fill, col, lw = 0.3) list(
  geom_rect(data = df, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2), fill = fill, colour = col, linewidth = lw, inherit.aes = FALSE),
  geom_text(data = df, aes(x = x, y = y + h / 2 - 2.2, label = title), fontface = "bold", size = 2.7, inherit.aes = FALSE),
  geom_text(data = df, aes(x = x, y = y + (h - 10) / 2 - 1.6, label = sub), size = 2.3, lineheight = 0.9, colour = "grey25", inherit.aes = FALSE))
ar <- arrow(length = unit(1.6, "mm"), type = "closed")
top_lv <- 30 + 12 / 2; bot_lv <- 30 - 12 / 2
pa <- ggplot() +
  box(steps, "#F1EFE8", "#5F5E5A") +
  box(levels_df[!levels_df$hl, ], "#EEEDFE", "#534AB7") + box(levels_df[levels_df$hl, ], "#CECBF6", "#26215C", lw = 0.9) +
  box(remeasure, "#F1EFE8", "#5F5E5A") +
  geom_segment(data = data.frame(x = XC[1:3] + 10, xend = XC[2:4] - 10.2, y = 52), aes(x = x, xend = xend, y = y, yend = y), arrow = ar, linewidth = 0.3) +
  geom_segment(aes(x = XC[4], xend = XC[4], y = 47, yend = 40), linewidth = 0.3) +
  geom_segment(aes(x = XC[1], xend = XC[4], y = 40, yend = 40), linewidth = 0.3) +
  geom_segment(data = data.frame(x = XC), aes(x = x, xend = x, y = 40, yend = top_lv + 0.2), arrow = ar, linewidth = 0.3) +
  geom_segment(data = data.frame(x = XC), aes(x = x, xend = x, y = bot_lv, yend = 19.5), linewidth = 0.3) +
  geom_segment(aes(x = XC[1], xend = XC[4], y = 19.5, yend = 19.5), linewidth = 0.3) +
  geom_segment(aes(x = 50.5, xend = 50.5, y = 19.5, yend = 15.7), arrow = ar, linewidth = 0.3) +
  geom_path(data = data.frame(x = c(39, 1.5, 1.5, 4.3), y = c(11, 11, 52, 52)), aes(x = x, y = y), arrow = ar, linewidth = 0.35, colour = "grey35") +
  annotate("text", x = 20.5, y = 13, label = "next iteration", size = 2.3, colour = "grey35") +
  geom_text(data = case, aes(x = x, y = y, label = txt), size = 2.2, fontface = "italic", colour = "#3C3489", lineheight = 0.9, hjust = 0.5) +
  coord_cartesian(xlim = c(0, 100), ylim = c(6, 59), expand = FALSE) + theme_void(base_size = 8) +
  labs(title = "a  The ignorance audit loop, with the seafloor-TOC case (italic)") +
  theme(plot.title = element_text(face = "bold", size = 9))

fig <- (pa / pb) + plot_layout(heights = c(1.25, 0.7))
ggsave(file.path(DIRS$compare, "fig_ignorance_audit_loop_budget.png"), fig, width = 200, height = 195, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(DIRS$compare, "fig_ignorance_audit_loop_budget.pdf"), fig, width = 200, height = 195, units = "mm", bg = "white")
msg("audit figure written (%s)", paste(avail, collapse = ", "))
