## ================================================================
## fig1_data_information_plane.R — Figure 1 of the manuscript, from a script
## ================================================================
## Until MSv9 this figure was a drawing without a generator, which the pre-submission review flagged
## (minor point 17: wasted space on the right, unlabelled crossing arrows, no link to the levels of
## ignorance). This script rebuilds it in ggplot2:
##   - the four quadrants of the data-information plane (Q1-Q4);
##   - the three moves between them, each labelled, and the move that goes backwards (no-analogue);
##   - the level of ignorance (L0-L5) that dominates in each quadrant, which ties Figure 1 to Table 1.
## Text sits in the upper part of each quadrant and the arrows in the free bands between, so that
## nothing overlaps. Output: outputs/fig1_data_information_plane.png and .pdf, 180 x 115 mm, 300 dpi.
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" figures_src/fig1_data_information_plane.R
suppressPackageStartupMessages(library(ggplot2))
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/Gustavo/OneDrive/Documentos/Ignorance_MS")
OUT  <- file.path(ROOT, "figures_src", "outputs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FILL <- c(Q1 = "#DCE6F1", Q2 = "#E2EFE2", Q3 = "#FCE7CF", Q4 = "#EFEFEF")
INK  <- c(Q1 = "#1F3864", Q2 = "#2E6B45", Q3 = "#B5651D", Q4 = "#4D4D4D")
GREEN <- "#2E8B6A"; BLUE <- "#2E5496"; RED <- "#B23A48"

quad <- data.frame(
  q     = c("Q2", "Q1", "Q4", "Q3"),
  xmin  = c(0, 5, 0, 5), xmax = c(5, 10, 5, 10),
  ymin  = c(5, 5, 0, 0),  ymax = c(10, 10, 5, 5),
  title = c("Q2 \u00b7 Knowledge without data", "Q1 \u00b7 Well constrained",
            "Q4 \u00b7 Blank", "Q3 \u00b7 Data without knowledge"),
  sub   = c("theory or local knowledge carry it", "data and knowledge agree",
            "no data, no model, no question", "measured, but not answering"),
  lev   = c("gap at L3: unmeasured here and now", "gap at L4, then the L5 floor",
            "gap at L0\u2013L1: no model, no question", "gap at L1\u2013L2: unnamed or dispersed"),
  stringsAsFactors = FALSE)
quad$lx <- quad$xmin + 0.3
# Q3 carries its heading a little lower, to leave the foot of the diagonal arrow clear
quad$ty <- quad$ymax - ifelse(quad$q == "Q3", 1.15, 0.5)
# the level line sits at the foot of each quadrant, except in Q1, where the diagonal arrow arrives
quad$ly <- ifelse(quad$q == "Q1", quad$ymax - 1.5, quad$ymin + 0.35)

p <- ggplot() +
  geom_rect(data = quad, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = q),
            colour = "white", linewidth = 1.2) +
  scale_fill_manual(values = FILL, guide = "none") +
  geom_text(data = quad, aes(lx, ty, label = title, colour = q),
            hjust = 0, fontface = "bold", size = 3.4, show.legend = FALSE) +
  scale_colour_manual(values = INK, guide = "none") +
  geom_text(data = quad, aes(lx, ty - 0.5, label = sub), hjust = 0, size = 2.8, colour = "grey25") +
  geom_text(data = quad, aes(lx, ly, label = lev), hjust = 0, size = 2.6,
            fontface = "italic", colour = "grey35") +

  ## the moves between quadrants, in the bands the text leaves free
  # Q4 -> Q3: more data, no knowledge gained
  annotate("segment", x = 0.9, xend = 7.2, y = 1.35, yend = 1.35, linewidth = 0.8, colour = BLUE,
           linetype = "22", arrow = arrow(length = unit(2.4, "mm"), type = "closed")) +
  annotate("text", x = 4.0, y = 1.7, label = "\"more data\" without knowledge gain",
           size = 2.9, colour = BLUE) +
  # Q4 -> Q2: knowledge held without local data
  annotate("segment", x = 3.9, xend = 3.9, y = 1.9, yend = 6.4, linewidth = 0.9, colour = GREEN,
           arrow = arrow(length = unit(2.6, "mm"), type = "closed")) +
  annotate("text", x = 3.65, y = 4.1, label = "knowledge without local data",
           size = 2.9, colour = GREEN, angle = 90) +
  # Q3 -> Q1: modelling and integration. The arrow leaves the Q1/Q3 boundary, and a crossbar
  # perpendicular to it (on paper, which is why the offsets differ in x and y) marks that base.
  annotate("segment", x = 5.0 - 0.29, xend = 5.0 + 0.29, y = 5.0 + 0.50, yend = 5.0 - 0.50,
           linewidth = 1.0, colour = GREEN) +
  annotate("segment", x = 5.0, xend = 6.9, y = 5.0, yend = 7.7, linewidth = 1.0, colour = GREEN,
           arrow = arrow(length = unit(2.8, "mm"), type = "closed")) +
  annotate("text", x = 6.35, y = 6.2, label = "modelling and integration",
           size = 3.0, colour = GREEN, fontface = "bold", angle = 42) +
  # Q1 -> Q2: transfer without revalidation
  annotate("segment", x = 9.5, xend = 3.5, y = 8.1, yend = 8.1, linewidth = 0.9, colour = RED,
           linetype = "42", arrow = arrow(length = unit(2.6, "mm"), type = "closed")) +
  annotate("text", x = 2.4, y = 8.45, label = "transfer without revalidation (no-analogue risk)",
           size = 2.8, colour = RED) +

  ## axes drawn inside the panel, so no margin is wasted
  annotate("segment", x = -0.12, xend = -0.12, y = 0, yend = 10.0, linewidth = 0.6, colour = "grey30",
           arrow = arrow(length = unit(2.4, "mm"), type = "closed")) +
  annotate("segment", x = 0, xend = 10.0, y = -0.12, yend = -0.12, linewidth = 0.6, colour = "grey30",
           arrow = arrow(length = unit(2.4, "mm"), type = "closed")) +
  annotate("text", x = -0.45, y = 5, label = "information relative to the conceptual model",
           angle = 90, size = 3.1, fontface = "bold", colour = "grey20") +
  annotate("text", x = 5, y = -0.5, label = "data availability",
           size = 3.1, fontface = "bold", colour = "grey20") +
  coord_cartesian(xlim = c(-0.75, 10.02), ylim = c(-0.85, 10.02), expand = FALSE, clip = "off") +
  theme_void(base_size = 9) + theme(plot.margin = margin(3, 3, 3, 3))

ggsave(file.path(OUT, "fig1_data_information_plane.png"), p, width = 180, height = 115, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig1_data_information_plane.pdf"), p, width = 180, height = 115, units = "mm", bg = "white")
cat("Figure 1 written to", OUT, "\n")
