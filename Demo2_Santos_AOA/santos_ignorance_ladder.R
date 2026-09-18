## ================================================================
## santos_ignorance_ladder.R — panel (a) of Figure 4: the Santos Basin case on the ignorance ladder
## ================================================================
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_ignorance_ladder.R
## Rebuilds, as code, the panel that the MSv5 package held only as an image
## (outputs/santos_ignorance_ladder_MSv5.png, kept for provenance). Same content and layout, with one
## change decided by the author on 18 Sep 2026:
##   L3 (sampling gaps) ACHIEVED -> PARTIAL. The stratified design covers 25-2400 m in two surveys, but
##   18.6% of the basin (mostly deeper than the design) lies in environmental extrapolation, and a new
##   survey still reduces ignorance about twice as much as new stations (santos_sampling_priority.R,
##   santos_audit_T4.R).
## Output: outputs/santos_ignorance_ladder.png (1600 x 800 px, as the original)
suppressPackageStartupMessages(library(grid))
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs")

rows <- data.frame(
  level  = c("L0 \u00b7 Meta-ignorance", "L1 \u00b7 Specified ignorance", "L2 \u00b7 Structural / integration",
             "L3 \u00b7 Sampling gaps", "L4 \u00b7 Model uncertainty", "L5 \u00b7 Noise floor"),
  sub    = c("a model exists", "named unknowns", "dispersed knowledge", "unmeasured here/now",
             "quantified error", "irreducible only"),
  did    = c("Conceptual model in place: drivers \u2192 environment \u2192 biota",
             "44 environmental + 12 biological essential variables named, each with a question",
             "Datasets integrated (water, sediment, organic, biota); little local knowledge to add",
             "198 samples \u00b7 25\u20132400 m \u00b7 2019 & 2021; deep basin and new intervals unsampled",
             "Quantified per variable \u2014 environment R\u00b2 0.10\u20130.98; biota R\u00b2 0.38\u20130.84 (macro > meio)",
             "Reached for physics (R\u00b2>0.9) & best macrofauna; meiofaunal richness still reducible"),
  status = c("ACHIEVED", "ACHIEVED", "ACHIEVED", "PARTIAL", "ACHIEVED", "PARTIAL"),
  accent = c("#1F4E5F", "#2E6478", "#4F8194", "#8FB1BC", "#BDD4D9", "#12303A"),
  stringsAsFactors = FALSE)
PILL <- c(ACHIEVED = "#2E8B6A", PARTIAL = "#C77F32")
DARK <- "#12303A"; TEAL <- "#1F4E5F"

png(file.path(OUT, "santos_ignorance_ladder.png"), width = 1600, height = 800, res = 100, type = "cairo")
grid.newpage()
pushViewport(viewport(xscale = c(0, 1600), yscale = c(0, 800)))
# layout is written in pixels measured from the top, as on the original image; Y() flips to grid's upward axis
Y <- function(y) 800 - y
R <- function(x0, y0, x1, y1, fill, col = NA, r = 0) {          # any two opposite corners
  cx <- (x0 + x1) / 2; cy <- Y((y0 + y1) / 2); w <- abs(x1 - x0); h <- abs(y1 - y0)
  if (r > 0) grid.roundrect(unit(cx, "native"), unit(cy, "native"), unit(w, "native"), unit(h, "native"),
                            r = unit(r * 0.72, "points"), gp = gpar(fill = fill, col = col))   # r in px at 100 dpi
  else grid.rect(unit(cx, "native"), unit(cy, "native"), unit(w, "native"), unit(h, "native"),
                 gp = gpar(fill = fill, col = col))
}
Tx <- function(label, x, y, size, face = "plain", col = "#1B1B1B", hjust = 0)
  grid.text(label, unit(x, "native"), unit(Y(y), "native"), hjust = hjust, vjust = 0.5,
            gp = gpar(fontsize = size, fontface = face, col = col))

Tx("The Santos Basin case on the ignorance ladder \u2014 levels reached", 70, 35, 18.5, "bold", TEAL)
R(70, 116, 440, 72, DARK);  Tx("Level", 86, 94, 14, "bold", "white")
R(452, 116, 1152, 72, DARK); Tx("What the Santos case did (R\u00b2 as the ignorance metric)", 468, 94, 14, "bold", "white")
R(1200, 116, 1530, 72, DARK); Tx("Status", 1216, 94, 14, "bold", "white")

for (i in seq_len(nrow(rows))) {
  top <- 126 + (i - 1) * 86; bot <- top + 82; mid <- (top + bot) / 2
  R(70, bot, 1530, top, if (i %% 2) "#F3F7F8" else "white", "#E3E9EB")
  R(70, bot, 82, top, rows$accent[i])
  Tx(rows$level[i], 96, mid - 13, 15, "bold", "#12303A")
  Tx(rows$sub[i], 96, mid + 17, 11, col = "#5B6B70")
  Tx(rows$did[i], 468, mid, 12.5, col = "#1B1B1B")
  R(1230, mid + 21, 1500, mid - 21, PILL[[rows$status[i]]], NA, r = 21)
  Tx(rows$status[i], 1365, mid, 13, "bold", "white", hjust = 0.5)
}
R(70, 764, 1530, 668, "#E8F0F2", "#5C7F8A", r = 8)
Tx("Ignorance concentrates at the environmental tier (R\u00b2 down to 0.10) but propagates little upward:",
   800, 697, 13.5, "bold", TEAL, hjust = 0.5)
Tx("biological predictions are comparable (mean R\u00b2 0.66) and less extreme (floor 0.38). Macrofauna (mean R\u00b2 0.78) > meiofauna (0.57).",
   800, 726, 12.5, col = "#1B1B1B", hjust = 0.5)
invisible(dev.off())
cat("ladder written:", file.path(OUT, "santos_ignorance_ladder.png"), "\n")
