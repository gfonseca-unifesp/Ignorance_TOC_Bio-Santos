## ================================================================
## santos_figure4_compose.R — F6.4 of ROADMAP_MSv6: Figure 4 of the manuscript, from its three panels
## ================================================================
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_figure4_compose.R
## Panels (all read from outputs/):
##   a  santos_ignorance_ladder.png       — santos_ignorance_ladder.R (the MSv5 image is kept as
##      santos_ignorance_ladder_MSv5.png for provenance; L3 changed from ACHIEVED to PARTIAL)
##   b  santos_ignorance_aoa_map_v2.png   — santos_ignorance_aoa_map_v2.R
##   c  fig_santos_sampling_priority.png  — santos_sampling_priority.R
## Layout: a above b on the left, c on the right at the same height. Output: outputs/fig4_santos_composite.png
suppressPackageStartupMessages(library(magick))
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS")
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs")

a <- image_read(file.path(OUT, "santos_ignorance_ladder.png"))
b <- image_read(file.path(OUT, "santos_ignorance_aoa_map_v2.png"))
c <- image_read(file.path(OUT, "fig_santos_sampling_priority.png"))

label <- function(img, txt) image_annotate(img, txt, size = round(image_info(img)$width / 40),
                                           weight = 700, gravity = "northwest", location = "+20+10", color = "black")
W <- image_info(b)$width
a <- image_scale(a, as.character(W))
left  <- image_append(c(label(a, "a"), label(b, "b")), stack = TRUE)
H <- image_info(left)$height
right <- image_border(label(image_scale(c, paste0("x", H)), "c"), "white", "40x0")
fig <- image_background(image_append(c(left, right)), "white")
fig <- image_scale(fig, "5400")                     # ~46 cm wide at 300 dpi; journals rescale
image_write(fig, file.path(OUT, "fig4_santos_composite.png"), format = "png", density = 300)
cat("Figure 4 written:", paste(image_info(fig)[, c("width", "height")], collapse = " x "), "px\n")
