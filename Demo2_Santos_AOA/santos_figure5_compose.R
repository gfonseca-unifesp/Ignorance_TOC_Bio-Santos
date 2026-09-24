## ================================================================
## santos_figure5_compose.R — MSv9: Figure 5 of the manuscript, from two panels
## ================================================================
## The pre-submission review (minor point 20) found the three-panel Figure 5 of MSv8 too dense, and
## its ladder panel inconsistent with the "ordered descent" language (point M2). The ladder becomes a
## table in the main text, and the figure keeps the two panels that carry data:
##   a  santos_ignorance_aoa_map_v2.png   — santos_ignorance_aoa_map_v2.R (propagated applicability)
##   b  fig_santos_sampling_priority.png  — santos_sampling_priority.R (where and when to sample next)
## Output: outputs/fig5_santos_composite.png; the MSv8 composite is kept.
##   "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" Demo2_Santos_AOA/santos_figure5_compose.R
suppressPackageStartupMessages(library(magick))
ROOT <- Sys.getenv("IGNORANCE_ROOT", "C:/Users/Gustavo/OneDrive/Documentos/Ignorance_MS")
OUT  <- file.path(ROOT, "Demo2_Santos_AOA", "outputs")

TAG <- Sys.getenv("SANTOS_FIG_TAG", "")     # MSv10 (F3.1): usa o painel b regerado, sem sobrescrever
a <- image_read(file.path(OUT, "santos_ignorance_aoa_map_v2.png"))
b <- image_read(file.path(OUT, paste0("fig_santos_sampling_priority", TAG, ".png")))
label <- function(img, txt) image_annotate(img, txt, size = round(image_info(img)$width / 40),
                                           weight = 700, gravity = "northwest", location = "+20+10", color = "black")
# MSv10 (F3.1): larguras 1.3 : 1, para o painel b deixar de ficar espremido
W <- 3600; wa <- round(W * 1.3 / 2.3); wb <- W - wa
A <- label(image_scale(a, as.character(wa)), "a")
Bp <- label(image_scale(b, as.character(wb)), "b")
H <- max(image_info(A)$height, image_info(Bp)$height)
pad <- function(im) image_extent(im, paste0(image_info(im)$width, "x", H), gravity = "north", color = "white")
fig <- image_append(c(pad(A), pad(Bp)))
fig <- image_background(fig, "white")
image_write(fig, file.path(OUT, paste0("fig5_santos_composite", TAG, ".png")), format = "png", density = 300)
cat("Figure 5 written:", paste(image_info(fig)[, c("width", "height")], collapse = " x "), "px\n")
