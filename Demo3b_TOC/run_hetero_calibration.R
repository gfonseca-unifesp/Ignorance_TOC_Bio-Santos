# =============================================================================
# run_hetero_calibration.R — recalibrate the ignorance maps of the clean chain with heteroscedastic error models
#   (C4, C5: + depth and |latitude|; selection by stratum calibration), then refresh the outputs that read them.
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" run_hetero_calibration.R
# Predictions and withheld residuals are unchanged (same RF seeds); only expected error and its map change.
# =============================================================================
setwd(file.path(Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS"), "Demo3b_TOC"))  # F0.2
run <- function(script, iter, strict = TRUE) {
  Sys.setenv(DEMO3B_ITER = iter)
  message("\n==== ", script, " [", iter, "] ", format(Sys.time(), "%H:%M:%S"), " ====")
  if (strict) return(invisible(source(script, local = new.env())))
  tryCatch(source(script, local = new.env()), error = function(e) message("!!!! FAILED ", script, " [", iter, "]: ", conditionMessage(e)))
}
Sys.setenv(CMP_FROM = "iter1c_clean", CMP_TO = "iter2c_clean", CAL_SELECT = "strata_rule",
           CAL_SEL_ITERS = "iter1c_clean,iter2c_clean")
run("R/07_calibration.R", "iter1c_clean")
message("HETERO ITER1C DONE ", format(Sys.time(), "%H:%M:%S"))
run("R/07_calibration.R", "iter2c_clean")
message("HETERO ITER2C DONE ", format(Sys.time(), "%H:%M:%S"))
for (s in c("R/08z_calibration_selection.R", "R/08y_magnitude_by_stratum.R", "R/08_compare_iterations.R",
            "R/09_audit_figure.R", "R/09b_toc_map_iter2.R", "R/08f_localization_by_scale.R", "R/10_clean_numbers.R"))
  run(s, "iter2c_clean", strict = FALSE)
message("HETERO CHAIN DONE ", format(Sys.time(), "%H:%M:%S"))
