# =============================================================================
# run_clean_chain.R — the clean chain (provenance check, ANALYTIC_LOG): steps 03 onward for
#   iter1c_clean (base model with supply; distance to large land masses; bottom layers smoothed in deep water)
#   iter2c_clean (+ sample-based lithology)
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" run_clean_chain.R
# Core steps stop the chain on error; later analyses are logged and skipped on error.
# =============================================================================
setwd(file.path(Sys.getenv("IGNORANCE_ROOT", "C:/Users/fonse/OneDrive/Documentos/Ignorance_MS"), "Demo3b_TOC"))  # F0.2
run <- function(script, iter, strict = TRUE) {
  Sys.setenv(DEMO3B_ITER = iter)
  message("\n==== ", script, " [", iter, "] ", format(Sys.time(), "%H:%M:%S"), " ====")
  if (strict) return(invisible(source(script, local = new.env())))
  tryCatch(source(script, local = new.env()), error = function(e) message("!!!! FAILED ", script, " [", iter, "]: ", conditionMessage(e)))
}
Sys.setenv(CMP_FROM = "iter1c_clean", CMP_TO = "iter2c_clean")
core <- c("R/03_model.R", "R/04_predict_aoa.R", "R/05_transfer.R", "R/06_figures.R", "R/07_calibration.R")

run("R/02e_clean_stacks.R", "iter2c_clean")
for (s in core) run(s, "iter1c_clean")
message("CLEAN ITER1C CORE DONE ", format(Sys.time(), "%H:%M:%S"))
for (s in core) run(s, "iter2c_clean")
message("CLEAN ITER2C CORE DONE ", format(Sys.time(), "%H:%M:%S"))

run("R/07b_diagnose.R", "iter1c_clean", strict = FALSE)        # diagnosis H1-H4 of the clean base model
for (s in c("R/08_compare_iterations.R", "R/08k_iteration_bootstrap.R", "R/09_audit_figure.R", "R/09b_toc_map_iter2.R",
            "R/08o_prediction_artefact.R", "R/08d_lithology_edges.R", "R/08e_edges_realised_error.R", "R/08f_localization_by_scale.R"))
  run(s, "iter2c_clean", strict = FALSE)
message("CLEAN COMPARISONS DONE ", format(Sys.time(), "%H:%M:%S"))
run("R/08g_regionalization.R", "iter1c_clean", strict = FALSE)  # decision rule by region
run("R/08h_targeted_sampling.R", "iter1c_clean", strict = FALSE) # collection strategies
run("R/08l_collection_by_region.R", "iter1c_clean", strict = FALSE)
run("R/08i_ignorance_spaces.R", "iter1c_clean", strict = FALSE)
run("R/08n_space_env_structure.R", "iter1c_clean", strict = FALSE)
run("R/08g_regionalization.R", "iter2c_clean", strict = FALSE)  # decision rule after texture (was pending)
message("CLEAN CHAIN DONE ", format(Sys.time(), "%H:%M:%S"))
