# =============================================================================
# Demo 3b — seafloor TOC ignorance map (R >= 4.6; run from this folder)
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" run_all.R
# Requires the Demo 3a predictor cache (Demo3a_AOA/R/02_predictors.R) and the
# NN-TOC label and texture files in <cache>/toc (see 01_data.R and 02b_texture.R).
# The analysis is iterative (see ANALYTIC_LOG.md); each iteration writes to outputs/<iteration>/.
# =============================================================================
run <- function(script, iter) {
  Sys.setenv(DEMO3B_ITER = iter)
  message("\n==== ", script, " [", iter, "] ====")
  source(script, local = new.env())
}
shared <- c("R/01_data.R", "R/01b_nntoc_composition.R", "R/02_predictors.R", "R/02b_texture.R", "R/02c_supply.R")   # data and predictor stacks
per_iter <- c("R/03_model.R",        # spatial CV, ffs, random forest regression
              "R/04_predict_aoa.R",  # prediction, DI, AOA, ignorance map
              "R/05_transfer.R",     # leave-region-out: error vs DI
              "R/06_figures.R",      # main figure of the iteration
              "R/07_calibration.R")  # nested calibration of the ignorance map

for (s in shared) run(s, "iter1_baseline")
for (s in per_iter) run(s, "iter1_baseline")          # iteration 1: baseline predictors
run("R/07b_diagnose.R", "iter1_baseline")             # why the map failed to locate error (before iteration 2)
for (s in per_iter) run(s, "iter1b_supply")           # iteration 1b: + carbon supply to the seafloor
run("R/03b_texture_screen.R", "iter2_texture")        # paired test of H3: does texture reduce error (spatial and region CV)?
for (s in per_iter) run(s, "iter2_texture")           # iteration 2: + sediment texture
run("R/08_compare_iterations.R", "iter2_texture")     # knowledge gain between iterations
run("R/08b_ridge_check.R", "iter2_texture")           # where expected error rose: ridges / relief classes
run("R/08c_ridge_mechanism.R", "iter2_texture")       # roughness deciles, texture outside sampled range, zooms
run("R/08d_lithology_edges.R", "iter2_texture")       # ridge-axis lithology class and polygon boundaries
run("R/08e_edges_realised_error.R", "iter2_texture")  # realised vs expected error change at those features
run("R/08f_localization_by_scale.R", "iter2_texture") # does the ignorance map locate error at 1-10 deg and region scale?
for (it in c("iter1_baseline", "iter2_texture")) {
  run("R/07c_learned_ignorance.R", it)                # H6: learned error model vs DI-based maps
  run("R/08g_regionalization.R", it)                  # decision rule: collect data, regionalize, or new drivers
}
run("R/08h_targeted_sampling.R", "iter1_baseline")    # H7: directed vs random addition of data within regions
run("R/08h_targeted_sampling.R", "iter1b_supply")     # H7 repeated with the supply predictors
run("R/08i_ignorance_spaces.R", "iter1_baseline")     # the same ignorance in geographic and environmental space
run("R/08j_arctic_check.R", "iter1b_supply")          # why supply layers worsened transfer into the Siberian Arctic
run("R/08l_collection_by_region.R", "iter1b_supply")  # per region: path (08g) + how to collect (08h)
run("R/01c_nntoc_qc_emulation.R", "iter1b_supply")    # NN-TOC composition by emulating its published QC
# presentation sequence: base model iter1b (supply) -> iter2b (texture on top of supply)
Sys.setenv(CMP_FROM = "iter1b_supply", CMP_TO = "iter2b_supply_texture")
run("R/07b_diagnose.R", "iter1b_supply")
run("R/02d_supply_texture.R", "iter2b_supply_texture")
for (s in per_iter) run(s, "iter2b_supply_texture")
for (s in c("R/08_compare_iterations.R", "R/08d_lithology_edges.R", "R/08e_edges_realised_error.R",
            "R/08f_localization_by_scale.R", "R/08k_iteration_bootstrap.R", "R/09_audit_figure.R")) run(s, "iter2b_supply_texture")
run("R/09_audit_figure.R", "iter2_texture")           # conceptual figure: audit loop + ignorance budget
# provenance check of the integrated layers (08q-08u), then the clean chain used by the report and manuscript:
for (s in c("R/08q_layer_sensitivity.R", "R/08r_wedge_attribution.R", "R/08s_wedge_neutralise.R", "R/08t_distance_to_land.R", "R/08u_clean_candidate.R"))
  run(s, "iter2b_supply_texture")
source("run_clean_chain.R")                           # iter1c_clean -> iter2c_clean, steps 02e and 03 onward
