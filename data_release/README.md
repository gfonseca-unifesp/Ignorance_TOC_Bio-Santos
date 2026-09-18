# Data release — *Measuring and reducing ignorance in data- and AI-rich environmental sciences*

G. Fonseca, Marine Science Institute, Federal University of São Paulo (UNIFESP) · gfonseca@unifesp.br

This package holds the derived data behind the two case studies: the cells a model was trained on, the
gridded predictions and ignorance layers, the environmental-space sample, the ranked sampling sites, and
the realised-versus-expected errors that calibrate the ignorance map. It is built by
`Demo3b_TOC/R/11_export_data.R` and `Demo2_Santos_AOA/santos_export_data.R`, and checked by
`verify_release.R` in this folder.

**Coordinate reference system:** EPSG:4326 (WGS 84). Case 1 grids are 0.1°; Case 2 grids are 0.02°.
**Response units:** Case 1, TOC in weight % of dry sediment; the model works on `log10(TOC % + 0.01)`,
so errors and expected errors are in log₁₀ units (0.30 ≈ a factor of 2).

---

## Case 1 — global seafloor organic carbon

| File | Rows / layers | Contents |
|---|---|---|
| `case1_toc_training_cells.csv` | 12,944 | the cells the model was trained on: `cell` (raster cell id), `x_lon`, `y_lat`, `n_records` (measurements aggregated), `n_sites` (distinct 0.01° sites), `y` = log10(TOC % + 0.01), `toc_median` (%), `y_sd` (within-cell SD of `y`, the noise floor), `mask`, the 30 predictor layers of the clean stack, `block` (1,000-km spatial CV block, Equal Earth), `fold` (1–5), `region` (the 8 k-means regions) |
| `case1_grid_iter2c.tif` | 7 layers, 3600 × 1800 | `toc_pct`, `toc_log10`, `DI`, `DI_over_threshold` (1 = edge of the area of applicability), `inside_AOA` (1/0), `expected_RMSE_calibrated_log10`, `distance_nearest_observation_km` |
| `case1_grid_iter2c.csv.gz` | 4,294,424 | the same seven layers per ocean cell, with `cell`, `lon`, `lat` and `region` |
| `case1_calibration_outer_points.csv` | 12,944 | nested leave-one-region-out predictions: `region`, `lon`, `lat`, `obs`, `pred`, `e2` (squared error), `DI`, `AOA`, `dist_km`, `depth`, and the expected squared error of each candidate error model (`exp_C0`–`exp_C5`; the map uses C3) |
| `case1_fig_S10_environmental_space.csv.gz` | 60,000 | area-weighted ocean sample in the dissimilarity space: `PC1`, `PC2`, `DInorm`, `lon`, `lat` |
| `case1_fig_S10_environmental_space_loadings.csv` | 23 | the PCA loadings and the dissimilarity weight of each predictor |
| `case1_sampling_priority_sites.csv` | 50 | ranked places for the next samples (representative novelty), with depth, predicted TOC, dissimilarity before and at selection, gain, expected error, distance to the nearest observation, region, decision path and recommended strategy, and two flags |
| `case1_sampling_priority_uncertainty_sites.csv` | 50 | the same ranking under the model-uncertainty criterion (sensitivity) |
| `case1_sampling_priority_curve.csv` | 51 | mean dissimilarity of the ocean against the number of sites added, for novelty, uncertainty and 20 random sets |
| `case1_toc_excluded_northsea_rule.md` | — | the exclusion rule and its counts |

**Regions** are k-means clusters of the *observation* cells. An ocean cell far from every observation is
attached to the nearest centroid and belongs to that region **nominally only**; the sampling-priority
table flags those sites (`flag_outside_region_support`, more than 1,000 km from any observation).

## Case 2 — deep-sea benthos of the Santos Basin

| File | Rows / layers | Contents |
|---|---|---|
| `case2_station_ignorance.csv` | 198 | 99 stations × 2 surveys (2019, 2021): `station`, `lon`, `lat`, `depth_m`, `campaign` (1/2), `survey`, `set` (160 training / 38 held-out test), `ignorance` (mean standardised out-of-sample \|observed − predicted\| over the 12 indicators), `bio_DI_over_threshold`, `bio_models_inside_AOA` (0–12) |
| `case2_grid_layers.tif` | 8 layers, 0.02° | for each survey: `bio_DInorm` (biological tier), `env_DInorm` (environmental tier, mean over the 44 models), `env_edge_class` (1 inside, 2 transition, 3 extrapolation) and its sieved version used for display |
| `case2_sampling_priority_A_stations.csv` | 30 | layer A: stations drawn at random within the stratified design (10 per depth zone), ≥ 20 km apart |
| `case2_sampling_priority_B_sites.csv` | 20 | layer B: ranked sites that shrink the extrapolated area, with the area remaining after each site |
| `case2_sampling_priority_curve.csv` | 21 | % of the basin in extrapolation against the number of layer-B sites, ranked and random |
| `case2_sampling_priority_summary.csv` | 11 | the numbers quoted in the text, including the reproduction check |

**Read layers A and B apart.** Layer A is the strategy the audit supports for reducing *error* (no guide
beat random stations at this level of ignorance). Layer B pursues a different objective — extending the
applicability of the published models — which was **not** tested against error.

---

## Provenance, licence and what not to redistribute

- **Case 1 labels.** Aggregated from the NN-TOC label compilation (Zenodo
  [10.5281/zenodo.11186224](https://doi.org/10.5281/zenodo.11186224), CC-BY 4.0), which itself compiles
  Seiter et al. (2004; PANGAEA 10.1594/PANGAEA.199835), MOSAIC (van der Voort et al. 2021; Paradis et al.
  2023), Romankevich et al. (2009) and regional sets. **The raw NN-TOC records are not redistributed
  here**: this package holds only the 0.1° aggregates used by the model. Cite the NN-TOC DOI for the labels.
- **Case 1 predictors.** Bio-ORACLE v3 (Assis et al. 2024) and the NN-TOC texture and supply layers; the
  stacks themselves are not redistributed.
- **Case 2.** Models and biological data are those of Fonseca et al. (2026), distributed through iMESC;
  only the derived ignorance and applicability layers are released here.
- **Licence of this package:** CC-BY 4.0, inherited from the NN-TOC labels.

## Checking the package

```
Rscript verify_release.R
```

It reads every file and checks the row and layer counts stated above.
