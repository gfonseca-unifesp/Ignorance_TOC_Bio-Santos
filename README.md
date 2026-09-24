# Ignorance_TOC_Bio-Santos

Code and data for *Measuring and reducing ignorance in data- and AI-rich environmental sciences*
(G. Fonseca, Marine Science Institute, Federal University of São Paulo — manuscript submitted to
*Ecological Monographs*).

The repository holds what is needed to reproduce the results: the analysis code of both case studies
and the derived data package.

- **Case 1, global:** a re-analysis of a published machine-learning map of seafloor total organic carbon
  (TOC) — spatial validation, dissimilarity index and area of applicability, an ignorance map calibrated
  against realised error, a region-by-region decision rule, and where the next samples should be taken.
- **Case 2, regional:** an audit of a hybrid structural random-forest model of deep-sea benthos in the
  Santos Basin (Brazil) — applicability propagated through the model hierarchy, and where and when to
  sample next.

## Contents

| Folder | Contents |
|---|---|
| `data_release/` | **The data package.** Training cells, the global grid of predictions, dissimilarity, applicability, calibrated expected error and distance to data (GeoTIFF and table), the nested calibration points, the environmental-space sample, the ranked sampling sites of both cases, and the Santos Basin station and grid layers. `data_release/README.md` is the variable dictionary (units, CRS, provenance, licence); `verify_release.R` checks every file |
| `Demo3b_TOC/` | Case 1 code. `R/` holds the pipeline in order: data and quality control (`01*`), predictors (`02*`), model and spatial cross-validation (`03*`), prediction and applicability (`04`), transfer (`05`), calibration of the ignorance map (`07*`), the analyses behind each result (`08*`), figures (`06`, `09*`), all quoted numbers (`10_clean_numbers.R`) and the data package (`11_export_data.R`). `run_all.R` and `run_clean_chain.R` run the chain; `data/toc_cells.csv` is the aggregated input (12,944 cells of 0.1°) |
| `Demo2_Santos_AOA/` | Case 2 code: applicability of the published models (`santos_aoa.R`), the ignorance × applicability map (`santos_ignorance_aoa_map_v2.R`), the audit (`santos_audit*.R`), the structure of the error, the sampling priority (`santos_sampling_priority.R`), when to sample (`santos_when.R`), the intervals of the where-or-when test (`santos_T4b_intervals*.R`), Figure 5 (`santos_ignorance_ladder.R`, `santos_figure4_compose.R`) and the data package (`santos_export_data.R`) |
| `figures_src/` | The script that draws Figure 1 of the paper, the data-information plane |
| `Demo3a_AOA/R/` | The two scripts that build the shared predictor cache used by Case 1 |
| `sessionInfo_MSv6.txt` | R and package versions of the last full run |

## Reproducing

R ≥ 4.5 with `CAST (>= 1.1.1)`, `caret`, `ranger`, `randomForest`, `terra`, `sf`, `FNN`, `scam`, `ggplot2`,
`patchwork`, `tidyterra`, `rnaturalearth`, `ggnewscale`, `magick`, `dplyr` and `sessioninfo`.

Copy `.Renviron.example` to `.Renviron` in the repository root **and** in each case folder (R reads it from
the working directory) and set:

```
IGNORANCE_ROOT   the repository root
DEMO3A_CACHE     the predictor cache (Bio-ORACLE v3 layers and the NN-TOC label and texture files)
```

Build the cache with `Demo3a_AOA/R/02_predictors.R` and `Demo3b_TOC/R/02*.R`, then run
`Demo3b_TOC/run_all.R` from `Demo3b_TOC/`. For Case 2, run the scripts in `Demo2_Santos_AOA/` from the
repository root. `data_release/verify_release.R` checks the data package.

## Not included, and where to get it

- **Raw TOC labels:** the NN-TOC compilation, Zenodo [10.5281/zenodo.11186224](https://doi.org/10.5281/zenodo.11186224)
  (CC-BY 4.0; Parameswaran et al. 2025, *Geoscientific Model Development* 18:2521–2544). Cite that DOI.
- **Predictor layers:** Bio-ORACLE v3 (Assis et al. 2024) and the NN-TOC texture and supply layers; the cache
  (~4 GB) is rebuilt by the scripts above.
- **Santos Basin models and biological data:** Fonseca et al. (2026), distributed through iMESC; the scripts
  read them from the iMESC savepoints.

## Licence and citation

Code under the MIT licence (`LICENSE`). The data package inherits CC-BY 4.0 from the NN-TOC labels.
A versioned archive of this repository with a DOI is on Zenodo: [DOI to be added on release].
