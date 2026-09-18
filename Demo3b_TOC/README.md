# Demo 3b — Global ignorance map of seafloor total organic carbon (TOC)

This demo applies the ignorance framework to a variable whose prediction error can be **measured**. The presence-background SDM of Demo 3a (*Pocillopora damicornis*) could not show that error grows with dissimilarity: AUC in presence-background data reflects class overlap and case-mix, not model knowledge. Every TOC record is a measurement, so error is observed directly.

Global TOC maps are routinely produced with machine learning (Lee et al. 2019, kNN; Atwood et al. 2020, RF; NN-TOC, GMD 2025, deep learning) and used in policy debates (e.g., the carbon impacts of bottom trawling). NN-TOC validated with a random 85:15 split. It used no spatial CV, no test of uncertainty against realised error, and no area of applicability.

## Run

```
"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" run_all.R
```

Requires the Demo 3a predictor cache (`Demo3a_AOA/R/02_predictors.R`) and the NN-TOC label files in `<cache>/toc`.

| Step | Script | What it does |
|---|---|---|
| 1 | `R/01_data.R` | Clean labels, exclude the North Sea gridded product, snap to ocean, aggregate to 0.1° cells; noise floor |
| 2 | `R/02_predictors.R` | Extra Bio-ORACLE v3 layers (chl, near-bottom phytoplankton, bottom current, slope, TRI, TPI), distance to coast, log transforms |
| 3 | `R/03_model.R` | RF regression of log10(TOC + 0.01): 1000-km spatial block CV, `ffs`, random-CV twin, geodist diagnostic |
| 4 | `R/04_predict_aoa.R` | Global prediction, DI, AOA, ignorance map (expected RMSE from `CAST::errorProfiles`) |
| 5 | `R/05_transfer.R` | Leave-region-out: error vs DI in withheld regions; calibration of the ignorance map (expected vs realised RMSE) |
| 6 | `R/06_figures.R` | Main figure |

## Data and QC decisions

Source: NN-TOC label compilation, extracted by HTTP range requests from Zenodo `10.5281/zenodo.11186224` (CC-BY 4.0; only `data/raw/labels/toc_continentalshelves.csv` and `toc_deep.csv`). It compiles Seiter et al. (2004), Romankevich et al. (2009), van der Voort et al. (2021) / MOSAIC, Paradis et al. (2023) and regional sets.

- **Units** are consistent. At 3,874 coordinates shared with Seiter et al. (2004), 94.8% of values are identical (median ratio 1.00). Different medians between the shelf and deep sets reflect sampling composition (sandy shelves vs upwelling margins), not units.
- **North Sea regional set excluded.** It is attributed to "Wenyen Zhang, personal communication, 2023, HEREON", and the records are a gridded product:
  - regular 0.0298° grid;
  - 98.5% of coordinates with 5 decimals;
  - 89% of TOC values with ≥ 5 decimals;
  - within-cell SD 0.034 log10 units (0.12 elsewhere).

  Keeping it would train on another model's output and give the North Sea 40% of all cells. The rule drops TOC values with ≥ 4 decimals inside the North Sea box (84,265 records) and keeps 1,743 conventionally reported samples.
- **Aggregation**: one observation per 0.1° cell (mean of log10 TOC). Result: 25,884 records → 12,944 cells.
- **Noise floor**: SD of log10 TOC across distinct sites (0.01°) within cells with ≥ 3 sites; median 0.123 (453 cells). This is the variability no model at 0.1° can remove (level L5 in the manuscript).

## Key references

- Meyer & Pebesma (2021) *Methods Ecol Evol* 12:1620–1633 (AOA, error profiles).
- Seiter, Hensen, Schröter & Zabel (2004) *Deep-Sea Res I* 51:2001–2026; PANGAEA doi:10.1594/PANGAEA.199835.
- Lee, Wood & Phrampus (2019) *Global Biogeochem Cycles* 33:37–46.
- Atwood et al. (2020) *Front Mar Sci* 7:165.
- NN-TOC v1: *Geosci Model Dev* 18:2521 (2025); data doi:10.5281/zenodo.11186224.
- van der Voort et al. (2021) *ESSD* 13:2135; Paradis et al. (2023) *ESSD* 15:4105 (MOSAIC).
