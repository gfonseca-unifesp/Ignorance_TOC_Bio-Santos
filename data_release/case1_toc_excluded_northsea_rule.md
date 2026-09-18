# Case 1 — the North Sea exclusion rule

The NN-TOC label compilation (Zenodo 10.5281/zenodo.11186224, CC-BY 4.0) holds 110,149 records of total
organic carbon in surface sediments. One regional set inside the North Sea box is a **gridded product**,
not a set of measurements, and is attributed to a personal communication:

- coordinates on a regular 0.0298 deg grid; 98.5% of coordinates with 5 decimals;
- 89% of TOC values with 5 or more decimals;
- within-cell standard deviation 0.034 log10 units, against 0.12 elsewhere.

**Rule applied here:** drop records inside the North Sea box whose TOC value carries 4 or more decimal
places. It removes **84,265 records** and keeps **1,743** conventionally reported samples from that box.
Keeping the product would train the model on another model's output and give the North Sea 40% of all cells.

After the rule, 25,884 measurements aggregate into **12,944 cells** of 0.1 deg (mean of log10 TOC per cell).

Source: Demo3b_TOC/R/01_data.R and Demo3b_TOC/README.md.
