# Occupancy-Based Building Clustering

Replication repository for the occupancy clustering analysis presented in:

**Schaffer M, Vera-Valdés JE, Marszal-Pomianowska A (submitted) Clustering Large-Scale Residential Binary Occupancy Time Series Data, Nordic Energy Informatics Academy Conference 2026**

The script clusters buildings by their hourly occupancy patterns
using three complementary methods and compares the resulting groupings.

If you use this code, please cite the above, mentioned publication. 

## Methods overview

| Step | Method | R package | Reference |
|------|--------|-----------|-----------|
| Categorical functional data analysis | CFDA | `cfda` | Preda et al. (2021) |
| Model-based sequence clustering | ClickClust | `ClickClust` | Melnykov (2016) |
| Functional co-clustering | FunLBM | `funLBM` | Bouveyron et al. (2018) |


## Requirements

- **R** ≥ 4.2.0  
- **renv** for reproducible package management (see [renv documentation](https://rstudio.github.io/renv/))

All package versions are recorded in `renv.lock`. No manual package
installation is required beyond the steps in [Getting started](#getting-started).

## Getting started

### 1 — Clone the repository

```bash
git clone https://github.com/<owner>/<repo>.git
cd <repo>
```

### 2 — Restore the R environment

Open R in the project root and restore all packages at the exact versions
recorded in `renv.lock`:

```r
install.packages("renv")   # only needed if renv is not yet installed
renv::restore()
```

`renv` will download and install every package into an isolated project
library. This does not affect your system-wide R installation.


## Data
The data used in our publication is not included in this repository.
The data is the result of the method presented in: 

**Schaffer M, Vera-Valdés JE, Marszal-Pomianowska A (2026) Non-intrusive hourly occupancy detection in residential buildings using remotely readable water meter data: Validation and large-scale analysis. Build Environ 288:113917. https://doi.org/10.1016/j.buildenv.2025.113917**

The format of the data as expected by the code is described below. 


| Column | Type | Description |
|--------|------|-------------|
| `bldg_id` | integer | Unique building identifier |
| `time` | POSIXct (tz = `Europe/Copenhagen`) | Timestamp at hourly resolution |
| `occ_estimated` | logical / 0–1 | Binary estimated occupancy indicator |


## Output

All figures are written to `plots/` as PDF files.

| File | Contents | Paper reference |
|------|----------|-----------------|
| `01_hourly_occupancy.pdf` | Annual occupancy heatmap for all buildings | Not shown |
| `02_hourly_occupancy_cfda.pdf` | Per-cluster heatmaps — CFDA | Figure 1|
| `03_hourly_occupancy_lbm.pdf` | Per-cluster heatmaps — FunLBM |Figure 2 |
| `04_lbm_combined.pdf` | FunLBM time-cluster calendar + mean profiles | Figure 3|
| `05_cls_comparison.pdf` | Cross-method cluster comparison (bubble plots) |Figure 4 |

Intermediate clustering objects are cached in `data/cluster/` as `.RDS` files
so that the computationally expensive fitting steps (CFDA encoding, ClickClust,
FunLBM) can be skipped on subsequent runs by reading directly from disk.


## Computational notes

- The CFDA encoding step uses **16 cores** (`nCores = 16`).  
- The ClickClust and FunLBM grid searches use **8 cores** via `furrr` /
  `future`.  
- The FunLBM basis selection step uses **5 cores**.  

Adjust the `workers` arguments and `nCores` to match your hardware before
running.


## References

Bouveyron C, Bozzi L, Jacques J, Jollois FX (2018) The functional latent block model for the co-clustering of electricity consumption curves. J R Stat Soc Ser C Appl Stat 67:897–915. https://doi.org/10.1111/rssc.12260
  
Melnykov V (2016) ClickClust: An R package for model-based clustering of categorical sequences. J Stat Softw 74:. https://doi.org/10.18637/jss.v074.i09
  
Preda C, Grimonprez Q, Vandewalle V (2021) Categorical functional data analysis. The cfda r package. Mathematics 9:1–31. https://doi.org/10.3390/math9233074
  
