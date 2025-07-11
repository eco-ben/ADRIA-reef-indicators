# ADRIA-reef-indicators

Time series analyses for GBR health using ReefModEngine ecosystem model.
- Looking for reefs that decline early compared to their surrounding reefs. These could be considered as bellwether reefs that indicate oncoming decline.
- Looking for reefs that have a delay in decline compared to their surrounding reefs. These could be considered as more resilient to the surrounding stressors.

Previously this study used ADRIA and CoralBlox for ecological modelling. This has moved to focus on
ReefModEngine instead, scripts used for ADRIA modelling have been archived to `src/ADRIA` folder.

## Structure

``` code
ADRIA-reef-indicators/
├─ src/           # Analysis scripts
├─ outputs/
├─ figs/
├─ data/
├─ .gitignore
├─ config.toml    # configuration file for ADRIA
├─ Project.toml   # Julia project spec
└─ README.md      # this file
```

### Analysis scripts - `src/RME`

- `rme_1_1_initial_runs_remote_desktop.jl` : Script for running RME scenarios on the remote desktop (takes 3 days to run 1200 scenarios).
- `rme_1_2_processing_run_results.jl` : Script for taking the RME results for each GCM and combining them into a single ensemble dataset, then calculating median cover and taxa evenness.
- `rme_2_Subregion_correlation.jl` : Lagged correlation analysis at for bioregions.
- `rme_3_analysis_context_layers.jl` : Addition of context layers to the correlation results for statistical analysis.
- `rme_4_exploratory_analysis.jl` : Exploratory analysis of bellwether reefs.

### Data Requirements

- `raw scenario data` : Raw scenario data resulting from `rme_1_1_initial_runs_remote_desktop.jl` can be found at
- `processed median data` : Processed data resulting from `rme_1_2_processing_run_results.jl` can be found at https://aimsgovau.sharepoint.com/:f:/s/ADRIADevelopment/EqJeFYVKvadKop2h2sVoUpEBDqlEzVdUDq66hK8KqNvPiw?e=ptVG3c
- `canonical-reefs gpkg` : Canonical reefs geopackage can be found at _ or created by running _
- `GBRMPA_Reefal_bioregions.gpkg` : https://aimsgovau.sharepoint.com/:u:/s/ADRIADevelopment/EefReOCyAwVCnhKHLdTPeuYBIiNfoTucyMcV4MgtU9PVbA?e=Yso4uC

## Methods

### Domain

Using GBR-wide domain from ReefModEngine data: "rme_ml_2024_06_13".

Currently using RCP 45 from domain.

#### Includes:

- Initial coral cover data from RME
- DHW/environmental data from RME
- Reef spatial data from RME
- Connectivity data from RME

### Scales

Timeseries are grouped by reefs into their bioregion (GBRMPA reefal bioregion data).
The timeseries are then compared for each reef within these groups with the overall group median timeseries.

### Lagged Correlation Analyses

Lagged correlation analysis is performed between each reef in a group and the group's median time series.
A strong-positive correlation at positive lags means a reef's changes are early compared to group-median. A strong-positive correlation at negative lags means a reef's changes are delayed compared to the group-median.

## Results

Results from lagged correlation analyses are used to explore characteristics of reefs that occur early/late compared to their surrounding group.
Reefs with timeseries that occur early/late at both subregion/bioregion scales are considered to be candidates for further analyses.

These candidate reefs are then explored for their common characteristics that differ from other surrounding reefs, such as:

- Strength/number of outgoing and incoming larval connections
- Connectivity centrality (eigenvector)
- Initial coral cover values
- Mean DHW thermal stress experienced
- Correlation of cover/evenness timeseries with DHW timeseries (susceptibility)
- Crown of Thorns Starfish Mortality