# ADRIA-reef-indicators

Time series analyses for GBR resilience using ADRIA - CoralBlox reef ecosystem modelling suite.
- Clustering reef timeseries to identify what characteristics define resilient clusters.
- Assessing the characteristics that influence positive carbonate budget maintainence into the future.

## Structure

Repository structure:
``` code
ADRIA-reef-indicators/
├─ src/     # Analysis scripts
├─ outputs/     # ADRIA results sets and analysis / figure outputs
├─ data/
├─ paper/   # Reproducible quarto manuscript via the rendering of .qmd file
├─ .gitignore
├─ config.toml  # configuration file for project
├─ Project.toml     # Julia project spec
├─ Manifest.toml    # Julia project spec
└─ README.md    # this file
```

The ADRIADomain `GBR_2024_10_15_HighResCoralStress` is required to run ADRIA in this project. The path to this domain folder should be specified in the config.toml file.

The following files should be located in the data folder:
``` code
ADRIA-reef-indicators/data/
├─ GBRMPA_Management_Areas.gpkg    # File containing polygons of GBRMPA management areas
├─ GBRMPA_Reef_Features.gpkg    # File containing polygons of GBRMPA reef feature data
└─ GBRMPA_Reefal_Bioregions.gpkg   # File containing polygons of GBRMPA reefal bioregions.
```

The following structure and files should be maintained in the outputs folder:
``` code
ADRIA-reef-indicators/outputs/
└───ADRIA_results
    └───HighResCoralStress
        ├───figs    # Directory directly containing general GCM-wide and paper-methods figures, as well as GCM subdirectories
        │   ├───ACCESS-CM2
        │   │   ├───bioregion
        │   │   ├───gbr
        │   │   └───management_area
        │   ├───ACCESS-ESM1-5
        │   │   ├───bioregion
        │   │   ├───gbr
        │   │   └───management_area
        │   ├───EC-Earth3-Veg
        │   │   ├───bioregion
        │   │   ├───gbr
        │   │   └───management_area
        │   ├───GFDL-CM4
        │   │   ├───bioregion
        │   │   ├───gbr
        │   │   └───management_area
        │   └───NorESM2-MM
        │       ├───bioregion
        │       ├───gbr
        │       └───management_area
        ├───processed_model_outputs     # Folder containing scenario-median timeseries arrays for each GCM
        │       ├───median_cover_ACCESS-CM2.nc
        │       ├───median_cover_ACCESS-ESM1-5.nc
        │       ├───median_cover_EC-Earth3-Veg.nc
        │       ├─── median_cover_GFDL-CM4.nc
        │       └───median_cover_NorESM2-MM.nc
        └───analysis_context_layers_carbonate.gpkg      # Geopackage containing the required timeseries clustering, connectivity and carbonate budget values for each reef
```

### Analysis scripts - `src/ADRIA/HighResCoralStress/`

Note that scripts 1, 2 and 3 are only required to be run if creating model outputs from scratch. 
If median cover timeseries and `analysis_context_layers_carbonate.gpkg` are already available in the relevant output structure, then scripts 4 and 5 can be run immediately.
If figures are already available, then the paper.qmd document can be rendered immediately.

Workflow:
```{mermaid}
flowchart LR
    A[Start] --> B{All publication figures produced?}
    B -->|Yes| C[Render paper .qmd file]
    B -->|No| D{median_cover_GCM*.nc timeseries and analysis_context_layers_carbonate.gpkg available?}
    D -->|Yes| E[Run scripts 4 and 5]
    E --> C
    D -->|No| F[Run scripts 1 and 2]
    F --> E
    E --> C
```

- `1_GCM_model_runs.jl` : Script to run ADRIA - CoralBlox with the required ADRIADomain and produce processed median coral cover outputs.
- `2_timeseries_clustering_absolute.jl` : Perform timeseries clustering analysis on reef-area-scaled absolute coral cover. Analysis is performed at bioregion, management area and GBR wide scales. Cluster assignments for each reef are saved as `analysis_context_layers_carbonate.jl`
- `3_collating_context_layers.jl` : Attach the required values for each reef for analysis factors including connectivity metrics, reef DHW levels and perform carbonate budget analysis for a range of carbonate budget live coral cover thresholds.
- `4_analysis_plots_clustering_carbonate.jl` : Create required results for paper using timeseries data and analysis factors collated in `3_*.jl`. Additionally, calculate the proportion of reefs that change cluster assignment across the GCM levels and the proportion of reefs with depths of 1-10m that occur in low and medium coral cover clusters.
- `5_publication_plots.jl` : Create extra required plots for publication such as a context map and diagram of GCM Equilibium Climate Sensitivity values.

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