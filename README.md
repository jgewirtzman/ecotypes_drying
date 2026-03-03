# Ecotypes Drying Experiment

Analysis code and data for: **"Soil drying induces widespread productivity loss but unequal climate vulnerability among ecotypes of a foundational Arctic sedge"** (Gewirtzman et al.).

## Project Overview

This study examines how *Eriophorum vaginatum* tussocks from three populations along a latitudinal gradient in northern Alaska (Sagwon, Toolik Lake, Coldfoot) respond to experimental water table manipulation (Wet, Dry, Deep) in a common garden at Toolik Field Station.

### Measurements

- Soil moisture dynamics
- Plant phenology and growth
- Canopy structure (NDVI, LAI)
- Leaf gas exchange (Amax, stomatal conductance)
- A/Ci curve parameters (Vcmax, Jmax)
- Pre-dawn water potential
- Carbon isotopes (d13C) and water use efficiency
- Productivity relationships

## Repository Structure

```
R/                                  # Analysis scripts (run in numerical order)
  00_setup.R                        # Packages, constants, themes
  01_soil_moisture.R                # Soil moisture calibration & dynamics
  02_phenology.R                    # Tiller phenology mixed models & random forest
  03_canopy.R                       # LAI/NDVI analysis
  04_gas_exchange.R                 # Amax, gs, iWUE mixed models
  05_water_potential.R              # Pre-dawn water potential
  06_aci_curves.R                   # A/Ci curve analysis (uses pre-computed fits)
  06b_aci_refit.R                   # A/Ci refit: Tcorrect sensitivity analysis
  07_isotopes.R                     # d13C analysis
  08_combined_physiology.R          # Combined physiology figures
  09_productivity.R                 # NDVI-productivity relationships
  10_model_diagnostics.R            # Model diagnostic plots
  11_wue_new.R                      # WUE analysis (revised)
  12_wue_new_2.R                    # WUE analysis (alternative)
  13_senescence_dates.R             # Senescence date extraction
  14_productivity2.R                # Productivity (LAI-based)
  ecotypes_map.R                    # Study site map
  treatment_timeline.R              # Experimental timeline figure
  main_workflow.R                   # Runs all scripts in sequence
  utils/
    themes.R                        # ggplot themes
    stats.R                         # Statistical helper functions
    plotting.R                      # Plotting helper functions

data/
  raw/                              # Original data files
    gwc_calibration.csv             # Gravimetric water content calibration
    Ecotypes2017_Drying_*.csv       # Soil moisture, phenology, NDVI, Amax, etc.
    A_Ci Ouputs.csv                 # Original A/Ci curve fits
    A_Ci Ouputs_Corrected.csv       # A/Ci fits with metadata
    Ecotypes2017_Drying_A-Ci.csv    # Raw LI-6400 A/Ci curve data

output/
  figures/                          # Generated PDF/PNG figures
  tables/                           # Model summaries, CSV tables
```

## Usage

1. Open `ecotypes_drying.Rproj` in RStudio
2. Run `source("R/00_setup.R")` to install/load all packages
3. Run scripts in numerical order (01-14), or use `R/main_workflow.R`

Script `06b_aci_refit.R` is a standalone A/Ci sensitivity analysis and can be run independently after `00_setup.R`.

## Key R Packages

Installed automatically by `00_setup.R`:
- **tidyverse** -- data wrangling and visualization
- **lme4 / lmerTest** -- linear mixed-effects models
- **emmeans** -- estimated marginal means and contrasts
- **randomForest** -- variable importance
- **patchwork** -- multi-panel figures
- **plantecophys** -- A/Ci curve fitting (used in 06b_aci_refit.R)

## Contact

Jon Gewirtzman
Yale University
jonathan.gewirtzman@yale.edu
