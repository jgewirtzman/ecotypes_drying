# Ecotypes drying experiment

Analysis code and data for examining Eriophorum vaginatum tussock responses to experimental water treatments across a latitudinal gradient in northern Alaska.

## Project Overview

This repository contains the analysis code for a study examining how Arctic tundra plants respond to different soil moisture conditions. The experiment collected tussocks from three sites along a latitudinal gradient (Sagwon, Toolik Lake, and Coldfoot) and subjected them to three water treatments (Wet, Dry, and Deep) in a common garden experiment.

### Measurements analyzed:
- Soil moisture dynamics
- Plant phenology and growth
- Canopy development (NDVI/LAI)
- Gas exchange (Amax, gs)
- Water potential
- A/Ci curve parameters
- Carbon isotopes (δ13C)
- Productivity relationships

## Repository Structure

```
project/
├── R/                          # Analysis scripts
│   ├── 00_setup.R             # Setup: packages, settings
│   ├── 01_soil_moisture.R     # Soil moisture analysis
│   ├── 02_phenology.R         # Phenology analysis
│   ├── 03_canopy.R            # LAI/NDVI analysis
│   ├── 04_gas_exchange.R      # Photosynthesis analysis
│   ├── 05_water_potential.R   # Water potential
│   ├── 06_aci_curves.R        # A/Ci curves
│   ├── 07_isotopes.R          # δ13C analysis
│   ├── 08_combined_phys.R     # Combined physiology plots
│   ├── 09_productivity.R      # NDVI-GWC relationships
│   ├── 10_model_diag.R        # Model diagnostics
│   └── utils/                 # Utility functions
│       ├── themes.R           # ggplot themes
│       ├── stats.R            # Statistical functions
│       └── plotting.R         # Plotting functions
├── data/
│   ├── raw/                   # Original data files
│   └── processed/             # Cleaned data
└── output/
    ├── figures/               # Generated figures
    └── tables/                # Generated tables
```

## Setup and Usage

1. Clone this repository:
```bash
git clone https://github.com/yourusername/arctic-grass-water.git
cd arctic-grass-water
```

2. Open the R project file (`arctic-grass-water.Rproj`) in RStudio

3. Install required packages:
```r
source("R/00_setup.R")  # This will install and load all necessary packages
```

4. Place raw data files in `data/raw/`:
- `gwc_calibration.csv`
- `Ecotypes2017_Drying_Soil Moisture.csv`
- `Ecotypes2017_Drying_NDVI.csv`
- `Ecotypes2017_Drying_Phenology.csv`
- `Ecotypes2017_Drying_Amax.csv`
- `Ecotypes2017_Drying_Water Potential.csv`
- `A_Ci Outputs_Corrected.csv`
- `Ecotypes2017_Drying_d13C.csv`

5. Run analysis scripts in numerical order (00-10)

## Required R Packages

Core packages (installed automatically by 00_setup.R):
- tidyverse (data manipulation and visualization)
- lme4, lmerTest (mixed models)
- emmeans (estimated marginal means)
- ggridges (ridge plots)
- patchwork (combining plots)
- and others (see 00_setup.R for complete list)

## Output

Analysis outputs are organized in:
- `output/figures/`: PDF figures of all analyses
- `output/tables/`: Statistical summaries and model results

## Contact

Jon Gewirtzman
Yale University  
jonathan.gewirtzman@yale.edu


## Citation

[citation information when available]
