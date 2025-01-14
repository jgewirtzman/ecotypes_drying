# 00_setup.R
# Get all loaded packages
loaded_packages <- names(sessionInfo()$otherPkgs)

# Detach each one
if(!is.null(loaded_packages)) {
  sapply(paste0("package:", loaded_packages), detach, character.only = TRUE, unload = TRUE)
}

# Check and install pacman if not already available
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

# Load required packages
pacman::p_load(
  tidyverse,      # Data manipulation and visualization
  nlme,           # Mixed models
  lme4,           # Mixed models
  lmerTest,       # P-values for mixed models
  emmeans,        # Estimated marginal means
  multcomp,       # Multiple comparisons
  ggridges,       # Ridge plots
  viridis,        # Color palettes
  patchwork,      # Combining plots
  boot,           # Bootstrap analyses
  car,           # Statistical functions
  pracma,        # Numerical analysis
  randomForest,  # Random forest models
  EnvStats,      # Environmental statistics
  ggpmisc,
  broom,
  gt,
  broom.mixed,
  ggrepel
)

# Global constants
POPULATION_LEVELS <- c("SG", "TL", "CF")
POPULATION_LABELS <- c("Sagwon", "Toolik", "Coldfoot")
TREATMENT_LEVELS <- c("Wet", "Dry", "Deep")
TREATMENT_COLORS <- c(
  "Wet" = "#1f77b4", 
  "Dry" = "#8c564b", 
  "Deep" = "#d62728"
)
custom_colors_pop <- c("Coldfoot" = "#3E6347",  
                       "Toolik" = "#90C49B", 
                       "Sagwon" = "#79A9C8")  

# Source utility functions
source("R/utils/themes.R")
source("R/utils/stats.R")
source("R/utils/plotting.R")

# Create directories if they don't exist
dirs <- c(
  "data/raw",
  "data/processed",
  "output/figures",
  "output/tables"
)

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}