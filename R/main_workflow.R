# main.R
# Master script to run all analyses in the correct order

# Create required directories
dirs <- c(
  "data/raw",
  "data/processed", 
  "output/figures",
  "output/tables"
)

for (dir in dirs) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Analysis Pipeline
scripts <- c(
  "00_setup.R",
  "01_soil_moisture.R", 
  "02_phenology.R",
  "03_canopy.R",
  "04_gas_exchange.R",
  "05_water_potential.R", 
  "06_aci_curves.R",
  "07_isotopes.R",
  "08_combined_physiology.R",
  "09_productivity.R",
  "10_model_diagnostics.R"
)

# Run each script
for (script in scripts) {
  message("\nRunning ", script)
  source(file.path("R", script))
}

message("\nAnalysis complete")