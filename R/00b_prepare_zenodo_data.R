# 00b_prepare_zenodo_data.R
# Prepares clean processed data files for Zenodo upload.
# Reads from data/raw/, writes to data/processed/.
# Run standalone (does NOT require 00_setup.R).

library(readr)
library(dplyr)

outdir <- "data/processed"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Helper: read CSV with BOM handling, drop all-NA columns and all-NA rows
read_clean <- function(path) {
  df <- read_csv(path, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
  # Drop columns that are entirely NA (trailing empties from Excel)
  df <- df %>% select(where(~ !all(is.na(.))))
  # Drop unnamed index columns (often read as "...1" or empty-name)
  bad_cols <- grepl("^\\.\\.\\.\\d+$", names(df)) | names(df) == ""
  df <- df[, !bad_cols, drop = FALSE]
  # Drop rows where ALL non-identifier columns are NA
  df <- df %>% filter(rowSums(!is.na(pick(everything()))) > 0)
  as.data.frame(df)
}

# Helper: standardize the common identifier columns
clean_ids <- function(df) {
  # Rename standard columns to snake_case
  renames <- c(
    "Population"     = "population",
    "Treatment"      = "treatment_code",
    "Treatment Name" = "treatment",
    "TreatmentName"  = "treatment",
    "Replicate"      = "replicate",
    "Tag"            = "tag",
    "Date"           = "date",
    "DOY"            = "doy"
  )
  for (old in names(renames)) {
    if (old %in% names(df)) {
      names(df)[names(df) == old] <- renames[old]
    }
  }
  df
}

cat("Preparing clean data files for Zenodo...\n\n")

# ============================================================================
# 1. SOIL MOISTURE (+ GWC calibration)
# ============================================================================
cat("1. Soil moisture...\n")

# Build calibration model from calibration data
cal <- read_clean("data/raw/gwc_calibration.csv")
names(cal) <- c("date", "pot", "mass_g", "weight_minus_pot", "per_1", "per_2", "per_3", "per_mean", "gwc")
# Extract empty pot weights
empty_pots <- cal %>% filter(grepl("Empty Pot", date))
empty_weights <- setNames(empty_pots$mass_g, empty_pots$pot)
# Process calibration
cal <- cal %>%
  filter(!grepl("Empty Pot", date)) %>%
  mutate(
    wet_soil_mass = mass_g - empty_weights[pot],
  ) %>%
  group_by(pot) %>%
  mutate(
    dry_soil_mass = min(wet_soil_mass),
    water_content = wet_soil_mass - dry_soil_mass,
    gwc_calc = water_content / dry_soil_mass
  ) %>%
  ungroup() %>%
  na.omit()
cal_model <- lm(gwc_calc ~ poly(per_mean, 3), data = cal)

# Read soil moisture data
sm <- read_clean("data/raw/Ecotypes2017_Drying_Soil Moisture.csv")
sm <- clean_ids(sm)
# The raw file has a PER column (and maybe GWC already); rename PER
if ("PER" %in% names(sm)) names(sm)[names(sm) == "PER"] <- "per"
if ("GWC" %in% names(sm)) names(sm)[names(sm) == "GWC"] <- "gwc_raw"
# Compute GWC from calibration
sm$per_mean <- sm$per  # match calibration model variable name
sm$gwc <- predict(cal_model, newdata = sm)
sm$per_mean <- NULL
# Select final columns
sm <- sm %>%
  filter(!is.na(treatment)) %>%
  select(population, treatment, replicate, tag, date, doy, per, gwc)
write_csv(sm, file.path(outdir, "soil_moisture.csv"))
cat("   ->", nrow(sm), "rows\n")

# ============================================================================
# 2. PHENOLOGY
# ============================================================================
cat("2. Phenology...\n")

phen <- read_clean("data/raw/Ecotypes2017_Drying_Phenology.csv")
phen <- clean_ids(phen)

# Fix duplicate column name: readr renames duplicates to "L1 Total Length...8" and
# "L1 Total Length...10". The second one should be "L2 Total Length".
raw_names <- names(phen)
if (any(grepl("L1 Total Length\\.\\.\\.10", raw_names))) {
  names(phen)[grepl("L1 Total Length\\.\\.\\.10", names(phen))] <- "L2 Total Length"
  names(phen)[grepl("L1 Total Length\\.\\.\\.8", names(phen))]  <- "L1 Total Length"
} else {
  # Fallback for base R read
  l1_idx <- which(raw_names == "L1 Total Length")
  if (length(l1_idx) == 2) names(phen)[l1_idx[2]] <- "L2 Total Length"
}

# Rename leaf columns to snake_case
leaf_renames <- c(
  "L1 Total Length" = "l1_total_length", "L1 Green Length" = "l1_green_length",
  "L2 Total Length" = "l2_total_length", "L2 Green Length" = "l2_green_length",
  "L3 Total Length" = "l3_total_length", "L3 Green Length" = "l3_green_length",
  "L4 Total Length" = "l4_total_length", "L4 Green Length" = "l4_green_length",
  "L5 Total Length" = "l5_total_length", "L5 Green Length" = "l5_green_length",
  "L6 Total Length" = "l6_total_length", "L6 Green Length" = "l6_green_length",
  "Tiller Total Green Length" = "tiller_total_green_length",
  "Average Green Leaf Length"  = "avg_green_leaf_length",
  "Tiller % Green"            = "tiller_pct_green"
)
for (old in names(leaf_renames)) {
  if (old %in% names(phen)) {
    names(phen)[names(phen) == old] <- leaf_renames[old]
  }
}

phen <- phen %>%
  filter(!is.na(treatment)) %>%
  select(population, treatment, replicate, tag, date, doy,
         l1_total_length, l1_green_length, l2_total_length, l2_green_length,
         l3_total_length, l3_green_length, l4_total_length, l4_green_length,
         l5_total_length, l5_green_length, l6_total_length, l6_green_length,
         tiller_total_green_length, avg_green_leaf_length, tiller_pct_green)
write_csv(phen, file.path(outdir, "phenology.csv"))
cat("   ->", nrow(phen), "rows\n")

# ============================================================================
# 3. NDVI / CANOPY
# ============================================================================
cat("3. NDVI / canopy...\n")

ndvi <- read_clean("data/raw/Ecotypes2017_Drying_NDVI.csv")
ndvi <- clean_ids(ndvi)
names(ndvi)[names(ndvi) == "NDVI"] <- "ndvi"

# Add derived canopy metrics (same formulas as 03_canopy.R)
ndvi <- ndvi %>%
  filter(!is.na(treatment)) %>%
  mutate(
    lai = 0.03 * exp(7.65 * ndvi),
    biomass_g = 0.0256 * exp(5.32 * ndvi)
  ) %>%
  select(population, treatment, replicate, tag, date, doy, ndvi, lai, biomass_g)
write_csv(ndvi, file.path(outdir, "ndvi.csv"))
cat("   ->", nrow(ndvi), "rows\n")

# ============================================================================
# 4. GAS EXCHANGE (Amax)
# ============================================================================
cat("4. Gas exchange...\n")

gas <- read_clean("data/raw/Ecotypes2017_Drying_Amax.csv")
gas <- clean_ids(gas)

# Rename measurement columns to snake_case with units
gas_renames <- c(
  "Obs"    = "obs",
  "Photo"  = "photo_umol_m2_s",
  "Cond"   = "cond_mol_m2_s",
  "Ci"     = "ci_umol_mol",
  "Trmmol" = "trans_mmol_m2_s",
  "VpdL"   = "vpd_kpa",
  "Tleaf"  = "tleaf_c",
  "Tair"   = "tair_c",
  "CO2R"   = "co2_ref_umol_mol",
  "CO2S"   = "co2_sample_umol_mol",
  "PARi"   = "par_umol_m2_s",
  "Press"  = "press_kpa"
)
for (old in names(gas_renames)) {
  if (old %in% names(gas)) {
    names(gas)[names(gas) == old] <- gas_renames[old]
  }
}

gas <- gas %>%
  filter(!is.na(treatment)) %>%
  select(population, treatment, replicate, tag, date, doy, obs,
         photo_umol_m2_s, cond_mol_m2_s, ci_umol_mol, trans_mmol_m2_s,
         vpd_kpa, tleaf_c, tair_c, co2_ref_umol_mol, co2_sample_umol_mol,
         par_umol_m2_s, press_kpa)
write_csv(gas, file.path(outdir, "gas_exchange.csv"))
cat("   ->", nrow(gas), "rows\n")

# ============================================================================
# 5. WATER POTENTIAL
# ============================================================================
cat("5. Water potential...\n")

wpot <- read_clean("data/raw/Ecotypes2017_Drying_Water Potential.csv")
wpot <- clean_ids(wpot)
names(wpot)[names(wpot) == "Pressure (Bar)"] <- "pressure_bar"

wpot <- wpot %>%
  filter(!is.na(pressure_bar)) %>%
  mutate(
    pressure_mpa = -pressure_bar * 0.1  # Convert bar to MPa, flip sign (predawn)
  ) %>%
  select(population, treatment, replicate, tag, date, doy, pressure_bar, pressure_mpa)
write_csv(wpot, file.path(outdir, "water_potential.csv"))
cat("   ->", nrow(wpot), "rows\n")

# ============================================================================
# 6. A/Ci CURVES (raw LI-6400 data)
# ============================================================================
cat("6. A/Ci curves (raw)...\n")

aci <- read_clean("data/raw/Ecotypes2017_Drying_A-Ci.csv")
aci <- clean_ids(aci)

# Rename columns
aci_renames <- c(
  "Unique ID" = "curve_id",
  "Obs"       = "obs",
  "Photo"     = "photo_umol_m2_s",
  "Cond"      = "cond_mol_m2_s",
  "Ci"        = "ci_umol_mol",
  "Trmmol"    = "trans_mmol_m2_s",
  "VpdL"      = "vpd_kpa",
  "Tleaf"     = "tleaf_c",
  "TBlk"      = "tblock_c",
  "Tair"      = "tair_c",
  "CO2R"      = "co2_ref_umol_mol",
  "CO2S"      = "co2_sample_umol_mol",
  "PARi"      = "par_umol_m2_s",
  "Press"     = "press_kpa",
  "Area"      = "leaf_area_cm2"
)
for (old in names(aci_renames)) {
  if (old %in% names(aci)) {
    names(aci)[names(aci) == old] <- aci_renames[old]
  }
}

aci <- aci %>%
  filter(!is.na(treatment)) %>%
  select(population, treatment, replicate, tag, date, curve_id, obs,
         photo_umol_m2_s, cond_mol_m2_s, ci_umol_mol, trans_mmol_m2_s,
         vpd_kpa, tleaf_c, tblock_c, tair_c,
         co2_ref_umol_mol, co2_sample_umol_mol, par_umol_m2_s,
         press_kpa, leaf_area_cm2)
write_csv(aci, file.path(outdir, "aci_curves.csv"))
cat("   ->", nrow(aci), "rows\n")

# ============================================================================
# 7. A/Ci FITTED PARAMETERS (both Tcorrect methods)
# ============================================================================
cat("7. A/Ci fitted parameters...\n")

# Read the Tcorrect=TRUE values from original outputs
orig <- read_clean("data/raw/A_Ci Ouputs_Corrected.csv")
orig <- clean_ids(orig)
names(orig)[names(orig) == "Unique.ID"] <- "curve_id"
names(orig)[names(orig) == "Vcmax"]    <- "vcmax_25c"
names(orig)[names(orig) == "Jmax"]     <- "jmax_25c"
names(orig)[names(orig) == "Rd"]       <- "rd_25c"
names(orig)[names(orig) == "Vcmax_SE"] <- "vcmax_25c_se"
names(orig)[names(orig) == "Jmax_SE"]  <- "jmax_25c_se"
names(orig)[names(orig) == "Rd_SE"]    <- "rd_25c_se"

# Read the Tcorrect=FALSE values from the refit output
refit_file <- "output/tables/aci_refit_validation.csv"
if (file.exists(refit_file)) {
  refit <- read_csv(refit_file, show_col_types = FALSE)
  # This file has Vcmax_refit (Tcorrect=TRUE) - we need Tcorrect=FALSE from
  # the summary tables instead. Use the per-treatment summary to get curve-level.
}

# Better approach: read the refit script's per-curve output if available,
# or re-derive from the existing tables. The refit_validation.csv has the
# Tcorrect=TRUE values. For Tcorrect=FALSE, read from summary_by_pop_treatment.
# Actually, the cleanest approach: run fitacis here if plantecophys is available.
has_plantecophys <- requireNamespace("plantecophys", quietly = TRUE)

if (has_plantecophys) {
  cat("   Refitting with Tcorrect=FALSE for per-curve values...\n")
  library(plantecophys)
  aci_fit <- read_csv("data/raw/Ecotypes2017_Drying_A-Ci.csv", show_col_types = FALSE) %>%
    rename(Curve = `Unique ID`) %>%
    dplyr::select(Curve, Photo, Ci, Tleaf, PARi, Press) %>%
    as.data.frame()
  fits_uncorr <- tryCatch(
    fitacis(aci_fit, "Curve", Tcorrect = FALSE, quiet = TRUE),
    error = function(e) {
      cat("   fitacis error:", e$message, "\n")
      NULL
    }
  )
  if (!is.null(fits_uncorr)) {
    coefs_uncorr <- as.data.frame(coef(fits_uncorr))
    # Curve ID is in the "Curve" column (not rownames)
    names(coefs_uncorr)[names(coefs_uncorr) == "Curve"]    <- "curve_id"
    names(coefs_uncorr)[names(coefs_uncorr) == "Vcmax"]    <- "vcmax_meas_temp"
    names(coefs_uncorr)[names(coefs_uncorr) == "Jmax"]     <- "jmax_meas_temp"
    names(coefs_uncorr)[names(coefs_uncorr) == "Rd"]       <- "rd_meas_temp"
    names(coefs_uncorr)[names(coefs_uncorr) == "Vcmax_SE"] <- "vcmax_meas_temp_se"
    names(coefs_uncorr)[names(coefs_uncorr) == "Jmax_SE"]  <- "jmax_meas_temp_se"
    names(coefs_uncorr)[names(coefs_uncorr) == "Rd_SE"]    <- "rd_meas_temp_se"
  }
} else {
  cat("   plantecophys not available; Tcorrect=FALSE columns will be NA\n")
  fits_uncorr <- NULL
}

# Build the combined parameters table
params <- orig %>%
  select(curve_id, population, treatment, replicate, tag, date,
         vcmax_25c, jmax_25c, rd_25c, vcmax_25c_se, jmax_25c_se, rd_25c_se)

if (!is.null(fits_uncorr)) {
  uncorr_cols <- coefs_uncorr %>%
    select(curve_id, vcmax_meas_temp, jmax_meas_temp, rd_meas_temp,
           vcmax_meas_temp_se, jmax_meas_temp_se, rd_meas_temp_se)
  params <- params %>% left_join(uncorr_cols, by = "curve_id")
} else {
  params$vcmax_meas_temp <- NA
  params$jmax_meas_temp  <- NA
  params$rd_meas_temp    <- NA
  params$vcmax_meas_temp_se <- NA
  params$jmax_meas_temp_se  <- NA
  params$rd_meas_temp_se    <- NA
}

write_csv(params, file.path(outdir, "aci_parameters.csv"))
cat("   ->", nrow(params), "rows\n")

# ============================================================================
# 8. CARBON ISOTOPES (d13C)
# ============================================================================
cat("8. Carbon isotopes...\n")

iso <- read_clean("data/raw/Ecotypes2017_Drying_d13C.csv")
iso <- clean_ids(iso)
names(iso)[names(iso) == "Sample ID"]                 <- "sample_id"
names(iso)[names(iso) == "d13C"]                      <- "d13c_permil"
names(iso)[names(iso) == "Actual (mg) Mass for d13C"] <- "sample_mass_mg"

iso <- iso %>%
  filter(!is.na(d13c_permil)) %>%
  select(population, treatment, replicate, sample_id, d13c_permil, sample_mass_mg)
write_csv(iso, file.path(outdir, "carbon_isotopes.csv"))
cat("   ->", nrow(iso), "rows\n")

# ============================================================================
# Summary
# ============================================================================
cat("\n=== Done ===\n")
cat("Files written to", outdir, ":\n")
for (f in list.files(outdir, pattern = "\\.csv$")) {
  sz <- file.size(file.path(outdir, f))
  cat(sprintf("  %-25s %6.1f KB\n", f, sz / 1024))
}
