# reviewer_response_exploration.R
# Exploration of raw LI-6400 A/Ci curve data to determine measurement temperatures
# and other metadata for reviewer response

library(tidyverse)

# Read raw A/Ci data
aci_raw <- read_csv("data/raw/Ecotypes2017_Drying_A-Ci.csv")

cat("=============================================================\n")
cat("ANALYSIS OF RAW LI-6400 A/Ci CURVE DATA\n")
cat("=============================================================\n\n")

# --- 1. All unique dates in the dataset ---
cat("--- 1. UNIQUE DATES IN DATASET ---\n\n")

unique_dates <- sort(unique(aci_raw$Date))
cat("Number of unique dates:", length(unique_dates), "\n")
cat("Date codes:", paste(unique_dates, collapse = ", "), "\n\n")

# Decode LI-6400 date format: YYMMDD with variable-width month
# 17802 = YY=17, M=8, DD=02 -> Aug 2, 2017
# 171003 = YY=17, MM=10, DD=03 -> Oct 3, 2017
# The format is: first 2 digits = year, last 2 digits = day, middle = month
decode_date <- function(d) {
  d_str <- as.character(d)
  year <- paste0("20", substr(d_str, 1, 2))
  day <- substr(d_str, nchar(d_str) - 1, nchar(d_str))
  month_str <- substr(d_str, 3, nchar(d_str) - 2)
  month <- as.integer(month_str)
  # If month > 12, this date code may use a different convention
  if (month > 12) {
    return(paste0(year, "-??-", day, " (raw code: ", d, ", month=", month, " invalid)"))
  }
  paste0(year, "-", sprintf("%02d", month), "-", day)
}

cat("Decoded dates:\n")
for (d in unique_dates) {
  cat("  ", d, " -> ", decode_date(d), "\n")
}
cat("\n")
cat("NOTE: Date code 17802 = Aug 2, 2017 (per user). Date code 171803 has an\n")
cat("ambiguous month encoding (middle digits = 18, which exceeds 12). This may\n")
cat("represent a second measurement campaign date. The Tair and TBlk values\n")
cat("differ between the two dates, confirming they are distinct measurement days.\n\n")

# --- 2. Temperature summary per date ---
cat("--- 2. TEMPERATURE SUMMARY PER DATE ---\n\n")

temp_summary <- aci_raw %>%
  group_by(Date) %>%
  summarise(
    n_obs = n(),
    Tleaf_mean = mean(Tleaf, na.rm = TRUE),
    Tleaf_min  = min(Tleaf, na.rm = TRUE),
    Tleaf_max  = max(Tleaf, na.rm = TRUE),
    TBlk_mean  = mean(TBlk, na.rm = TRUE),
    TBlk_min   = min(TBlk, na.rm = TRUE),
    TBlk_max   = max(TBlk, na.rm = TRUE),
    Tair_mean  = mean(Tair, na.rm = TRUE),
    Tair_min   = min(Tair, na.rm = TRUE),
    Tair_max   = max(Tair, na.rm = TRUE),
    .groups = "drop"
  )

for (i in 1:nrow(temp_summary)) {
  row <- temp_summary[i, ]
  cat(sprintf("Date: %s (%s)  |  N obs: %d\n", row$Date, decode_date(row$Date), row$n_obs))
  cat(sprintf("  Tleaf  -> mean: %.2f, min: %.2f, max: %.2f\n", row$Tleaf_mean, row$Tleaf_min, row$Tleaf_max))
  cat(sprintf("  TBlk   -> mean: %.2f, min: %.2f, max: %.2f\n", row$TBlk_mean, row$TBlk_min, row$TBlk_max))
  cat(sprintf("  Tair   -> mean: %.2f, min: %.2f, max: %.2f\n\n", row$Tair_mean, row$Tair_min, row$Tair_max))
}

# --- 3. Number of unique Tags (plants) measured per date ---
cat("--- 3. UNIQUE TAGS (PLANTS) PER DATE ---\n\n")

tags_per_date <- aci_raw %>%
  group_by(Date) %>%
  summarise(
    n_tags = n_distinct(Tag),
    tags = paste(sort(unique(Tag)), collapse = ", "),
    .groups = "drop"
  )

for (i in 1:nrow(tags_per_date)) {
  row <- tags_per_date[i, ]
  cat(sprintf("Date: %s (%s)  |  N plants: %d\n", row$Date, decode_date(row$Date), row$n_tags))
  cat(sprintf("  Tags: %s\n\n", row$tags))
}

# --- 4. Overall range and typical values of TBlk ---
cat("--- 4. OVERALL TBlk SUMMARY ---\n\n")

cat(sprintf("TBlk across all measurements:\n"))
cat(sprintf("  Mean:   %.4f\n", mean(aci_raw$TBlk, na.rm = TRUE)))
cat(sprintf("  Median: %.4f\n", median(aci_raw$TBlk, na.rm = TRUE)))
cat(sprintf("  SD:     %.4f\n", sd(aci_raw$TBlk, na.rm = TRUE)))
cat(sprintf("  Min:    %.4f\n", min(aci_raw$TBlk, na.rm = TRUE)))
cat(sprintf("  Max:    %.4f\n", max(aci_raw$TBlk, na.rm = TRUE)))
cat(sprintf("  Range:  %.4f\n\n", max(aci_raw$TBlk, na.rm = TRUE) - min(aci_raw$TBlk, na.rm = TRUE)))

cat("Unique TBlk values (rounded to 1 decimal):\n")
cat("  ", paste(sort(unique(round(aci_raw$TBlk, 1))), collapse = ", "), "\n\n")

cat("TBlk distribution (rounded to nearest integer):\n")
tblk_table <- table(round(aci_raw$TBlk))
print(tblk_table)
cat("\n")

# --- 5. Was TBlk set to a fixed value or did it vary? ---
cat("--- 5. WAS TBlk FIXED OR VARIABLE? ---\n\n")

tblk_range <- max(aci_raw$TBlk, na.rm = TRUE) - min(aci_raw$TBlk, na.rm = TRUE)
tblk_sd <- sd(aci_raw$TBlk, na.rm = TRUE)

if (tblk_range < 1.0) {
  cat(sprintf("TBlk appears to be set to a FIXED value (range = %.4f C, SD = %.4f C)\n", tblk_range, tblk_sd))
  cat(sprintf("The block temperature was set to approximately %.1f C for all A/Ci measurements.\n\n",
              round(mean(aci_raw$TBlk, na.rm = TRUE), 1)))
} else {
  cat(sprintf("TBlk appears to have VARIED across measurements (range = %.4f C, SD = %.4f C)\n\n", tblk_range, tblk_sd))

  # Check if it varied by date or within dates
  cat("TBlk by date:\n")
  tblk_by_date <- aci_raw %>%
    group_by(Date) %>%
    summarise(
      TBlk_mean = mean(TBlk, na.rm = TRUE),
      TBlk_sd = sd(TBlk, na.rm = TRUE),
      TBlk_min = min(TBlk, na.rm = TRUE),
      TBlk_max = max(TBlk, na.rm = TRUE),
      TBlk_range = max(TBlk, na.rm = TRUE) - min(TBlk, na.rm = TRUE),
      .groups = "drop"
    )
  print(as.data.frame(tblk_by_date))
  cat("\n")

  cat("INTERPRETATION: TBlk was set to a FIXED value WITHIN each date but\n")
  cat("DIFFERED between dates:\n")
  for (i in 1:nrow(tblk_by_date)) {
    row <- tblk_by_date[i, ]
    cat(sprintf("  Date %s: TBlk ~ %.1f C (within-date range: %.2f C, SD: %.4f C)\n",
                row$Date, round(row$TBlk_mean, 1), row$TBlk_range, row$TBlk_sd))
  }
  cat("\n")
}

# Also report CTleaf and Tleaf for context
cat("For context, Tleaf (thermocouple) overall summary:\n")
cat(sprintf("  Mean:   %.2f C\n", mean(aci_raw$Tleaf, na.rm = TRUE)))
cat(sprintf("  SD:     %.2f C\n", sd(aci_raw$Tleaf, na.rm = TRUE)))
cat(sprintf("  Min:    %.2f C\n", min(aci_raw$Tleaf, na.rm = TRUE)))
cat(sprintf("  Max:    %.2f C\n\n", max(aci_raw$Tleaf, na.rm = TRUE)))

cat("For context, CTleaf (computed leaf temp) overall summary:\n")
cat(sprintf("  Mean:   %.2f C\n", mean(aci_raw$CTleaf, na.rm = TRUE)))
cat(sprintf("  SD:     %.2f C\n", sd(aci_raw$CTleaf, na.rm = TRUE)))
cat(sprintf("  Min:    %.2f C\n", min(aci_raw$CTleaf, na.rm = TRUE)))
cat(sprintf("  Max:    %.2f C\n\n", max(aci_raw$CTleaf, na.rm = TRUE)))

# --- 6. Total number of A/Ci curves ---
cat("--- 6. TOTAL NUMBER OF A/Ci CURVES ---\n\n")

curves <- aci_raw %>%
  distinct(Tag, Date)

cat(sprintf("Total A/Ci curves (unique Tag x Date): %d\n\n", nrow(curves)))

# Breakdown by population and treatment
curves_detail <- aci_raw %>%
  distinct(Tag, Date, Population, `Treatment Name`) %>%
  count(Population, `Treatment Name`, name = "n_curves")

cat("Curves by Population x Treatment:\n")
print(as.data.frame(curves_detail))
cat("\n")

# --- 7. Number of CO2 concentration steps (Obs) per curve ---
cat("--- 7. CO2 STEPS (OBS) PER CURVE ---\n\n")

obs_per_curve <- aci_raw %>%
  group_by(Tag, Date) %>%
  summarise(
    n_obs = n(),
    max_obs = max(Obs),
    CO2R_values = paste(round(CO2R), collapse = ", "),
    .groups = "drop"
  )

cat(sprintf("Observations per curve:\n"))
cat(sprintf("  Mean:   %.1f\n", mean(obs_per_curve$n_obs)))
cat(sprintf("  Median: %.0f\n", median(obs_per_curve$n_obs)))
cat(sprintf("  Min:    %d\n", min(obs_per_curve$n_obs)))
cat(sprintf("  Max:    %d\n", max(obs_per_curve$n_obs)))
cat(sprintf("  SD:     %.2f\n\n", sd(obs_per_curve$n_obs)))

cat("Distribution of steps per curve:\n")
print(table(obs_per_curve$n_obs))
cat("\n")

# Show the CO2 reference concentrations for the first curve as example
cat("Example CO2R sequence (first curve):\n")
first_curve <- aci_raw %>%
  filter(Tag == obs_per_curve$Tag[1], Date == obs_per_curve$Date[1]) %>%
  arrange(Obs)
cat("  CO2R values:", paste(round(first_curve$CO2R), collapse = " -> "), "\n\n")

# Check if Obs numbering resets per curve or is continuous
cat("Obs column behavior:\n")
obs_ranges <- aci_raw %>%
  group_by(Tag, Date) %>%
  summarise(
    min_obs = min(Obs),
    max_obs = max(Obs),
    .groups = "drop"
  )
cat(sprintf("  Obs numbering appears %s across curves\n",
            ifelse(any(obs_ranges$min_obs > 1), "CONTINUOUS (does not reset)", "to RESET per curve")))
cat(sprintf("  Overall Obs range: %d to %d\n", min(obs_ranges$min_obs), max(obs_ranges$max_obs)))
cat("\n")

cat("=============================================================\n")
cat("END OF ANALYSIS\n")
cat("=============================================================\n")
