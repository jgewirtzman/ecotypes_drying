# 01_soil_moisture.R
source("R/00_setup.R")

# Read and process calibration data
curve <- read_csv("data/raw/gwc_calibration.csv")

# Extract empty pot weights
empty_pot_data <- curve[grep("Empty Pot", curve$Date), ]
empty_pot_weights <- setNames(empty_pot_data$`Mass (g)`, empty_pot_data$`Pot #`)

# Process calibration data
curve <- curve %>%
  filter(!grepl("Empty Pot", Date)) %>%
  mutate(
    `Wet Soil Mass` = `Mass (g)` - empty_pot_weights[`Pot #`]
  ) %>%
  group_by(`Pot #`) %>%
  mutate(
    `Dry Soil Mass` = min(`Wet Soil Mass`),
    `Water Content` = `Wet Soil Mass` - `Dry Soil Mass`,
    `Gravimetric Water Content` = `Water Content` / `Dry Soil Mass`
  ) %>%
  ungroup() %>%
  na.omit()

# Fit calibration model
calibration_model <- lm(`Gravimetric Water Content` ~ poly(`PER (mean)`, 3), data = curve)

# Define the formula for the polynomial regression
formula <- y ~ poly(x, 3)

# Plot the calibration curve with the equation and R²
calibration_plot <- ggplot(curve, aes(x = `PER (mean)`, y = `Gravimetric Water Content`)) +
  geom_point(aes(color = `Pot #`)) +
  stat_smooth(method = "lm", formula = formula, se = TRUE, color = "black") +
  stat_poly_eq(
    formula = formula,
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")),
    parse = TRUE,
    label.x.npc = "right",
    label.y.npc = "top"
  ) +
  scale_color_viridis_d() +
  labs(
    title = "Soil Moisture Calibration Curve",
    x = "PER (mean)",
    y = "Gravimetric Water Content (GWC)"
  ) +
  theme_minimal()

# Display the plot
calibration_plot


# Save calibration plot
ggsave("output/figures/calibration_curve.pdf", calibration_plot, width = 6, height = 6, units="in", dpi=300)

# Read and process soil moisture data
soil_gwc <- read_csv("data/raw/Ecotypes2017_Drying_Soil Moisture.csv")

# Rename PER column to match calibration data
names(soil_gwc)[names(soil_gwc) == "PER"] <- "PER (mean)"

# Predict GWC using calibration model
soil_gwc$GWC <- predict(calibration_model, newdata = soil_gwc)

# Clean NA values if present
soil_gwc <- soil_gwc %>%
  filter(!is.na(`Treatment Name`))

soil_gwc <- soil_gwc %>%
  filter(!is.na(`Treatment Name`)) %>%  # Remove NA treatments
  mutate(
    Population = factor(Population, 
                        levels = c("SG", "TL", "CF"), 
                        labels = c("Sagwon", "Toolik", "Coldfoot")),
    `Treatment Name` = factor(`Treatment Name`, 
                              levels = c("Wet", "Dry", "Deep"))
  ) %>%
  droplevels()  # Drop any unused factor levels


# Fit mixed model and get summaries
lmm_gwc <- lmer(GWC ~ `Treatment Name` * Population + (1 | Tag), data = soil_gwc)

# Print model summaries
cat("Mixed Model Summary:\n")
print(summary(lmm_gwc))

cat("\nANOVA Table:\n")
anova_lmm <- anova(lmm_gwc)
print(anova_lmm)

# Calculate EMMs
gwc_emm <- emmeans(lmm_gwc, ~ `Treatment Name` | Population)
emm_df <- as.data.frame(gwc_emm)

# Get pairwise contrasts
contrast_output <- contrast(gwc_emm, interaction = "pairwise")

# Save statistical results
sink("output/tables/soil_moisture_stats.txt")
cat("Mixed Model Summary:\n")
print(summary(lmm_gwc))
cat("\nANOVA Table:\n")
print(anova_lmm)
cat("\nEstimated Marginal Means:\n")
print(summary(gwc_emm))
cat("\nPairwise Contrasts:\n")
print(contrast_output)
sink()

# Calculate raw means
raw_means <- soil_gwc %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_gwc = mean(GWC, na.rm = TRUE) * 100,
    se_gwc = sd(GWC, na.rm = TRUE) / sqrt(n()) * 100,
    .groups = "drop"
  )

treatment_means <- soil_gwc %>%
  group_by(`Treatment Name`) %>%
  summarize(
    mean_gwc = mean(GWC, na.rm = TRUE) * 100,
    se_gwc = sd(GWC, na.rm = TRUE) / sqrt(n()) * 100
  )

# Time series plot
ts_plot <- ggplot(soil_gwc, aes(x = DOY, y = GWC * 100, color = `Treatment Name`)) +
  geom_line(aes(group = Tag), alpha = 0.3, size = 0.5) +
  stat_summary(fun = mean, geom = "line", size = 1.5) +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", 
               aes(fill = `Treatment Name`), color = NA, alpha = 0.3) +
  scale_color_manual(values = c("Wet" = "#1f77b4", "Dry" = "#8c564b", "Deep" = "#d62728"),
                     name = "Treatment") +
  scale_fill_manual(values = c("Wet" = "#1f77b4", "Dry" = "#8c564b", "Deep" = "#d62728"),
                    guide = "none") +
  labs(x = "Day of Year", y = "GWC (%)") +
  theme_classic() +
  facet_wrap(~ Population, ncol = 1)

# Bar plot for raw means
bar_plot_raw <- ggplot(raw_means, aes(x = `Treatment Name`, y = mean_gwc, fill = `Treatment Name`)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_gwc - se_gwc, ymax = mean_gwc + se_gwc), width = 0.2) +
  scale_fill_manual(values = c("Wet" = "#1f77b4", "Dry" = "#8c564b", "Deep" = "#d62728"),
                    guide = "none") +
  labs(y = "Mean GWC (%)", x = "Treatment") +
  facet_wrap(~ Population, ncol = 1) +
  theme_classic() +
  theme(legend.position = "none")

# Combine plots with shared legend
combined_plot <- ts_plot + bar_plot_raw + 
  plot_layout(ncol = 2, widths = c(2, 1)) &
  theme(legend.position = "top")

# Save plots
ggsave("output/figures/soil_moisture_combined.pdf", combined_plot, width = 8, height = 6, units="in", dpi=600)
#ggsave("output/figures/soil_moisture_timeseries.pdf", ts_plot, width = 8, height = 8)
#ggsave("output/figures/soil_moisture_bars.pdf", bar_plot_raw, width = 4, height = 8)

# Save summary statistics
write_csv(raw_means, "output/tables/soil_moisture_raw_means.csv")
write_csv(treatment_means, "output/tables/soil_moisture_treatment_means.csv")


# Calculate means per population and treatment (which defines the block)
block_means <- soil_gwc %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_gwc = mean(GWC, na.rm = TRUE) * 100,
    sd_gwc = sd(GWC, na.rm = TRUE) * 100,
    se_gwc = sd(GWC, na.rm = TRUE) / sqrt(n()) * 100,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(Block = paste(Population, `Treatment Name`, sep = "_"))

# Calculate treatment differences within each population
treatment_diffs <- block_means %>%
  group_by(Population) %>%
  arrange(Population, `Treatment Name`) %>%
  summarize(
    Wet_minus_Dry = mean_gwc[`Treatment Name` == "Wet"] - mean_gwc[`Treatment Name` == "Dry"],
    Wet_minus_Deep = mean_gwc[`Treatment Name` == "Wet"] - mean_gwc[`Treatment Name` == "Deep"],
    Dry_minus_Deep = mean_gwc[`Treatment Name` == "Dry"] - mean_gwc[`Treatment Name` == "Deep"],
    .groups = "drop"
  )

# Calculate population differences within each treatment
pop_diffs <- block_means %>%
  group_by(`Treatment Name`) %>%
  arrange(`Treatment Name`, Population) %>%
  summarize(
    Sagwon_minus_Toolik = mean_gwc[Population == "Sagwon"] - mean_gwc[Population == "Toolik"],
    Sagwon_minus_Coldfoot = mean_gwc[Population == "Sagwon"] - mean_gwc[Population == "Coldfoot"],
    Toolik_minus_Coldfoot = mean_gwc[Population == "Toolik"] - mean_gwc[Population == "Coldfoot"],
    .groups = "drop"
  )

# Save results
write_csv(block_means, "output/tables/soil_moisture_block_means.csv")
write_csv(treatment_diffs, "output/tables/treatment_differences.csv")
write_csv(pop_diffs, "output/tables/population_differences.csv")

# Print summary
cat("\nBlock Means (Population × Treatment):\n")
print(block_means)

cat("\nTreatment Differences within Population:\n")
print(treatment_diffs)

cat("\nPopulation Differences within Treatment:\n")
print(pop_diffs)



# Calculate treatment differences within each population (absolute and %)
treatment_diffs <- block_means %>%
  group_by(Population) %>%
  arrange(Population, `Treatment Name`) %>%
  summarize(
    # Absolute differences
    Wet_minus_Dry = mean_gwc[`Treatment Name` == "Wet"] - mean_gwc[`Treatment Name` == "Dry"],
    Wet_minus_Deep = mean_gwc[`Treatment Name` == "Wet"] - mean_gwc[`Treatment Name` == "Deep"],
    Dry_minus_Deep = mean_gwc[`Treatment Name` == "Dry"] - mean_gwc[`Treatment Name` == "Deep"],
    # Percent differences (relative to baseline)
    Wet_minus_Dry_pct = (mean_gwc[`Treatment Name` == "Wet"] - mean_gwc[`Treatment Name` == "Dry"]) / mean_gwc[`Treatment Name` == "Dry"] * 100,
    Wet_minus_Deep_pct = (mean_gwc[`Treatment Name` == "Wet"] - mean_gwc[`Treatment Name` == "Deep"]) / mean_gwc[`Treatment Name` == "Deep"] * 100,
    Dry_minus_Deep_pct = (mean_gwc[`Treatment Name` == "Dry"] - mean_gwc[`Treatment Name` == "Deep"]) / mean_gwc[`Treatment Name` == "Deep"] * 100,
    .groups = "drop"
  )

# Calculate population differences within each treatment (absolute and %)
pop_diffs <- block_means %>%
  group_by(`Treatment Name`) %>%
  arrange(`Treatment Name`, Population) %>%
  summarize(
    # Absolute differences
    Sagwon_minus_Toolik = mean_gwc[Population == "Sagwon"] - mean_gwc[Population == "Toolik"],
    Sagwon_minus_Coldfoot = mean_gwc[Population == "Sagwon"] - mean_gwc[Population == "Coldfoot"],
    Toolik_minus_Coldfoot = mean_gwc[Population == "Toolik"] - mean_gwc[Population == "Coldfoot"],
    # Percent differences (relative to baseline)
    Sagwon_minus_Toolik_pct = (mean_gwc[Population == "Sagwon"] - mean_gwc[Population == "Toolik"]) / mean_gwc[Population == "Toolik"] * 100,
    Sagwon_minus_Coldfoot_pct = (mean_gwc[Population == "Sagwon"] - mean_gwc[Population == "Coldfoot"]) / mean_gwc[Population == "Coldfoot"] * 100,
    Toolik_minus_Coldfoot_pct = (mean_gwc[Population == "Toolik"] - mean_gwc[Population == "Coldfoot"]) / mean_gwc[Population == "Coldfoot"] * 100,
    .groups = "drop"
  )

# Save results
write_csv(block_means, "output/tables/soil_moisture_block_means.csv")
write_csv(treatment_diffs, "output/tables/treatment_differences.csv")
write_csv(pop_diffs, "output/tables/population_differences.csv")

# Print summary
cat("\nBlock Means (Population × Treatment):\n")
print(block_means)

cat("\nTreatment Differences within Population:\n")
print(treatment_diffs)

cat("\nPopulation Differences within Treatment:\n")
print(pop_diffs)


# Calculate CV for treatments (across all populations)
treatment_cv <- soil_gwc %>%
  group_by(`Treatment Name`) %>%
  summarize(mean_gwc = mean(GWC, na.rm = TRUE) * 100) %>%
  summarize(
    grand_mean = mean(mean_gwc),
    sd = sd(mean_gwc),
    cv = (sd / grand_mean) * 100
  )

cat("\nCoefficient of Variation for Treatments:\n")
print(treatment_cv)

# Calculate CV for populations (across all treatments)
population_cv <- soil_gwc %>%
  group_by(Population) %>%
  summarize(mean_gwc = mean(GWC, na.rm = TRUE) * 100) %>%
  summarize(
    grand_mean = mean(mean_gwc),
    sd = sd(mean_gwc),
    cv = (sd / grand_mean) * 100
  )

cat("\nCoefficient of Variation for Populations:\n")
print(population_cv)

# Alternative: Calculate CV from the block_means you already have
treatment_cv_from_blocks <- block_means %>%
  group_by(`Treatment Name`) %>%
  summarize(mean_gwc = mean(mean_gwc)) %>%
  summarize(
    grand_mean = mean(mean_gwc),
    sd = sd(mean_gwc),
    cv = (sd / grand_mean) * 100
  )

population_cv_from_blocks <- block_means %>%
  group_by(Population) %>%
  summarize(mean_gwc = mean(mean_gwc)) %>%
  summarize(
    grand_mean = mean(mean_gwc),
    sd = sd(mean_gwc),
    cv = (sd / grand_mean) * 100
  )

cat("\nFrom block means:\n")
cat("Treatment CV:", round(treatment_cv_from_blocks$cv, 1), "%\n")
cat("Population CV:", round(population_cv_from_blocks$cv, 1), "%\n")