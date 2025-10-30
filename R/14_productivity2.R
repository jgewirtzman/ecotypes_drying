source("R/00_setup.R")

# Read and process NDVI data, then calculate LAI
ndvi <- read_csv("data/raw/Ecotypes2017_Drying_NDVI.csv") %>%
  # Remove columns with all NA values
  select_if(~ !all(is.na(.))) %>%
  # Remove rows with all NA values
  filter(rowSums(is.na(.)) < ncol(.)) %>%
  # Calculate LAI from NDVI
  mutate(LAI = 0.03 * exp(7.65 * NDVI))

# Calculate LAI AUC
lai_dat <- ndvi %>%
  group_by(Population, `Treatment Name`, Tag) %>%
  summarize(
    lai_auc = pracma::trapz(DOY[!is.na(LAI)], LAI[!is.na(LAI)]),
    .groups = "drop"
  )

# Read and process soil moisture data for AUC
soil_gwc <- read_csv("data/raw/Ecotypes2017_Drying_Soil Moisture.csv") %>%
  select_if(~!all(is.na(.))) %>%
  filter(rowSums(is.na(.)) < ncol(.))

# Calculate GWC AUC
gwc_dat <- soil_gwc %>%
  group_by(Population, `Treatment Name`, Tag) %>%
  summarize(
    gwc_auc = pracma::trapz(DOY[!is.na(GWC)], GWC[!is.na(GWC)]),
    .groups = "drop"
  )

# Merge LAI and GWC AUC datasets
lai_gwc_merged <- left_join(lai_dat, gwc_dat, 
                            by = c("Population", "Treatment Name", "Tag"))

# Update population names and order
lai_gwc_merged <- lai_gwc_merged %>%
  mutate(
    Population = factor(Population,
                        levels = POPULATION_LEVELS,
                        labels = POPULATION_LABELS),
    `Treatment Name` = factor(`Treatment Name`, levels = TREATMENT_LEVELS)
  )

# Calculate z-scores
lai_gwc_merged <- lai_gwc_merged %>%
  mutate(
    gwc_auc_z = scale(gwc_auc),
    lai_auc_z = scale(lai_auc)
  )

# Custom colors
custom_colors <- c(
  "Coldfoot" = "#3E6347",  
  "Toolik" = "#90C49B", 
  "Sagwon" = "#79A9C8"
)

# Create productivity plot with LAI
productivity_plot <- ggplot(lai_gwc_merged, 
                            aes(y = lai_auc_z, x = gwc_auc_z, 
                                color = Population, fill = Population)) +
  geom_point(size = 3.5, shape = 21, stroke = 1.2, alpha = 0.7) +
  geom_smooth(method = "lm", size = 1.5, 
              se = FALSE, fullrange = TRUE, alpha = 0.2) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  labs(
    x = "Growing Season Soil Moisture\n(GWC AUC Z-Score)",
    y = "Cumulative Leaf Area Index\n(LAI AUC Z-Score)",
    color = "Population",
    fill = "Population"
  ) +
  theme(
    legend.position = "top",
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  guides(color = guide_legend(title = "Population"))

# Save plot
ggsave("output/figures/productivity_relationship_lai.pdf", 
       productivity_plot, width = 6, height = 6, dpi=600)

# Fit models and save statistics
full_model <- lm(lai_auc ~ gwc_auc * Population * `Treatment Name`, 
                 data = lai_gwc_merged)
simple_model <- lm(lai_auc ~ gwc_auc, data = lai_gwc_merged)

# Save model summaries
sink("output/tables/productivity_models_lai.txt")
cat("Full Model ANOVA:\n")
print(anova(full_model))
cat("\nSimple Model Summary:\n")
print(summary(simple_model))
sink()

# Calculate summary statistics
productivity_summary <- lai_gwc_merged %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_lai_auc = mean(lai_auc, na.rm = TRUE),
    se_lai_auc = sd(lai_auc, na.rm = TRUE) / sqrt(n()),
    mean_gwc_auc = mean(gwc_auc, na.rm = TRUE),
    se_gwc_auc = sd(gwc_auc, na.rm = TRUE) / sqrt(n()),
    correlation = cor(lai_auc, gwc_auc, use = "complete.obs"),
    .groups = "drop"
  )

# Save summary statistics
write_csv(productivity_summary, "output/tables/productivity_summary_lai.csv")

# Print some basic diagnostics
cat("\nCorrelation between GWC and LAI by population:\n")
lai_gwc_merged %>%
  group_by(Population) %>%
  summarize(correlation = cor(gwc_auc, lai_auc, use = "complete.obs"))

# Calculate productivity differences between treatments
treatment_comparison <- lai_gwc_merged %>%
  group_by(`Treatment Name`) %>%
  summarize(
    mean_lai_auc = mean(lai_auc, na.rm = TRUE),
    se_lai_auc = sd(lai_auc, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_lai_auc))

print(treatment_comparison)

# Calculate percent reduction from highest to lowest treatment
if(nrow(treatment_comparison) >= 2) {
  highest_productivity <- treatment_comparison$mean_lai_auc[1]
  lowest_productivity <- treatment_comparison$mean_lai_auc[nrow(treatment_comparison)]
  
  percent_reduction <- ((highest_productivity - lowest_productivity) / highest_productivity) * 100
  
  cat("\nProductivity Comparison:\n")
  cat("Highest treatment mean LAI AUC:", round(highest_productivity, 2), "\n")
  cat("Lowest treatment mean LAI AUC:", round(lowest_productivity, 2), "\n")
  cat("Percent reduction:", round(percent_reduction, 2), "%\n")
}

# Also calculate by population and treatment
population_treatment_comparison <- lai_gwc_merged %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_lai_auc = mean(lai_auc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Population) %>%
  mutate(
    max_lai = max(mean_lai_auc),
    percent_of_max = (mean_lai_auc / max_lai) * 100,
    percent_reduction = 100 - percent_of_max
  ) %>%
  arrange(Population, desc(mean_lai_auc))

print("\nProductivity by Population and Treatment:")
print(population_treatment_comparison)
productivity_plot
