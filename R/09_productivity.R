# 09_productivity.R
source("R/00_setup.R")

# Read and process NDVI data
ndvi <- read_csv("data/raw/Ecotypes2017_Drying_NDVI.csv") %>%
  # Remove columns with all NA values
  select_if(~ !all(is.na(.))) %>%
  # Remove rows with all NA values
  filter(rowSums(is.na(.)) < ncol(.))

# Calculate NDVI AUC
ndvi_dat <- ndvi %>%
  group_by(Population, `Treatment Name`, Tag) %>%
  summarize(
    ndvi_auc = pracma::trapz(DOY[!is.na(NDVI)], NDVI[!is.na(NDVI)]),
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

# Merge NDVI and GWC AUC datasets
ndvi_gwc_merged <- left_join(ndvi_dat, gwc_dat, 
                             by = c("Population", "Treatment Name", "Tag"))

# Update population names and order
ndvi_gwc_merged <- ndvi_gwc_merged %>%
  mutate(
    Population = factor(Population,
                        levels = POPULATION_LEVELS,
                        labels = POPULATION_LABELS),
    `Treatment Name` = factor(`Treatment Name`, levels = TREATMENT_LEVELS)
  )

# Calculate z-scores
ndvi_gwc_merged <- ndvi_gwc_merged %>%
  mutate(
    gwc_auc_z = scale(gwc_auc),
    ndvi_auc_z = scale(ndvi_auc)
  )

# Custom colors
custom_colors <- c(
  "Coldfoot" = "#3E6347",  
  "Toolik" = "#90C49B", 
  "Sagwon" = "#79A9C8"
)

# Create productivity plot
productivity_plot <- ggplot(ndvi_gwc_merged, 
                            aes(y = ndvi_auc_z, x = gwc_auc_z, 
                                color = Population, fill = Population)) +
  geom_point(size = 3.5, shape = 21, stroke = 1.2, alpha = 0.7) +
  geom_smooth(method = "lm", size = 1.5, 
              se = TRUE, fullrange = TRUE, alpha = 0.2) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  labs(
    x = "Growing Season Soil Moisture\n(GWC AUC Z-Score)",
    y = "Cumulative Vegetation Greenness\n(NDVI AUC Z-Score)",
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
ggsave("output/figures/productivity_relationship.pdf", 
       productivity_plot, width = 6, height = 6, dpi=600)

# Fit models and save statistics
full_model <- lm(ndvi_auc ~ gwc_auc * Population * `Treatment Name`, 
                 data = ndvi_gwc_merged)
simple_model <- lm(ndvi_auc ~ gwc_auc, data = ndvi_gwc_merged)

# Save model summaries
sink("output/tables/productivity_models.txt")
cat("Full Model ANOVA:\n")
print(anova(full_model))
cat("\nSimple Model Summary:\n")
print(summary(simple_model))
sink()

# Calculate summary statistics
productivity_summary <- ndvi_gwc_merged %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_ndvi_auc = mean(ndvi_auc, na.rm = TRUE),
    se_ndvi_auc = sd(ndvi_auc, na.rm = TRUE) / sqrt(n()),
    mean_gwc_auc = mean(gwc_auc, na.rm = TRUE),
    se_gwc_auc = sd(gwc_auc, na.rm = TRUE) / sqrt(n()),
    correlation = cor(ndvi_auc, gwc_auc, use = "complete.obs"),
    .groups = "drop"
  )

# Save summary statistics
write_csv(productivity_summary, "output/tables/productivity_summary.csv")

# Print some basic diagnostics
cat("\nCorrelation between GWC and NDVI by population:\n")
ndvi_gwc_merged %>%
  group_by(Population) %>%
  summarize(correlation = cor(gwc_auc, ndvi_auc, use = "complete.obs"))





# Calculate productivity differences between treatments
treatment_comparison <- ndvi_gwc_merged %>%
  group_by(`Treatment Name`) %>%
  summarize(
    mean_ndvi_auc = mean(ndvi_auc, na.rm = TRUE),
    se_ndvi_auc = sd(ndvi_auc, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_ndvi_auc))

print(treatment_comparison)

# Calculate percent reduction from highest to lowest treatment
if(nrow(treatment_comparison) >= 2) {
  highest_productivity <- treatment_comparison$mean_ndvi_auc[1]
  lowest_productivity <- treatment_comparison$mean_ndvi_auc[nrow(treatment_comparison)]
  
  percent_reduction <- ((highest_productivity - lowest_productivity) / highest_productivity) * 100
  
  cat("\nProductivity Comparison:\n")
  cat("Highest treatment mean NDVI AUC:", round(highest_productivity, 2), "\n")
  cat("Lowest treatment mean NDVI AUC:", round(lowest_productivity, 2), "\n")
  cat("Percent reduction:", round(percent_reduction, 2), "%\n")
}

# Also calculate by population and treatment
population_treatment_comparison <- ndvi_gwc_merged %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_ndvi_auc = mean(ndvi_auc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Population) %>%
  mutate(
    max_ndvi = max(mean_ndvi_auc),
    percent_of_max = (mean_ndvi_auc / max_ndvi) * 100,
    percent_reduction = 100 - percent_of_max
  ) %>%
  arrange(Population, desc(mean_ndvi_auc))

print("\nProductivity by Population and Treatment:")
print(population_treatment_comparison)

# Save results
write_csv(treatment_comparison, "output/tables/treatment_productivity_comparison.csv")
write_csv(population_treatment_comparison, "output/tables/population_treatment_productivity.csv")
