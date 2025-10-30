library(ggplot2)
library(patchwork)  # <-- needed for `&`


# 08_combined_physiology.R
source("R/00_setup.R")

# Run previous analysis scripts if needed
source("R/05_water_potential.R")  # Creates wpot_plot
source("R/04_gas_exchange.R")     # Creates amax_timepoint and its analyses
source("R/06_aci_curves.R")       # Creates aci objects and models

# Define consistent theme settings
common_theme <- theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  )

# Remove legend for water potential plot
wpot_plot <- wpot_plot + theme(legend.position = "none")

# Modify A/Ci plot to show only shape (population) legend
aci_plot <- aci_plot & guides(color = "none", fill = "none")&theme(legend.position = "bottom")


# Create final combined plot using all the individual plots
final_plot <- ((wpot_plot | amax_plot) /
                 (cond_plot | aci_plot)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center"
    )
  )


# Save combined plot
ggsave("output/figures/combined_physiology.pdf", final_plot, width = 8, height = 8, dpi = 600)




# ============================================================================
# Detailed measurement summary with specific dates and sample sizes
# ============================================================================

cat("\n=== DETAILED PHYSIOLOGICAL MEASUREMENT SUMMARY ===\n\n")

# Water Potential
cat("WATER POTENTIAL:\n")
cat("Unique tussocks:", n_distinct(wpot$Tag), "\n")
cat("Total measurements:", nrow(wpot), "\n")
wpot_dates <- wpot %>% 
  distinct(Date) %>% 
  arrange(Date) %>%
  pull(Date)
cat("Measurement dates:", paste(wpot_dates, collapse = ", "), "\n")
cat("Number of dates:", length(wpot_dates), "\n\n")

# Amax (gas exchange)
cat("AMAX (GAS EXCHANGE):\n")
cat("Unique tussocks:", n_distinct(amax_timepoint$Tag), "\n")
cat("Total measurements:", nrow(amax_timepoint), "\n")
amax_dates <- amax_timepoint %>% 
  distinct(Date) %>% 
  arrange(Date) %>%
  pull(Date)
cat("Measurement dates:", paste(amax_dates, collapse = ", "), "\n")
cat("Number of dates:", length(amax_dates), "\n\n")

# A/Ci curves
cat("A/CI CURVES:\n")
cat("Unique tussocks:", n_distinct(aci$Tag), "\n")
cat("Total measurements:", nrow(aci), "\n")
aci_dates <- aci %>% 
  distinct(Date) %>% 
  arrange(Date) %>%
  pull(Date)
cat("Measurement dates:", paste(aci_dates, collapse = ", "), "\n")
cat("Number of dates:", length(aci_dates), "\n\n")

# Create summary dataframe
measurement_summary <- data.frame(
  Measurement = c("Water Potential", "Amax", "A/Ci Curves"),
  Unique_Tussocks = c(n_distinct(wpot$Tag), n_distinct(amax_timepoint$Tag), n_distinct(aci$Tag)),
  Total_Measurements = c(nrow(wpot), nrow(amax_timepoint), nrow(aci)),
  Number_of_Dates = c(length(wpot_dates), length(amax_dates), length(aci_dates)),
  Dates = c(
    paste(wpot_dates, collapse = ", "),
    paste(amax_dates, collapse = ", "),
    paste(aci_dates, collapse = ", ")
  )
)

print(measurement_summary)
write_csv(measurement_summary, "output/tables/detailed_measurement_summary.csv")