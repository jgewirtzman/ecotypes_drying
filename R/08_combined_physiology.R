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

# Remove legend only from water potential (Treatment shown via amax/cond at bottom)
wpot_plot <- wpot_plot + theme(legend.position = "none")

# scatter (panel f): Population shape legend inside the panel; color suppressed
# so Treatment is instead collected from amax/cond at the figure bottom.
scatter_fig5 <- scatter_plot +
  guides(
    color = "none",
    shape = guide_legend(ncol = 1)
  ) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.35, "cm")
  )

# Rebuild aci composite.
# guides = "keep" prevents the outer guides = "collect" from lifting the
# Population shape guide out of panel f.
aci_fig5 <- (vcmax_plot / jmax_plot | scatter_fig5) +
  plot_layout(guides = "keep")
aci_fig5 <- aci_fig5 & guides(fill = "none")

# Outer layout: collect Treatment color from amax/cond at the bottom.
final_plot <- ((wpot_plot | amax_plot) /
                 (cond_plot | aci_fig5)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center"
    )
  ) &
  theme(plot.tag = element_text(face = "plain"))


# Save combined plot
ggsave("output/figures/fig5_combined_physiology.pdf", final_plot, width = 9.5, height = 8, dpi = 600)
ggsave("output/figures/fig5_combined_physiology.png", final_plot, width = 9.5, height = 8, dpi = 300)




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