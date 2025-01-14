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
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

# Save combined plot
ggsave("output/figures/combined_physiology.pdf", final_plot, width = 8, height = 8, dpi = 600)
