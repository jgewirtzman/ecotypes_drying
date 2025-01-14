# 03_canopy.R
source("R/00_setup.R")

# ----------------------------------
# Step 1: Data Preparation
# ----------------------------------

# Read NDVI data
ndvi <- read_csv("data/raw/Ecotypes2017_Drying_NDVI.csv") %>%
  # Remove columns with all NA values
  dplyr::select_if(~ !all(is.na(.))) %>%
  # Remove rows with all NA values
  filter(rowSums(is.na(.)) < ncol(.)) %>%
  # Clean factor levels using the utility function
  clean_factor_levels() %>%
  mutate(
    Treatment = factor(`Treatment Name`, levels = c("Wet", "Dry", "Deep"))
  )

# Calculate LAI and Biomass
ndvi <- ndvi %>%
  mutate(
    LAI = 0.03 * exp(7.65 * NDVI),
    Biomass = 0.0256 * exp(5.32 * NDVI)
  )

# ----------------------------------
# Step 2: AUC Calculation
# ----------------------------------

calculate_auc <- function(data, metric) {
  data %>%
    group_by(Tag, Population, Treatment) %>%
    summarize(AUC = pracma::trapz(DOY, .data[[metric]]), .groups = "drop")
}

lai_auc_data <- calculate_auc(ndvi, "LAI")
ndvi_auc_data <- calculate_auc(ndvi, "NDVI")
biomass_auc_data <- calculate_auc(ndvi, "Biomass")

# ----------------------------------
# Step 3: Fit Quadratic Mixed Models
# ----------------------------------

fit_quadratic_model <- function(data, response_var) {
  formula <- as.formula(paste(response_var, "~ poly(DOY, 2) * Treatment * Population + (1 | Tag)"))
  model <- lmer(formula, data = data)
  return(model)
}

lai_model <- fit_quadratic_model(ndvi, "LAI")
ndvi_model <- fit_quadratic_model(ndvi, "NDVI")
biomass_model <- fit_quadratic_model(ndvi, "Biomass")

# ----------------------------------
# Step 4: ANOVA and Compact Letter Display (CLD)
# ----------------------------------

anova_cld_results <- function(model, response_var, auc_data) {
  cat("\nANOVA Results for", response_var, ":\n")
  print(anova(model))
  
  emm <- emmeans(model, ~ Treatment * Population)
  cld <- multcomp::cld(emm, Letters = letters)
  
  max_auc <- auc_data %>%
    group_by(Treatment, Population) %>%
    summarize(max_y = max(AUC, na.rm = TRUE), .groups = "drop")
  
  cld_data <- as.data.frame(cld) %>%
    left_join(max_auc, by = c("Treatment", "Population")) %>%
    mutate(y_position = max_y)
  
  return(cld_data)
}

lai_cld <- anova_cld_results(lai_model, "LAI", lai_auc_data)
ndvi_cld <- anova_cld_results(ndvi_model, "NDVI", ndvi_auc_data)
biomass_cld <- anova_cld_results(biomass_model, "Biomass", biomass_auc_data)

# ----------------------------------
# Step 5: Generate Violin Plots
# ----------------------------------

plot_auc_violin <- function(auc_data, cld_data, y_label) {
  # Define treatment colors
  treatment_colors <- c("Wet" = "#1f77b4", "Dry" = "#8c564b", "Deep" = "#d62728")
  
  # Plot with the exact aesthetics you requested
  ggplot(auc_data, aes(x = Treatment, y = AUC, fill = Treatment)) +
    geom_violin(
      alpha = 0.2, position = position_dodge(width = 0.8), scale = "width", color = NA
    ) +
    geom_point(
      aes(color = Treatment, fill = Treatment), 
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
      size = 3, shape = 16, alpha = 0.5
    ) +
    stat_summary(
      fun = mean, geom = "point", 
      position = position_dodge(width = 0.8), size = 5, shape = 21, 
      aes(fill = Treatment), color = "black"
    ) +
    stat_summary(
      fun.data = mean_cl_boot, geom = "errorbar", 
      position = position_dodge(width = 0.8), width = 0.2, color = "black"
    ) +
    geom_text(
      data = cld_data, aes(x = Treatment, y = y_position, label = .group), 
      size = 5, vjust = -0.5, hjust = 0.5, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = treatment_colors) +
    scale_color_manual(values = treatment_colors) +
    labs(x = "Treatment", y = y_label) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.title = element_blank()
    ) +
    expand_limits(y = max(auc_data$AUC, na.rm = TRUE) * 1.2) +
    facet_wrap(~ Population, nrow = 3)
}

lai_violin <- plot_auc_violin(lai_auc_data, lai_cld, "LAI AUC")
ndvi_violin <- plot_auc_violin(ndvi_auc_data, ndvi_cld, "NDVI AUC")
biomass_violin <- plot_auc_violin(biomass_auc_data, biomass_cld, "Biomass AUC")

# ----------------------------------
# Step 6: Generate Time-Series Plots
# ----------------------------------

plot_time_series <- function(data, y_var, y_label) {
  ggplot(data, aes(x = DOY, y = !!sym(y_var), color = Treatment)) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
    stat_summary(fun = mean, geom = "line", size = 1.2) +
    scale_color_manual(values = TREATMENT_COLORS) +
    labs(x = "Day of Year", y = y_label) +
    theme_classic() +
    facet_wrap(~ Population, ncol = 1) +
    theme(legend.position = "top")
}

lai_time_series <- plot_time_series(ndvi, "LAI", "Leaf Area Index")
ndvi_time_series <- plot_time_series(ndvi, "NDVI", "NDVI")
biomass_time_series <- plot_time_series(ndvi, "Biomass", "Biomass")

# ----------------------------------
# Step 7: Random Forest Importance Plots
# ----------------------------------

fit_rf <- function(data, metric) {
  rf_data <- data %>%
    dplyr::select(all_of(c(metric, "Population", "Treatment", "DOY"))) %>%
    drop_na()
  
  randomForest(as.formula(paste(metric, "~ .")), data = rf_data, importance = TRUE, ntree = 10000)
}

lai_rf <- fit_rf(ndvi, "LAI")
ndvi_rf <- fit_rf(ndvi, "NDVI")
biomass_rf <- fit_rf(ndvi, "Biomass")

plot_rf_importance <- function(rf_model, metric) {
  imp <- as.data.frame(randomForest::importance(rf_model))
  imp$Variable <- rownames(imp)
  
  ggplot(imp, aes(x = reorder(Variable, `%IncMSE`), y = `%IncMSE`)) +
    geom_col() +
    coord_flip() +
    labs(y = "%IncMSE", x = NULL, title = paste("Variable Importance:", metric)) +
    theme_classic()
}

lai_rf_plot <- plot_rf_importance(lai_rf, "LAI")
ndvi_rf_plot <- plot_rf_importance(ndvi_rf, "NDVI")
biomass_rf_plot <- plot_rf_importance(biomass_rf, "Biomass")

# Function to save RF model summaries to text files
save_rf_summary <- function(rf_model, file_name) {
  sink(file_name)
  cat("Random Forest Model Summary:\n")
  print(rf_model)
  cat("\nVariable Importance:\n")
  print(importance(rf_model))
  sink()
}

# Save RF summaries for LAI, NDVI, and Biomass
save_rf_summary(lai_rf, "output/tables/canopy_rf_lai_summary.txt")
save_rf_summary(ndvi_rf, "output/tables/canopy_rf_ndvi_summary.txt")
save_rf_summary(biomass_rf, "output/tables/canopy_rf_biomass_summary.txt")


# Save variable importance metrics as CSV files
write_csv(as.data.frame(importance(lai_rf)), "output/tables/canopy_rf_lai_importance.csv")
write_csv(as.data.frame(importance(ndvi_rf)), "output/tables/canopy_rf_ndvi_importance.csv")
write_csv(as.data.frame(importance(biomass_rf)), "output/tables/canopy_rf_biomass_importance.csv")



# ----------------------------------
# Step 8: Combine and Save Plots
# ----------------------------------

plot_metrics <- function(time_series_plot, violin_plot, importance_plot) {
  # First combine the top plots with 2:1 width ratio
  top_row <- time_series_plot + violin_plot + plot_layout(widths = c(2, 1))
  
  # Then combine top and bottom with the height ratio (4:1)
  combined <- top_row / importance_plot + plot_layout(heights = c(4, 1))
  
  return(combined)
}


# Generate final combined plots
final_lai_plot <- plot_metrics(lai_time_series, lai_violin, lai_rf_plot)
final_ndvi_plot <- plot_metrics(ndvi_time_series, ndvi_violin, ndvi_rf_plot)
final_biomass_plot <- plot_metrics(biomass_time_series, biomass_violin, biomass_rf_plot)

# Save the plots
ggsave("output/figures/canopy_lai.pdf", final_lai_plot, width = 8, height = 8, units = "in", dpi = 600)
ggsave("output/figures/canopy_ndvi.pdf", final_ndvi_plot, width = 10, height = 12, units = "in", dpi = 600)
ggsave("output/figures/canopy_biomass.pdf", final_biomass_plot, width = 10, height = 12, units = "in", dpi = 600)



# Calculate and save AUC values
auc_values <- ndvi %>%
  group_by(Population, Treatment, Tag) %>%
  summarise(
    ndvi_auc = pracma::trapz(DOY[!is.na(NDVI)], NDVI[!is.na(NDVI)]),
    lai_auc = pracma::trapz(DOY[!is.na(LAI)], LAI[!is.na(LAI)]),
    biomass_auc = pracma::trapz(DOY[!is.na(Biomass)], Biomass[!is.na(Biomass)]),
    .groups = "drop"
  )

write_csv(auc_values, "output/tables/canopy_auc_values.csv")


# ----------------------------------
# Step 9: Create and Save GT Table
# ----------------------------------

# Create fixed effects summary
fixed_effects_summary <- bind_rows(
  broom.mixed::tidy(lai_model) %>% mutate(Metric = "LAI"),
  broom.mixed::tidy(ndvi_model) %>% mutate(Metric = "NDVI"),
  broom.mixed::tidy(biomass_model) %>% mutate(Metric = "Biomass")
)

# Create a polished table using gt
fixed_effects_summary %>%
  dplyr::select(Metric, term, estimate, std.error, df, statistic, p.value) %>%
  arrange(Metric, term) %>%
  gt() %>%
  tab_header(title = "Fixed Effects Summary for Canopy Models") %>%
  fmt_number(columns = c(estimate, std.error, statistic, p.value), decimals = 3) %>%
  gtsave("output/tables/canopy_fixed_effects_summary.rtf")  # Save as PDF


# ----------------------------------
# Step 10: Save ANOVA Summary as a Polished Table
# ----------------------------------

# Perform Type III ANOVA on each model
lai_anova <- anova(lai_model, type = "III")
ndvi_anova <- anova(ndvi_model, type = "III")
biomass_anova <- anova(biomass_model, type = "III")

# Convert ANOVA results to data frames for gt
anova_summary <- bind_rows(
  as.data.frame(lai_anova) %>% mutate(Metric = "LAI"),
  as.data.frame(ndvi_anova) %>% mutate(Metric = "NDVI"),
  as.data.frame(biomass_anova) %>% mutate(Metric = "Biomass")
) %>%
  rownames_to_column(var = "Term") %>%
  dplyr::select(Metric, Term, `Sum Sq`, `Mean Sq`, NumDF, DenDF, `F value`, `Pr(>F)`)
  
# Save the ANOVA summary as a polished table using gt
anova_summary %>%
  gt() %>%
  tab_header(title = "Type III ANOVA Summary for Canopy Models") %>%
  fmt_number(columns = c(`Sum Sq`, `Mean Sq`, `F value`, `Pr(>F)`), decimals = 3) %>%
  cols_label(
    Term = "Effect",
    `Sum Sq` = "Sum of Squares",
    `Mean Sq` = "Mean Squares",
    `NumDF` = "Numerator DF",
    `DenDF` = "Denominator DF",
    `F value` = "F-Value",
    `Pr(>F)` = "p-Value"
  ) %>%
  gtsave("output/tables/canopy_anova_summary.rtf")

# ----------------------------------
# Step 11: Save Model Summaries as TXT
# ----------------------------------

# Open a connection to a text file
sink("output/tables/canopy_model_summaries.txt")

# Print summaries for each model
cat("LAI Model:\n")
print(summary(lai_model))
cat("\n")
cat("NDVI Model:\n")
print(summary(ndvi_model))
cat("\n")
cat("Biomass Model:\n")
print(summary(biomass_model))
cat("\n")

# Perform Type III ANOVA and print results
cat("Type III ANOVA for LAI Model:\n")
print(anova(lai_model, type = "III"))
cat("\n")
cat("Type III ANOVA for NDVI Model:\n")
print(anova(ndvi_model, type = "III"))
cat("\n")
cat("Type III ANOVA for Biomass Model:\n")
print(anova(biomass_model, type = "III"))
cat("\n")

# Close the text file connection
sink()

