# 02_phenology.R

# ----------------------------------
# Load setup and utilities
# ----------------------------------
source("R/00_setup.R")

# ----------------------------------
# Step 1: Read and Process Data
# ----------------------------------
pheno <- read_csv("data/raw/Ecotypes2017_Drying_Phenology.csv") %>%
  filter(
    !is.na(`Tiller Total Green Length`) |
      !is.na(`Tiller % Green`) |
      !is.na(`Average Green Leaf Length`)
  ) %>%
  rename(
    TillerTotalGreenLength = `Tiller Total Green Length`,
    TillerPercentGreen = `Tiller % Green`,
    AverageGreenLeafLength = `Average Green Leaf Length`
  ) %>%
  clean_factor_levels()

# ----------------------------------
# Step 2: Plot Phenology Metrics Over Time
# ----------------------------------
plot_pheno_metric <- function(data, metric, y_label) {
  ggplot(data, aes(x = DOY, y = !!sym(metric), color = TreatmentName)) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.3) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    scale_color_manual(values = TREATMENT_COLORS) +
    labs(x = "Day of Year", y = y_label, color = "Treatment") +
    theme_classic() +
    facet_wrap(~ Population, ncol = 3) +
    theme(legend.position = "top")
}

# Generate time series plots
tot_gl_plot <- plot_pheno_metric(pheno, "TillerTotalGreenLength", "Total Green Length (cm)")
pct_green_plot <- plot_pheno_metric(pheno, "TillerPercentGreen", "Percent Green (%)")
avg_green_plot <- plot_pheno_metric(pheno, "AverageGreenLeafLength", "Mean Green Length (cm)")

# ----------------------------------
# Step 3: Fit Mixed-Effects Models
# ----------------------------------
tiller_model <- lmer(TillerTotalGreenLength ~ `Treatment Name` * Population * poly(DOY, 2) + (1 | Tag), data = pheno)
percent_green_model <- lmer(TillerPercentGreen ~ `Treatment Name` * Population * poly(DOY, 2) + (1 | Tag), data = pheno)
avg_green_length_model <- lmer(AverageGreenLeafLength ~ `Treatment Name` * Population * poly(DOY, 2) + (1 | Tag), data = pheno)

# ----------------------------------
# Step 4: Save Model Summaries
# ----------------------------------
sink("output/tables/phenology_model_summaries.txt")

cat("Tiller Total Green Length Model:\n")
print(summary(tiller_model))
print(anova(tiller_model))

cat("\nPercent Green Model:\n")
print(summary(percent_green_model))
print(anova(percent_green_model))

cat("\nAverage Green Leaf Length Model:\n")
print(summary(avg_green_length_model))
print(anova(avg_green_length_model))

sink()

# ----------------------------------
# Step 5: Random Forest Models
# ----------------------------------
rf_data_total <- pheno %>%
  dplyr::select(TillerTotalGreenLength, Population, TreatmentName, DOY) %>%
  drop_na()

rf_data_percent <- pheno %>%
  dplyr::select(TillerPercentGreen, Population, TreatmentName, DOY) %>%
  drop_na()

rf_data_avg <- pheno %>%
  dplyr::select(AverageGreenLeafLength, Population, TreatmentName, DOY) %>%
  drop_na()

# Function to save RF model summary to a text file
save_rf_summary <- function(rf_model, file_name) {
  sink(file_name)
  cat("Random Forest Model Summary:\n")
  print(rf_model)
  cat("\nVariable Importance:\n")
  print(importance(rf_model))
  sink()
}

# Save RF model summaries to text files
set.seed(123)
rf_total <- randomForest(TillerTotalGreenLength ~ ., data = rf_data_total, importance = TRUE, ntree = 10000)
save_rf_summary(rf_total, "output/tables/pheno_rf_total_summary.txt")

rf_percent <- randomForest(TillerPercentGreen ~ ., data = rf_data_percent, importance = TRUE, ntree = 10000)
save_rf_summary(rf_percent, "output/tables/pheno_rf_percent_summary.txt")

rf_avg <- randomForest(AverageGreenLeafLength ~ ., data = rf_data_avg, importance = TRUE, ntree = 10000)
save_rf_summary(rf_avg, "output/tables/pheno_rf_avg_summary.txt")

write_csv(as.data.frame(importance(rf_total)), "output/tables/pheno_rf_total_importance.csv")
write_csv(as.data.frame(importance(rf_percent)), "output/tables/pheno_rf_percent_importance.csv")
write_csv(as.data.frame(importance(rf_avg)), "output/tables/pheno_rf_avg_importance.csv")

# ----------------------------------
# Step 6: Create Variable Importance Plots
# ----------------------------------
varImpPlot_to_gg <- function(rf_model) {
  imp <- as.data.frame(importance(rf_model))
  imp$Variable <- rownames(imp)
  imp$Variable <- gsub("TreatmentName", "Treatment", imp$Variable)
  
  ggplot(imp, aes(x = reorder(Variable, `%IncMSE`), y = `%IncMSE`)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(y = "%IncMSE", x = NULL) +
    theme_classic()
}

imp_total_plot <- varImpPlot_to_gg(rf_total)
imp_percent_plot <- varImpPlot_to_gg(rf_percent)
imp_avg_plot <- varImpPlot_to_gg(rf_avg)

# ----------------------------------
# Step 7: Combine Time Series Plots and Importance Plots
# ----------------------------------
final_plot <- wrap_plots(
  tot_gl_plot, imp_total_plot,
  pct_green_plot, imp_percent_plot,
  avg_green_plot, imp_avg_plot,
  ncol = 2,
  widths = c(4, 1),
  guides = "collect"
) &
  theme(legend.position = "top")

ggsave("output/figures/phenology_combined.pdf", final_plot, width = 8, height = 7, units = "in", dpi = 600)

# ----------------------------------
# Step 8: Save Fixed Effects Summary as a Polished Table
# ----------------------------------

# Extract fixed effects summaries
fixed_effects_summary <- bind_rows(
  broom.mixed::tidy(tiller_model) %>% mutate(Metric = "Tiller Total Green Length"),
  broom.mixed::tidy(percent_green_model) %>% mutate(Metric = "Tiller Percent Green"),
  broom.mixed::tidy(avg_green_length_model) %>% mutate(Metric = "Average Green Leaf Length")
)

# Save fixed effects as a .txt file
sink("output/tables/phenology_fixed_effects_summary.txt")
cat("Fixed Effects Summary for Phenology Models:\n")
print(fixed_effects_summary)
sink()

# Save a polished version using gt
fixed_effects_summary %>%
  dplyr::select(Metric, term, estimate, std.error, df, statistic, p.value) %>%
  arrange(Metric, term) %>%
  gt() %>%
  tab_header(title = "Fixed Effects Summary for Phenology Models") %>%
  fmt_number(columns = c(estimate, std.error, statistic, p.value), decimals = 3) %>%
  gtsave("output/tables/phenology_fixed_effects_summary.rtf")

fixed_effects_summary %>%
  dplyr::select(Metric, term, estimate, std.error, df, statistic, p.value) %>%
  arrange(Metric, term) %>%
  gt() %>%
  tab_header(title = "Fixed Effects Summary for Phenology Models") %>%
  fmt_number(columns = c(estimate, std.error, statistic, p.value), decimals = 3) %>%
  gtsave("output/tables/phenology_fixed_effects_summary.html")
