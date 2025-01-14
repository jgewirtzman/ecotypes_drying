# 04_gas_exchange.R
source("R/00_setup.R")

# Read Amax data
amax <- read_csv("data/raw/Ecotypes2017_Drying_Amax.csv")

# Calculate mean values for each measurement timepoint 
amax_timepoint <- amax %>%
  group_by(Population, `Treatment Name`, Tag, Date, DOY) %>%
  summarize(
    Amax_avg = mean(Photo, na.rm = TRUE),
    Cond_avg = mean(Cond, na.rm = TRUE),
    Ci_avg = mean(Ci, na.rm = TRUE),
    Trans_avg = mean(Trmmol, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Clean factor levels
  mutate(
    Population = factor(Population, 
                        levels = POPULATION_LEVELS, 
                        labels = POPULATION_LABELS),
    `Treatment Name` = factor(`Treatment Name`, 
                              levels = TREATMENT_LEVELS)
  )

# Fit mixed models
amax_model <- lmer(Amax_avg ~ Population * `Treatment Name` + (1 | Tag) + (1 | Date), 
                   data = amax_timepoint)
cond_model <- lmer(Cond_avg ~ Population * `Treatment Name` + (1 | Tag) + (1 | Date), 
                   data = amax_timepoint)
ci_model <- lmer(Ci_avg ~ Population * `Treatment Name` + (1 | Tag) + (1 | Date), 
                 data = amax_timepoint)

# Save model summaries
sink("output/tables/gas_exchange_models.txt")
cat("Amax Model Summary:\n")
print(summary(amax_model))
print(anova(amax_model))

cat("\nConductance Model Summary:\n")
print(summary(cond_model))
print(anova(cond_model))

cat("\nCi Model Summary:\n")
print(summary(ci_model))
print(anova(ci_model))
sink()

#function to create plot
create_phys_plot <- function(data, response_var, y_label) {
  # Fit model and calculate adjusted means
  model <- lmer(as.formula(paste(response_var, "~ Population * `Treatment Name` + (1 | Tag) + (1 | Date)")), 
                data = data)
  adjusted_means <- emmeans(model, ~ Population * `Treatment Name`)
  means_df <- as.data.frame(adjusted_means)
  
  # Get Tukey letters
  tukey_letters <- multcomp::cld(adjusted_means, Letters = letters)
  tukey_letters_df <- as.data.frame(tukey_letters)
  
  # Calculate dynamic y_position
  max_values <- data %>%
    group_by(Population, `Treatment Name`) %>%
    summarize(max_y = max(.data[[response_var]], na.rm = TRUE), .groups = "drop")
  
  tukey_letters_df <- tukey_letters_df %>%
    left_join(max_values, by = c("Population", "Treatment Name")) %>%
    mutate(y_position = max_y + 0.1 * abs(max_y))
  
  # Create plot
  ggplot(data, aes(x = Population, y = .data[[response_var]], fill = `Treatment Name`)) +
    geom_violin(alpha = 0.2, position = position_dodge(width = 0.8), scale = "width", color = NA) +
    geom_point(aes(color = `Treatment Name`),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               size = 3, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point",
                 aes(group = `Treatment Name`),
                 position = position_dodge(width = 0.8),
                 size = 5, shape = 21, color = "black") +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar",
                 aes(group = `Treatment Name`),
                 position = position_dodge(width = 0.8),
                 width = 0.2, color = "black") +
    geom_text(data = tukey_letters_df,
              aes(y = y_position, label = .group),
              position = position_dodge(width = 0.8),
              size = 5, color = "black") +
    scale_color_manual(values = TREATMENT_COLORS) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    labs(y = y_label) +
    theme_classic() +
    theme(legend.position = "top")
}


# Create plots
amax_plot <- create_phys_plot(amax_timepoint, "Amax_avg", 
                              expression(A[max]~"("*mu*mol~CO[2]~m^-2~s^-1*")"))
cond_plot <- create_phys_plot(amax_timepoint, "Cond_avg",
                              expression(g[s]~"("*mol~H[2]*O~m^-2~s^-1*")"))
ci_plot <- create_phys_plot(amax_timepoint, "Ci_avg",
                            expression(C[i]~"("*mu*mol~CO[2]~mol^-1*")"))

# Combine plots
final_plot <- (amax_plot / cond_plot / ci_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# Save plots
ggsave("output/figures/gas_exchange_combined.pdf", final_plot, width = 10, height = 12)

# Calculate and save summary statistics
gas_exchange_summary <- amax_timepoint %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    across(ends_with("_avg"),
           list(
             mean = ~mean(., na.rm = TRUE),
             se = ~sd(., na.rm = TRUE) / sqrt(n())
           )
    ),
    .groups = "drop"
  )

write_csv(gas_exchange_summary, "output/tables/gas_exchange_summary.csv")