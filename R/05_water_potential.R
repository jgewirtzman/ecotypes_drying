# 05_water_potential.R
source("R/00_setup.R")

# Read and process water potential data
wpot <- read_csv("data/raw/Ecotypes2017_Drying_Water Potential.csv") %>%
  filter(!is.na(`Pressure (Bar)`)) %>%
  mutate(
    Population = factor(Population, 
                        levels = POPULATION_LEVELS, 
                        labels = POPULATION_LABELS),
    `Treatment Name` = factor(`Treatment Name`, 
                              levels = TREATMENT_LEVELS),
    `Pressure (MPa)` = -`Pressure (Bar)` * 0.1  # Convert from Bar to MPa and flip sign
  )

# Fit mixed model
wpot_model <- lmer(`Pressure (MPa)` ~ Population * `Treatment Name` + (1 | Tag) + (1 | Date), 
                   data = wpot)

# Get EMMs
wpot_emm <- emmeans(wpot_model, ~ Population * `Treatment Name`)
emm_df <- as.data.frame(wpot_emm)

# Get letters for significant differences
wpot_cld <- multcomp::cld(wpot_emm, Letters = letters)
letters_df <- as.data.frame(wpot_cld)

# Save model summary
sink("output/tables/water_potential_model.txt")
cat("Water Potential Model Summary:\n")
print(summary(wpot_model))
print(anova(wpot_model))
cat("\nEstimated Marginal Means:\n")
print(wpot_emm)
cat("\nPairwise Comparisons:\n")
print(wpot_cld)
sink()

letters_df <- letters_df %>%
  left_join(wpot %>%
              group_by(Population, `Treatment Name`) %>%
              summarize(min_y = min(`Pressure (MPa)`, na.rm = TRUE), .groups = "drop"),
            by = c("Population", "Treatment Name")) %>%
  mutate(y_position = min_y - 0.1 * abs(min_y))

wpot_plot <- ggplot(wpot, aes(x = Population, y = `Pressure (MPa)`, 
                              color = `Treatment Name`, fill = `Treatment Name`)) +
  geom_violin(alpha = 0.2, position = position_dodge(width = 0.8), color = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 3, alpha = 0.5) +
  geom_point(data = emm_df,
             aes(y = emmean, group = `Treatment Name`),
             size = 5, shape = 21, color = "black",
             position = position_dodge(width = 0.8)) +
  geom_errorbar(data = emm_df,
                aes(y = emmean, 
                    ymin = emmean - SE, 
                    ymax = emmean + SE,
                    group = `Treatment Name`),
                position = position_dodge(width = 0.8),
                width = 0.2, color = "black") +
  geom_text(data = letters_df,
            aes(y = y_position, label = .group),
            position = position_dodge(width = 0.8),
            size = 5, color = "black") +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  labs(y = "Water potential (MPa)") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_y_reverse()



# Save plot
ggsave("output/figures/water_potential.pdf", wpot_plot, width = 8, height = 6)

# Calculate and save summary statistics
wpot_summary <- wpot %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_pressure = mean(`Pressure (MPa)`, na.rm = TRUE),
    se_pressure = sd(`Pressure (MPa)`, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write_csv(wpot_summary, "output/tables/water_potential_summary.csv")
