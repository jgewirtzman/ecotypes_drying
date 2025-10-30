# 06_aci_curves.R
source("R/00_setup.R")

# Read A/Ci data
aci <- read_csv("data/raw/A_Ci Ouputs_Corrected.csv") %>%
  mutate(
    Population = factor(Population, 
                        levels = POPULATION_LEVELS, 
                        labels = POPULATION_LABELS),
    `Treatment Name` = factor(`Treatment Name`, 
                              levels = TREATMENT_LEVELS)
  )

# Fit mixed models for Vcmax and Jmax
vcmax_model <- lmer(Vcmax ~ Population * `Treatment Name` + (1 | Tag) + (1 | Date), 
                    data = aci)
jmax_model <- lmer(Jmax ~ Population * `Treatment Name` + (1 | Tag) + (1 | Date), 
                   data = aci)

# Get EMMs
vcmax_emm <- emmeans(vcmax_model, ~ Population * `Treatment Name`)
jmax_emm <- emmeans(jmax_model, ~ Population * `Treatment Name`)

# Save model summaries
sink("output/tables/aci_models.txt")
cat("Vcmax Model Summary:\n")
print(summary(vcmax_model))
print(anova(vcmax_model))
cat("\nJmax Model Summary:\n")
print(summary(jmax_model))
print(anova(jmax_model))
sink()

# Function to create point plots
create_point_plot <- function(data, emmeans_df, y_var, y_lab) {
  ggplot(data, aes(x = Population, y = {{y_var}}, fill = `Treatment Name`)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
                size = 2, alpha = 0.6, shape = 21) +
    geom_point(data = emmeans_df, 
               aes(y = emmean, fill = `Treatment Name`),
               position = position_dodge(width = 0.8), 
               size = 3, shape = 21, color = "black") +
    geom_errorbar(data = emmeans_df, 
                  aes(y = emmean, 
                      ymin = emmean - SE, 
                      ymax = emmean + SE,
                      group = `Treatment Name`),
                  position = position_dodge(width = 0.8), 
                  width = 0.2) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    theme_classic() +
    labs(y = y_lab) +
    theme(legend.position = "top")
}

# Create Jmax vs Vcmax scatter plot
scatter_plot <- ggplot(aci, aes(x = Jmax, y = Vcmax)) +
  geom_point(aes(shape = Population, color = `Treatment Name`), 
             size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_shape_manual(values = c(16, 17, 18)) +
  theme_classic() +
  labs(x = expression(J[max]~(mu*mol~m^-2~s^-1)),
       y = expression(V[cmax]~(mu*mol~m^-2~s^-1)),
       shape = "Population") +
  theme(legend.position = "top")

# Create individual plots
vcmax_plot <- create_point_plot(
  aci, 
  as.data.frame(vcmax_emm),
  Vcmax,
  expression(V[cmax]~(mu*mol~m^-2~s^-1))
)

jmax_plot <- create_point_plot(
  aci,
  as.data.frame(jmax_emm),
  Jmax,
  expression(J[max]~(mu*mol~m^-2~s^-1))
)

# Fit and save Jmax-Vcmax relationship
jmax_vcmax_model <- lm(Jmax ~ Vcmax, data = aci)
sink("output/tables/jmax_vcmax_relationship.txt")
print(summary(jmax_vcmax_model))
sink()

library(patchwork)  # last

aci_plot <- (vcmax_plot / jmax_plot | scatter_plot) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "top"))


# Save plot
ggsave("output/figures/aci_combined.pdf", aci_plot, width = 12, height = 8)

# Calculate and save summary statistics
aci_summary <- aci %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    vcmax_mean = mean(Vcmax, na.rm = TRUE),
    vcmax_se = sd(Vcmax, na.rm = TRUE) / sqrt(n()),
    jmax_mean = mean(Jmax, na.rm = TRUE),
    jmax_se = sd(Jmax, na.rm = TRUE) / sqrt(n()),
    j_v_ratio_mean = mean(Jmax/Vcmax, na.rm = TRUE),
    j_v_ratio_se = sd(Jmax/Vcmax, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write_csv(aci_summary, "output/tables/aci_summary.csv")


# Count unique tussocks used for A/Ci measurements
n_tussocks_aci <- aci %>%
  summarize(n_unique_tussocks = n_distinct(Tag))

print(paste("Number of tussocks used for A/Ci curves:", n_tussocks_aci$n_unique_tussocks))

# Count tussocks per population and treatment for A/Ci measurements
aci_distribution <- aci %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(n_tussocks = n_distinct(Tag), .groups = "drop")

print(aci_distribution)


# Count measurements per tussock
measurements_per_tussock <- aci %>%
  group_by(Tag) %>%
  summarize(n_measurements = n_distinct(Date), .groups = "drop")

print(measurements_per_tussock)

# Total number of A/Ci measurements
total_measurements <- aci %>%
  group_by(Tag, Date) %>%
  summarize(n = n(), .groups = "drop") %>%
  nrow()

print(paste("Total measurements:", total_measurements))

# Count unique dates for A/Ci measurements
unique_dates <- aci %>%
  summarize(n_dates = n_distinct(Date))

print(unique_dates)
