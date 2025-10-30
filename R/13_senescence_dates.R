# ----------------------------------
# Calculate Senescence Onset from NDVI Data
# ----------------------------------

# Read NDVI data (if not already loaded)
if(!exists("ndvi")) {
  ndvi <- read_csv("data/raw/Ecotypes2017_Drying_NDVI.csv") %>%
    dplyr::select_if(~ !all(is.na(.))) %>%
    filter(rowSums(is.na(.)) < ncol(.)) %>%
    clean_factor_levels() %>%
    mutate(
      Treatment = factor(`Treatment Name`, levels = c("Wet", "Dry", "Deep"))
    )
}

# Method: Find peak NDVI for each plant, then first decline below threshold
senescence_after_peak_ndvi <- function(data, threshold_pct_of_peak = 0.95) {
  data %>%
    arrange(Tag, DOY) %>%
    group_by(Tag, Population, Treatment) %>%
    filter(!is.na(NDVI)) %>%
    mutate(
      max_ndvi = max(NDVI, na.rm = TRUE),
      DOY_at_peak = DOY[which.max(NDVI)][1],
      after_peak = DOY > DOY_at_peak,
      threshold_value = max_ndvi * threshold_pct_of_peak,
      below_threshold = NDVI < threshold_value
    ) %>%
    filter(after_peak, below_threshold) %>%
    slice(1) %>%
    dplyr::select(
      Tag, Population, Treatment, 
      DOY_senescence_onset = DOY,
      DOY_at_peak,
      max_ndvi,
      NDVI_at_onset = NDVI
    ) %>%
    ungroup()
}

# Calculate senescence onset at different thresholds
senescence_ndvi_95pct <- senescence_after_peak_ndvi(ndvi, threshold_pct_of_peak = 0.95)
senescence_ndvi_90pct <- senescence_after_peak_ndvi(ndvi, threshold_pct_of_peak = 0.90)
senescence_ndvi_85pct <- senescence_after_peak_ndvi(ndvi, threshold_pct_of_peak = 0.85)

# ----------------------------------
# Compare Thresholds
# ----------------------------------

cat("\n=== COMPARING NDVI THRESHOLD DEFINITIONS ===\n\n")

# Summarize each threshold
summary_95 <- senescence_ndvi_95pct %>%
  group_by(Population, Treatment) %>%
  summarise(
    mean_onset_95 = mean(DOY_senescence_onset, na.rm = TRUE),
    n_95 = n(),
    .groups = "drop"
  )

summary_90 <- senescence_ndvi_90pct %>%
  group_by(Population, Treatment) %>%
  summarise(
    mean_onset_90 = mean(DOY_senescence_onset, na.rm = TRUE),
    n_90 = n(),
    .groups = "drop"
  )

summary_85 <- senescence_ndvi_85pct %>%
  group_by(Population, Treatment) %>%
  summarise(
    mean_onset_85 = mean(DOY_senescence_onset, na.rm = TRUE),
    n_85 = n(),
    .groups = "drop"
  )

# Combine all thresholds
threshold_comparison_ndvi <- summary_95 %>%
  left_join(summary_90, by = c("Population", "Treatment")) %>%
  left_join(summary_85, by = c("Population", "Treatment")) %>%
  mutate(
    diff_95_to_90 = mean_onset_90 - mean_onset_95,
    diff_90_to_85 = mean_onset_85 - mean_onset_90
  )

cat("Mean senescence onset DOY by threshold:\n\n")
print(threshold_comparison_ndvi)

cat("\n\n=== OVERALL MEANS ACROSS ALL GROUPS ===\n\n")

overall_means_ndvi <- data.frame(
  Threshold = c("95%", "90%", "85%"),
  Mean_DOY = c(
    mean(senescence_ndvi_95pct$DOY_senescence_onset, na.rm = TRUE),
    mean(senescence_ndvi_90pct$DOY_senescence_onset, na.rm = TRUE),
    mean(senescence_ndvi_85pct$DOY_senescence_onset, na.rm = TRUE)
  ),
  Total_N = c(
    nrow(senescence_ndvi_95pct),
    nrow(senescence_ndvi_90pct),
    nrow(senescence_ndvi_85pct)
  )
)

print(overall_means_ndvi)

cat("\n=== DIFFERENCE IN DAYS BETWEEN THRESHOLDS ===\n")
cat(sprintf("95%% to 90%% threshold: %.1f days later\n", 
            overall_means_ndvi$Mean_DOY[2] - overall_means_ndvi$Mean_DOY[1]))
cat(sprintf("90%% to 85%% threshold: %.1f days later\n", 
            overall_means_ndvi$Mean_DOY[3] - overall_means_ndvi$Mean_DOY[2]))

# ----------------------------------
# Visualize Thresholds
# ----------------------------------

set.seed(456)
example_tags <- sample(unique(ndvi$Tag), 9)
example_ndvi <- ndvi %>% filter(Tag %in% example_tags)

threshold_comparison_long <- bind_rows(
  senescence_ndvi_95pct %>% mutate(threshold = "95%"),
  senescence_ndvi_90pct %>% mutate(threshold = "90%"),
  senescence_ndvi_85pct %>% mutate(threshold = "85%")
) %>%
  filter(Tag %in% example_tags)

threshold_vis_ndvi <- ggplot(example_ndvi, aes(x = DOY, y = NDVI)) +
  geom_line(aes(group = Tag)) +
  geom_point() +
  geom_vline(data = threshold_comparison_long %>% distinct(Tag, DOY_at_peak), 
             aes(xintercept = DOY_at_peak), 
             color = "green", linetype = "dashed", linewidth = 0.5) +
  geom_vline(data = threshold_comparison_long, 
             aes(xintercept = DOY_senescence_onset, color = threshold),
             linetype = "solid", linewidth = 0.7) +
  scale_color_manual(values = c("95%" = "#FFA500", "90%" = "#FF0000", "85%" = "#9370DB")) +
  facet_wrap(~ Tag, ncol = 3) +
  labs(
    title = "NDVI Senescence Detection at Different Thresholds",
    subtitle = "Green dashed = peak NDVI; colored lines = senescence onset",
    x = "Day of Year",
    y = "NDVI",
    color = "Threshold"
  ) +
  theme_classic() +
  theme(legend.position = "top")

print(threshold_vis_ndvi)

# ----------------------------------
# Final Analysis with 90% Threshold
# ----------------------------------

senescence_ndvi_final <- senescence_ndvi_90pct

cat("\n\n=== FINAL SENESCENCE ONSET RESULTS (NDVI, 90% Threshold) ===\n\n")

# Summary statistics
senescence_ndvi_summary <- senescence_ndvi_final %>%
  group_by(Population, Treatment) %>%
  summarise(
    mean_peak_DOY = mean(DOY_at_peak, na.rm = TRUE),
    mean_onset_DOY = mean(DOY_senescence_onset, na.rm = TRUE),
    mean_duration = mean(DOY_senescence_onset - DOY_at_peak, na.rm = TRUE),
    sd_onset_DOY = sd(DOY_senescence_onset, na.rm = TRUE),
    se_onset_DOY = sd_onset_DOY / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

print(senescence_ndvi_summary)

# Statistical models
senescence_ndvi_lm <- lm(DOY_senescence_onset ~ Treatment * Population, 
                         data = senescence_ndvi_final)

duration_ndvi_lm <- lm(I(DOY_senescence_onset - DOY_at_peak) ~ Treatment * Population,
                       data = senescence_ndvi_final)

cat("\n\n=== SENESCENCE ONSET MODEL (NDVI) ===\n")
print(summary(senescence_ndvi_lm))
print(anova(senescence_ndvi_lm))

cat("\n\n=== DURATION MODEL (NDVI) ===\n")
print(summary(duration_ndvi_lm))
print(anova(duration_ndvi_lm))

# ----------------------------------
# Plots
# ----------------------------------

senescence_ndvi_plot <- ggplot(senescence_ndvi_summary, 
                               aes(x = Population, y = mean_onset_DOY, 
                                   color = Treatment, group = Treatment)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = mean_onset_DOY - se_onset_DOY, 
                    ymax = mean_onset_DOY + se_onset_DOY),
                width = 0.2, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(x = "Population", 
       y = "Senescence Onset (DOY)", 
       color = "Treatment",
       title = "NDVI-Based Senescence Onset by Population and Treatment",
       subtitle = "Defined as decline to 90% of peak NDVI") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(senescence_ndvi_plot)

duration_ndvi_plot <- ggplot(senescence_ndvi_summary, 
                             aes(x = Population, y = mean_duration, 
                                 color = Treatment, group = Treatment)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(x = "Population", 
       y = "Days from Peak to Senescence Onset", 
       color = "Treatment",
       title = "Duration of Green Period (NDVI-based)") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(duration_ndvi_plot)

# ----------------------------------
# Threshold comparison plot
# ----------------------------------

threshold_long_ndvi <- threshold_comparison_ndvi %>%
  pivot_longer(
    cols = c(mean_onset_95, mean_onset_90, mean_onset_85),
    names_to = "threshold",
    values_to = "mean_onset"
  ) %>%
  mutate(
    threshold = case_when(
      threshold == "mean_onset_95" ~ "95%",
      threshold == "mean_onset_90" ~ "90%",
      threshold == "mean_onset_85" ~ "85%"
    ),
    threshold = factor(threshold, levels = c("95%", "90%", "85%"))
  )

threshold_comp_plot_ndvi <- ggplot(threshold_long_ndvi, 
                                   aes(x = threshold, y = mean_onset, 
                                       color = Treatment, group = Treatment)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = TREATMENT_COLORS) +
  facet_wrap(~ Population) +
  labs(
    title = "NDVI Senescence Onset Date by Threshold Choice",
    x = "Threshold (% of Peak NDVI)",
    y = "Mean Senescence Onset (DOY)",
    color = "Treatment"
  ) +
  theme_classic() +
  theme(legend.position = "top")

print(threshold_comp_plot_ndvi)

# ----------------------------------
# Save Results
# ----------------------------------

write_csv(senescence_ndvi_final, "output/tables/senescence_onset_ndvi_by_block.csv")
write_csv(senescence_ndvi_summary, "output/tables/senescence_ndvi_summary.csv")
write_csv(threshold_comparison_ndvi, "output/tables/threshold_comparison_ndvi.csv")

ggsave("output/figures/senescence_threshold_comparison_ndvi.pdf", threshold_vis_ndvi, 
       width = 10, height = 10, dpi = 600)
ggsave("output/figures/senescence_threshold_comp_plot_ndvi.pdf", threshold_comp_plot_ndvi,
       width = 10, height = 6, dpi = 600)
ggsave("output/figures/senescence_onset_ndvi.pdf", senescence_ndvi_plot, 
       width = 8, height = 6, dpi = 600)
ggsave("output/figures/senescence_duration_ndvi.pdf", duration_ndvi_plot, 
       width = 8, height = 6, dpi = 600)

# ----------------------------------
# Calculate Senescence Onset (After Peak Greenness)
# ----------------------------------

# Method 1: Find peak greenness for each plant, then first decline below threshold
senescence_after_peak <- function(data, threshold_pct_of_peak = 0.95) {
  data %>%
    arrange(Tag, DOY) %>%
    group_by(Tag, Population, TreatmentName) %>%
    filter(!is.na(TillerPercentGreen)) %>%
    mutate(
      max_green = max(TillerPercentGreen, na.rm = TRUE),
      DOY_at_peak = DOY[which.max(TillerPercentGreen)][1],
      # Only look at dates after peak
      after_peak = DOY > DOY_at_peak,
      # Threshold as percentage of peak (e.g., 95% of max)
      threshold_value = max_green * threshold_pct_of_peak,
      below_threshold = TillerPercentGreen < threshold_value
    ) %>%
    # Only consider observations after peak
    filter(after_peak, below_threshold) %>%
    # Get first date after peak that drops below threshold
    slice(1) %>%
    dplyr::select(
      Tag, Population, TreatmentName, 
      DOY_senescence_onset = DOY,
      DOY_at_peak,
      max_green,
      PercentGreen_at_onset = TillerPercentGreen
    ) %>%
    ungroup()
}

# Try different thresholds (% of peak greenness)
senescence_onset_95pct <- senescence_after_peak(pheno, threshold_pct_of_peak = 0.95)
senescence_onset_90pct <- senescence_after_peak(pheno, threshold_pct_of_peak = 0.90)
senescence_onset_80pct <- senescence_after_peak(pheno, threshold_pct_of_peak = 0.80)

# Method 2: Maximum rate of decline AFTER peak
senescence_max_decline_after_peak <- function(data) {
  data %>%
    arrange(Tag, DOY) %>%
    group_by(Tag, Population, TreatmentName) %>%
    filter(!is.na(TillerPercentGreen)) %>%
    mutate(
      DOY_at_peak = DOY[which.max(TillerPercentGreen)][1],
      after_peak = DOY > DOY_at_peak,
      decline_rate = c(NA, diff(TillerPercentGreen) / diff(DOY))
    ) %>%
    # Only look after peak and where declining
    filter(after_peak, !is.na(decline_rate), decline_rate < 0) %>%
    slice_min(decline_rate, n = 1) %>%  # Most negative = steepest decline
    dplyr::select(
      Tag, Population, TreatmentName, 
      DOY_senescence_onset = DOY,
      DOY_at_peak,
      max_decline_rate = decline_rate
    ) %>%
    ungroup()
}

senescence_onset_decline <- senescence_max_decline_after_peak(pheno)

# Method 3: Use fitted curves - find inflection point AFTER peak
senescence_from_fitted_curve <- function(model, data) {
  data_with_pred <- data %>%
    filter(!is.na(TillerPercentGreen)) %>%
    mutate(predicted = predict(model, newdata = .))
  
  data_with_pred %>%
    arrange(Tag, DOY) %>%
    group_by(Tag, Population, TreatmentName) %>%
    mutate(
      DOY_at_peak = DOY[which.max(predicted)][1],
      after_peak = DOY > DOY_at_peak,
      first_deriv = c(NA, diff(predicted) / diff(DOY)),
      second_deriv = c(NA, NA, diff(first_deriv, lag = 1) / diff(DOY, lag = 1)[-1])
    ) %>%
    # Look for inflection point after peak where it's declining
    filter(after_peak, !is.na(second_deriv), first_deriv < 0) %>%
    slice_min(second_deriv, n = 1) %>%
    dplyr::select(
      Tag, Population, TreatmentName,
      DOY_senescence_onset = DOY,
      DOY_at_peak
    ) %>%
    ungroup()
}

senescence_inflection <- senescence_from_fitted_curve(percent_green_model, pheno)

# ----------------------------------
# Visualize to check the approach makes sense
# ----------------------------------

# Plot greenness over time with peak and senescence onset marked
example_plants <- pheno %>%
  filter(Tag %in% sample(unique(pheno$Tag), 6))  # Sample 6 random plants

example_senescence <- senescence_onset_90pct %>%
  filter(Tag %in% example_plants$Tag)

ggplot(example_plants, aes(x = DOY, y = TillerPercentGreen)) +
  geom_line(aes(group = Tag)) +
  geom_point() +
  geom_vline(data = example_senescence, 
             aes(xintercept = DOY_at_peak), 
             color = "green", linetype = "dashed") +
  geom_vline(data = example_senescence, 
             aes(xintercept = DOY_senescence_onset), 
             color = "red", linetype = "dashed") +
  facet_wrap(~ Tag, ncol = 3) +
  labs(
    title = "Example: Peak (green) vs Senescence Onset (red)",
    subtitle = "Senescence defined as 90% of peak greenness",
    x = "Day of Year",
    y = "Percent Green (%)"
  ) +
  theme_classic()

# ----------------------------------
# Summary Statistics
# ----------------------------------

senescence_summary <- senescence_onset_90pct %>%
  group_by(Population, TreatmentName) %>%
  summarise(
    mean_peak_DOY = mean(DOY_at_peak, na.rm = TRUE),
    mean_onset_DOY = mean(DOY_senescence_onset, na.rm = TRUE),
    mean_duration = mean(DOY_senescence_onset - DOY_at_peak, na.rm = TRUE),
    sd_onset_DOY = sd(DOY_senescence_onset, na.rm = TRUE),
    n = n(),
    se_onset_DOY = sd_onset_DOY / sqrt(n),
    .groups = "drop"
  )

# ----------------------------------
# Statistical Analysis
# ----------------------------------

# Now we can fit a proper mixed-effects model
# Use linear model since each plant has only one senescence date
senescence_lm <- lm(DOY_senescence_onset ~ TreatmentName * Population, 
                    data = senescence_onset_90pct)

# Also analyze the duration from peak to senescence
duration_lm <- lm(I(DOY_senescence_onset - DOY_at_peak) ~ TreatmentName * Population,
                  data = senescence_onset_90pct)

# ----------------------------------
# Visualization
# ----------------------------------

senescence_plot <- ggplot(senescence_summary, 
                          aes(x = Population, y = mean_onset_DOY, 
                              color = TreatmentName, group = TreatmentName)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = mean_onset_DOY - se_onset_DOY, 
                    ymax = mean_onset_DOY + se_onset_DOY),
                width = 0.2, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(x = "Population", 
       y = "Senescence Onset (DOY)", 
       color = "Treatment",
       title = "Date of Senescence Onset by Population and Treatment",
       subtitle = "Defined as 90% of peak greenness") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Duration plot
duration_plot <- ggplot(senescence_summary, 
                        aes(x = Population, y = mean_duration, 
                            color = TreatmentName, group = TreatmentName)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(x = "Population", 
       y = "Days from Peak to Senescence Onset", 
       color = "Treatment",
       title = "Duration of Green Period") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----------------------------------
# Save Results
# ----------------------------------

write_csv(senescence_onset_90pct, "output/tables/senescence_onset_90pct_of_peak.csv")
write_csv(senescence_onset_decline, "output/tables/senescence_onset_max_decline.csv")
write_csv(senescence_summary, "output/tables/senescence_summary.csv")

sink("output/tables/senescence_model_summary.txt")
cat("Senescence Onset Model (90% of peak):\n")
print(summary(senescence_lm))
print(anova(senescence_lm))

cat("\n\nDuration from Peak to Senescence Model:\n")
print(summary(duration_lm))
print(anova(duration_lm))
sink()

ggsave("output/figures/senescence_onset_plot.pdf", senescence_plot, 
       width = 8, height = 6, units = "in", dpi = 600)
ggsave("output/figures/senescence_duration_plot.pdf", duration_plot, 
       width = 8, height = 6, units = "in", dpi = 600)

# Polished table
senescence_summary %>%
  gt() %>%
  tab_header(
    title = "Senescence Onset by Population and Treatment",
    subtitle = "Senescence defined as decline to 90% of peak greenness"
  ) %>%
  fmt_number(columns = c(mean_peak_DOY, mean_onset_DOY, mean_duration, 
                         sd_onset_DOY, se_onset_DOY), decimals = 1) %>%
  cols_label(
    Population = "Population",
    TreatmentName = "Treatment",
    mean_peak_DOY = "Peak DOY",
    mean_onset_DOY = "Onset DOY",
    mean_duration = "Duration (days)",
    sd_onset_DOY = "SD",
    n = "N",
    se_onset_DOY = "SE"
  ) %>%
  gtsave("output/tables/senescence_summary.html")
