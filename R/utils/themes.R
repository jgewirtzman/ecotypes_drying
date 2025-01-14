
# utils/stats.R
clean_factor_levels <- function(data) {
  data %>%
    mutate(
      Population = factor(Population, 
                          levels = POPULATION_LEVELS, 
                          labels = POPULATION_LABELS),
      `Treatment Name` = factor(`Treatment Name`, 
                                levels = TREATMENT_LEVELS)
    )
}

fit_mixed_model <- function(data, response_var, include_date = TRUE) {
  formula <- if(include_date) {
    as.formula(paste(response_var, "~ Population * `Treatment Name` + (1 | Tag) + (1 | Date)"))
  } else {
    as.formula(paste(response_var, "~ Population * `Treatment Name` + (1 | Tag)"))
  }
  
  lmer(formula, data = data)
}

calculate_auc <- function(data, value_col, time_col = "DOY") {
  data %>%
    group_by(Population, `Treatment Name`, Tag) %>%
    summarize(
      auc = pracma::trapz(.data[[time_col]], .data[[value_col]]),
      .groups = "drop"
    )
}

# utils/plotting.R
plot_treatment_response <- function(data, response_var, y_label) {
  ggplot(data, aes(x = Population, y = .data[[response_var]], fill = `Treatment Name`)) +
    geom_violin(alpha = 0.2, position = position_dodge(width = 0.8)) +
    geom_point(aes(color = `Treatment Name`),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               size = 3, alpha = 0.5) +
    scale_color_manual(values = TREATMENT_COLORS) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    labs(y = y_label) +
    theme_paper()
}

plot_time_series <- function(data, response_var, y_label) {
  ggplot(data, aes(x = DOY, y = .data[[response_var]], color = `Treatment Name`)) +
    geom_line(aes(group = Tag), alpha = 0.3, size = 0.5) +
    stat_summary(fun = mean, geom = "line", size = 1.5) +
    stat_summary(fun.data = mean_cl_boot, geom = "ribbon", 
                 aes(fill = `Treatment Name`), color = NA, alpha = 0.3) +
    scale_color_manual(values = TREATMENT_COLORS) +
    scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
    labs(x = "Day of Year", y = y_label) +
    theme_paper() +
    facet_wrap(~ Population, ncol = 1)
}