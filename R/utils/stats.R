# utils/stats.R

# Function to clean factor levels for Population and Treatment columns
clean_factor_levels <- function(data) {
  data %>%
    mutate(
      Population = factor(Population, levels = POPULATION_LEVELS, labels = POPULATION_LABELS),
      TreatmentName = factor(`Treatment Name`, levels = TREATMENT_LEVELS)
    )
}

# Function to fit a mixed-effects model
fit_mixed_model <- function(data, response_var, include_date = TRUE) {
  formula <- if (include_date) {
    as.formula(paste(response_var, "~ Population * TreatmentName + poly(DOY, 2) + (1 | Tag) + (1 | Date)"))
  } else {
    as.formula(paste(response_var, "~ Population * TreatmentName + poly(DOY, 2) + (1 | Tag)"))
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
