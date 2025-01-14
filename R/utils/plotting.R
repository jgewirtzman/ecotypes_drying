
# utils/plotting.R
# plot_treatment_response <- function(data, response_var, y_label) {
#   ggplot(data, aes(x = Population, y = .data[[response_var]], fill = `Treatment Name`)) +
#     geom_violin(alpha = 0.2, position = position_dodge(width = 0.8)) +
#     geom_point(aes(color = `Treatment Name`),
#                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
#                size = 3, alpha = 0.5) +
#     scale_color_manual(values = TREATMENT_COLORS) +
#     scale_fill_manual(values = TREATMENT_COLORS) +
#     labs(y = y_label) +
#     theme_paper()
# }
# 
# plot_time_series <- function(data, response_var, y_label) {
#   ggplot(data, aes(x = DOY, y = .data[[response_var]], color = `Treatment Name`)) +
#     geom_line(aes(group = Tag), alpha = 0.3, size = 0.5) +
#     stat_summary(fun = mean, geom = "line", size = 1.5) +
#     stat_summary(fun.data = mean_cl_boot, geom = "ribbon", 
#                  aes(fill = `Treatment Name`), color = NA, alpha = 0.3) +
#     scale_color_manual(values = TREATMENT_COLORS) +
#     scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
#     labs(x = "Day of Year", y = y_label) +
#     theme_paper() +
#     facet_wrap(~ Population, ncol = 1)
# }


# utils/plotting.R

# Function to plot treatment responses
plot_treatment_response <- function(data, response_var, y_label) {
  ggplot(data, aes(x = Population, y = .data[[response_var]], fill = TreatmentName)) +
    geom_violin(alpha = 0.3, position = position_dodge(width = 0.8)) +
    geom_jitter(aes(color = TreatmentName), position = position_jitterdodge(), size = 2, alpha = 0.6) +
    scale_fill_manual(values = TREATMENT_COLORS) +
    scale_color_manual(values = TREATMENT_COLORS) +
    labs(x = "Population", y = y_label, fill = "Treatment") +
    theme_paper()
}

# Function to plot phenology metrics over time
plot_pheno_metric <- function(data, metric, y_label) {
  ggplot(data, aes(x = DOY, y = .data[[metric]], color = TreatmentName)) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.3) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    scale_color_manual(values = TREATMENT_COLORS) +
    labs(x = "Day of Year", y = y_label, color = "Treatment") +
    theme_paper() +
    facet_wrap(~ Population, ncol = 3)
}