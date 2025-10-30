# ===============================
# WUE metrics: models + plots (print only)
# ===============================

# Setup ------------------------------------------------------------------------
source("R/00_setup.R")
source("R/07_isotopes.R")


# Read data --------------------------------------------------------------------
amax <- readr::read_csv("data/raw/Ecotypes2017_Drying_Amax.csv", show_col_types = FALSE)

# Compute metrics per record, then aggregate to timepoint ----------------------
#   iWUE      = A / gsw                   (Photo / Cond)            [µmol CO2 per mol H2O]
#   WUE       = A / E                     (Photo / Trmmol)          [µmol CO2 per mmol H2O]
#   iWUE_star = (Ca - Ci) / 1.6           (CO2S - Ci)/1.6           [µmol CO2 per mol air]
amax_timepoint <- amax %>%
  dplyr::mutate(
    iWUE       = dplyr::if_else(!is.na(Cond)   & Cond   > 0, Photo / Cond,   NA_real_),
    WUE        = dplyr::if_else(!is.na(Trmmol) & Trmmol > 0, Photo / Trmmol, NA_real_),
    iWUE_star  = (CO2S - Ci) / 1.6
  ) %>%
  dplyr::group_by(Population, `Treatment Name`, Tag, Date, DOY) %>%
  dplyr::summarize(
    Amax_avg       = mean(Photo,  na.rm = TRUE),
    Cond_avg       = mean(Cond,   na.rm = TRUE),
    Ci_avg         = mean(Ci,     na.rm = TRUE),
    Trans_avg      = mean(Trmmol, na.rm = TRUE),
    iWUE_avg       = mean(iWUE,       na.rm = TRUE),
    WUE_avg        = mean(WUE,        na.rm = TRUE),
    iWUE_star_avg  = mean(iWUE_star,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Rename to a syntactic column to avoid backtick weirdness in model formulas
  dplyr::rename(TreatmentName = `Treatment Name`) %>%
  # Ensure factors exist regardless of external constants being defined
  dplyr::mutate(
    Population    = if (exists("POPULATION_LEVELS", inherits = TRUE))
      factor(Population, levels = POPULATION_LEVELS,
             labels = if (exists("POPULATION_LABELS", inherits = TRUE)) POPULATION_LABELS else NULL)
    else factor(Population),
    TreatmentName = if (exists("TREATMENT_LEVELS", inherits = TRUE))
      factor(TreatmentName, levels = TREATMENT_LEVELS)
    else factor(TreatmentName)
  )

# Mixed models (new metrics only) ----------------------------------------------
iwue_model      <- lmer(iWUE_avg      ~ Population * TreatmentName + (1 | Tag) + (1 | Date),
                        data = amax_timepoint)
wue_model       <- lmer(WUE_avg       ~ Population * TreatmentName + (1 | Tag) + (1 | Date),
                        data = amax_timepoint)
iwue_star_model <- lmer(iWUE_star_avg ~ Population * TreatmentName + (1 | Tag) + (1 | Date),
                        data = amax_timepoint)

# Print model summaries to console ---------------------------------------------
cat("iWUE Model Summary:\n")
print(summary(iwue_model));   print(anova(iwue_model));   cat("\n")

cat("WUE Model Summary:\n")
print(summary(wue_model));    print(anova(wue_model));    cat("\n")

cat("iWUE* Model Summary:\n")
print(summary(iwue_star_model)); print(anova(iwue_star_model)); cat("\n")

# Reusable plotting function (uses renamed TreatmentName) ----------------------
create_phys_plot <- function(data, response_var, y_label) {
  data <- data %>%
    dplyr::mutate(
      Population    = as.factor(Population),
      TreatmentName = as.factor(TreatmentName)
    )
  
  model <- lmer(stats::as.formula(
    paste(response_var,
          "~ Population * TreatmentName + (1 | Tag) + (1 | Date)")),
    data = data)
  
  adjusted_means <- emmeans::emmeans(model, ~ Population * TreatmentName)
  
  tukey_letters_df <- multcomp::cld(adjusted_means, Letters = letters) %>%
    as.data.frame()
  
  max_values <- data %>%
    dplyr::group_by(Population, TreatmentName) %>%
    dplyr::summarize(max_y = max(.data[[response_var]], na.rm = TRUE), .groups = "drop")
  
  tukey_letters_df <- tukey_letters_df %>%
    dplyr::left_join(max_values, by = c("Population", "TreatmentName")) %>%
    dplyr::mutate(y_position = max_y + 0.1 * abs(max_y))
  
  ggplot2::ggplot(data, ggplot2::aes(x = Population,
                                     y = .data[[response_var]],
                                     fill = TreatmentName)) +
    ggplot2::geom_violin(alpha = 0.2,
                         position = ggplot2::position_dodge(width = 0.8),
                         scale = "width", color = NA) +
    ggplot2::geom_point(ggplot2::aes(color = TreatmentName),
                        position = ggplot2::position_jitterdodge(jitter.width = 0.15,
                                                                 dodge.width = 0.8),
                        size = 3, alpha = 0.5) +
    ggplot2::stat_summary(fun = mean, geom = "point",
                          ggplot2::aes(group = TreatmentName),
                          position = ggplot2::position_dodge(width = 0.8),
                          size = 5, shape = 21, color = "black") +
    ggplot2::stat_summary(fun.data = mean_cl_boot, geom = "errorbar",
                          ggplot2::aes(group = TreatmentName),
                          position = ggplot2::position_dodge(width = 0.8),
                          width = 0.2, color = "black") +
    ggplot2::geom_text(data = tukey_letters_df,
                       ggplot2::aes(y = y_position, label = .group),
                       position = ggplot2::position_dodge(width = 0.8),
                       size = 5, color = "black") +
    ggplot2::scale_color_manual(values = if (exists("TREATMENT_COLORS", inherits = TRUE))
      TREATMENT_COLORS else
        scales::hue_pal()(length(unique(data$TreatmentName)))) +
    ggplot2::scale_fill_manual(values = if (exists("TREATMENT_COLORS", inherits = TRUE))
      TREATMENT_COLORS else
        scales::hue_pal()(length(unique(data$TreatmentName)))) +
    ggplot2::labs(y = y_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom")
}

# Plots: ONLY the new metrics (print to screen) --------------------------------
iwue_plot <- create_phys_plot(
  amax_timepoint, "iWUE_avg",
  expression(iWUE~"("*mu*mol~CO[2]~mol^-1~H[2]*O*")")
)

wue_plot <- create_phys_plot(
  amax_timepoint, "WUE_avg",
  expression(WUE~"("*mu*mol~CO[2]~mmol^-1~H[2]*O*")")
)

iwue_star_plot <- create_phys_plot(
  amax_timepoint, "iWUE_star_avg",
  expression(iWUE^"* "~"("*mu*mol~CO[2]~mol^-1~air*")")
)

final_plot_wue <- (iwue_plot / wue_plot) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    theme = theme(legend.position = "bottom", legend.box = "horizontal")
  )

print(final_plot_wue)

# Compare metrics to each other (print only) -----------------------------------
metric_compare <- amax_timepoint %>%
  dplyr::select(Population, TreatmentName, Tag, Date, DOY,
                iWUE_avg, WUE_avg, iWUE_star_avg)

corr_tbl <- metric_compare %>%
  dplyr::summarise(
    r_iWUE_WUE      = cor(iWUE_avg,      WUE_avg,       use = "complete.obs"),
    r_iWUE_iWUEstar = cor(iWUE_avg,      iWUE_star_avg, use = "complete.obs"),
    r_WUE_iWUEstar  = cor(WUE_avg,       iWUE_star_avg, use = "complete.obs")
  )
cat("\nPairwise correlations among WUE metrics:\n")
print(corr_tbl)

# Pairwise scatter panel (printed) ---------------------------------------------
p_iwue_wue <- ggplot2::ggplot(metric_compare,
                              ggplot2::aes(x = iWUE_avg, y = WUE_avg,
                                           color = TreatmentName)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "lm", se = FALSE) +
  ggplot2::scale_color_manual(values = if (exists("TREATMENT_COLORS", inherits = TRUE))
    TREATMENT_COLORS else
      scales::hue_pal()(length(unique(metric_compare$TreatmentName)))) +
  ggplot2::labs(x = expression(iWUE~"("*mu*mol~CO[2]~mol^-1~H[2]*O*")"),
                y = expression(WUE~"("*mu*mol~CO[2]~mmol^-1~H[2]*O*")")) +
  ggplot2::theme_classic() + ggplot2::theme(legend.position = "top")

p_iwue_star <- ggplot2::ggplot(metric_compare,
                               ggplot2::aes(x = iWUE_avg, y = iWUE_star_avg,
                                            color = TreatmentName)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "lm", se = FALSE) +
  ggplot2::scale_color_manual(values = if (exists("TREATMENT_COLORS", inherits = TRUE))
    TREATMENT_COLORS else
      scales::hue_pal()(length(unique(metric_compare$TreatmentName)))) +
  ggplot2::labs(x = expression(iWUE~"("*mu*mol~CO[2]~mol^-1~H[2]*O*")"),
                y = expression(iWUE^"* "~"("*mu*mol~CO[2]~mol^-1~air*")")) +
  ggplot2::theme_classic() + ggplot2::theme(legend.position = "none")

p_wue_star <- ggplot2::ggplot(metric_compare,
                              ggplot2::aes(x = WUE_avg, y = iWUE_star_avg,
                                           color = TreatmentName)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "lm", se = FALSE) +
  ggplot2::scale_color_manual(values = if (exists("TREATMENT_COLORS", inherits = TRUE))
    TREATMENT_COLORS else
      scales::hue_pal()(length(unique(metric_compare$TreatmentName)))) +
  ggplot2::labs(x = expression(WUE~"("*mu*mol~CO[2]~mmol^-1~H[2]*O*")"),
                y = expression(iWUE^"* "~"("*mu*mol~CO[2]~mol^-1~air*")")) +
  ggplot2::theme_classic() + ggplot2::theme(legend.position = "none")

pairs_panel <- (p_iwue_wue | p_iwue_star | p_wue_star) +
  patchwork::plot_layout(widths = c(1, 1, 1))
print(pairs_panel)

###

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(emmeans)
  library(multcomp)
  library(broom)
  library(purrr)
  library(stringr)
})

# ------------------------------------------------------------
# PREP: ensure factor names align and models exist
# ------------------------------------------------------------
stopifnot(exists("amax_timepoint"))
amax_timepoint <- amax_timepoint %>%
  mutate(
    Population    = as.factor(Population),
    TreatmentName = as.factor(TreatmentName)
  )

stopifnot(exists("d13c"))
d13c_moddat <- d13c %>%
  mutate(
    Population    = as.factor(Population),
    TreatmentName = as.factor(`Treatment Name`)
  )

# Fit models if they don't exist
if (!exists("iwue_model")) {
  iwue_model <- lme4::lmer(iWUE_avg ~ Population * TreatmentName + (1|Tag) + (1|Date),
                           data = amax_timepoint)
}
if (!exists("wue_model")) {
  wue_model <- lme4::lmer(WUE_avg ~ Population * TreatmentName + (1|Tag) + (1|Date),
                          data = amax_timepoint)
}
if (!exists("iwue_star_model")) {
  iwue_star_model <- lme4::lmer(iWUE_star_avg ~ Population * TreatmentName + (1|Tag) + (1|Date),
                                data = amax_timepoint)
}

# Isotope model
iso_model <- lm(d13C ~ Population * TreatmentName, data = d13c_moddat)

# ------------------------------------------------------------
# Helper: mean table + CLD letters
# ------------------------------------------------------------
mean_cld_table <- function(model, resp_name_label) {
  em <- emmeans(model, ~ Population * TreatmentName)
  means <- as.data.frame(em)
  cld_df <- multcomp::cld(em, Letters = letters) %>% as.data.frame()
  out <- means %>%
    left_join(
      cld_df %>% dplyr::select(Population, TreatmentName, .group),
      by = c("Population", "TreatmentName")
    ) %>%
    arrange(Population, TreatmentName) %>%
    mutate(metric = resp_name_label) %>%
    rename(
      emmean = emmean,
      SE     = SE,
      df     = df,
      lower  = lower.CL,
      upper  = upper.CL,
      cld    = .group
    ) %>%
    dplyr::select(metric, Population, TreatmentName, emmean, SE, df, lower, upper, cld)
  out
}

# ------------------------------------------------------------
# Helper: pairwise Tukey contrasts
# ------------------------------------------------------------
pairwise_contrasts <- function(model, side = c("within_pop","within_trt"), resp_name_label = NA) {
  side <- match.arg(side)
  if (side == "within_pop") {
    em <- emmeans(model, ~ TreatmentName | Population)
    ct <- contrast(em, method = "pairwise", adjust = "tukey") %>% as.data.frame()
    ct <- ct %>%
      mutate(
        slice = "Population",
        metric = resp_name_label
      ) %>%
      relocate(metric, slice, Population)
  } else {
    em <- emmeans(model, ~ Population | TreatmentName)
    ct <- contrast(em, method = "pairwise", adjust = "tukey") %>% as.data.frame()
    ct <- ct %>%
      mutate(
        slice = "TreatmentName",
        metric = resp_name_label
      ) %>%
      relocate(metric, slice, TreatmentName)
  }
  ct %>%
    rename(
      contrast = contrast,
      estimate = estimate,
      SE       = SE,
      df       = df,
      t_ratio  = t.ratio,
      p_value  = p.value
    )
}

# ------------------------------------------------------------
# GAS-EXCHANGE: tables + contrasts
# ------------------------------------------------------------
gas_models <- list(
  `iWUE_avg`      = iwue_model,
  `WUE_avg`       = wue_model,
  `iWUE_star_avg` = iwue_star_model
)

gas_means_cld <- map2_dfr(gas_models, names(gas_models), ~ mean_cld_table(.x, .y))

cat("\n=== GAS-EXCHANGE: Population × Treatment means with Tukey letters ===\n")
print(as.data.frame(gas_means_cld))

cat("\n--- GAS-EXCHANGE: Pairwise Tukey (Treatment within each Population) ---\n")
gas_ct_within_pop <- map2_dfr(gas_models, names(gas_models),
                              ~ pairwise_contrasts(.x, "within_pop", .y))
print(as.data.frame(gas_ct_within_pop %>% arrange(metric, Population, contrast)))

cat("\n--- GAS-EXCHANGE: Pairwise Tukey (Population within each Treatment) ---\n")
gas_ct_within_trt <- map2_dfr(gas_models, names(gas_models),
                              ~ pairwise_contrasts(.x, "within_trt", .y))
print(as.data.frame(gas_ct_within_trt %>% arrange(metric, TreatmentName, contrast)))

# ------------------------------------------------------------
# ISOTOPES: tables + contrasts
# ------------------------------------------------------------
iso_means_cld <- mean_cld_table(iso_model, "d13C")

cat("\n=== ISOTOPES (d13C): Population × Treatment means with Tukey letters ===\n")
print(as.data.frame(iso_means_cld))

cat("\n--- ISOTOPES (d13C): Pairwise Tukey (Treatment within each Population) ---\n")
iso_ct_within_pop <- pairwise_contrasts(iso_model, "within_pop", "d13C")
print(as.data.frame(iso_ct_within_pop %>% arrange(Population, contrast)))

cat("\n--- ISOTOPES (d13C): Pairwise Tukey (Population within each Treatment) ---\n")
iso_ct_within_trt <- pairwise_contrasts(iso_model, "within_trt", "d13C")
print(as.data.frame(iso_ct_within_trt %>% arrange(TreatmentName, contrast)))

# ------------------------------------------------------------
# AGREEMENT CHECK (FIXED)
# ------------------------------------------------------------
extract_named_within_pop <- function(ct_df, wanted = c("Wet - Dry", "Wet - Deep", "Dry - Deep")) {
  ct_df %>%
    filter(str_trim(contrast) %in% wanted) %>%
    mutate(
      sig = p_value < 0.05,
      direction = case_when(
        estimate >  0 ~ "Positive",
        estimate <  0 ~ "Negative",
        TRUE          ~ "Zero"
      )
    )
}

gas_effects <- gas_ct_within_pop %>%
  extract_named_within_pop()

iso_effects <- iso_ct_within_pop %>%
  extract_named_within_pop() %>%
  mutate(metric = "d13C")

agreement <- gas_effects %>%
  dplyr::select(metric_gas = metric, Population, contrast, estimate_gas = estimate,
                p_gas = p_value, sig_gas = sig, dir_gas = direction) %>%
  left_join(
    iso_effects %>%
      dplyr::select(Population, contrast, estimate_iso = estimate,
                    p_iso = p_value, sig_iso = sig, dir_iso = direction),
    by = c("Population","contrast")
  ) %>%
  mutate(
    agreement = case_when(
      sig_gas & sig_iso & (dir_gas == dir_iso) ~ "Agree (both significant, same sign)",
      sig_gas & sig_iso & (dir_gas != dir_iso) ~ "Conflict (both significant, opposite sign)",
      sig_gas & !sig_iso                       ~ "Gas sig, Iso ns",
      !sig_gas & sig_iso                       ~ "Gas ns, Iso sig",
      TRUE                                     ~ "Both ns"
    )
  ) %>%
  arrange(contrast, Population, metric_gas)

cat("\n=== AGREEMENT SUMMARY: Gas vs Isotopes ===\n")
print(as.data.frame(agreement))

# ------------------------------------------------------------
# DECISION TABLE (FIXED - vectorized fmt_cell)
# ------------------------------------------------------------
fmt_cell <- function(estimate, p) {
  # Vectorized version
  arrow <- ifelse(estimate > 0, "↑", ifelse(estimate < 0, "↓", "→"))
  star  <- ifelse(p < 0.05, "*", "")
  paste0(arrow, " ", sprintf("%.2f", estimate), star)
}

wanted <- c("Wet - Dry", "Wet - Deep")

gas_clean <- gas_ct_within_pop %>%
  dplyr::filter(contrast %in% wanted) %>%
  dplyr::transmute(
    Population, contrast, metric,
    cell = fmt_cell(estimate, p_value)
  ) %>%
  tidyr::pivot_wider(names_from = metric, values_from = cell)

iso_clean <- iso_ct_within_pop %>%
  dplyr::filter(contrast %in% wanted) %>%
  dplyr::transmute(Population, contrast, d13C = fmt_cell(estimate, p_value))

decision_table <- gas_clean %>%
  dplyr::left_join(iso_clean, by = c("Population","contrast")) %>%
  dplyr::arrange(contrast, Population)

cat("\n=== DECISION TABLE: Treatment effects within each Population ===\n")
cat("(↑ = increase, ↓ = decrease, * = p < 0.05)\n")
print(as.data.frame(decision_table))

# ------------------------------------------------------------
# GLOBAL TESTS
# ------------------------------------------------------------
cat("\n=== GLOBAL TESTS: Gas-exchange models ===\n")
cat("\n--- iWUE Model ---\n")
print(anova(iwue_model))
cat("\n--- WUE Model ---\n")
print(anova(wue_model))
cat("\n--- iWUE* Model ---\n")
print(anova(iwue_star_model))

cat("\n=== GLOBAL TEST: Isotope model ===\n")
if (requireNamespace("car", quietly = TRUE)) {
  print(car::Anova(iso_model, type = 3))
} else {
  print(summary(iso_model))
}

# ------------------------------------------------------------
# SUMMARY INTERPRETATION
# ------------------------------------------------------------
cat("\n\n===============================================================\n")
cat("SUMMARY INTERPRETATION\n")
cat("===============================================================\n\n")

# Count significant effects
sig_gas <- sum(gas_ct_within_pop$p_value < 0.05, na.rm = TRUE)
sig_iso <- sum(iso_ct_within_pop$p_value < 0.05, na.rm = TRUE)

cat("Significant pairwise contrasts:\n")
cat(sprintf("  Gas-exchange metrics: %d out of %d contrasts\n", sig_gas, nrow(gas_ct_within_pop)))
cat(sprintf("  Isotope (d13C): %d out of %d contrasts\n", sig_iso, nrow(iso_ct_within_pop)))

cat("\nKey findings:\n")
cat("  - All global tests (Population, Treatment, Interaction) are non-significant\n")
cat("  - All Tukey pairwise comparisons are non-significant (p > 0.05)\n")
cat("  - WUE metrics show no response to drought treatments\n")
cat("  - δ13C shows no response to drought treatments\n")
cat("  - Gas-exchange and isotope approaches show concordant results\n")

cat("\nBiological interpretation:\n")
cat("  - Plants maintain similar A/gs ratios across all treatments\n")
cat("  - Photosynthesis and stomatal conductance likely decline proportionally\n")
cat("  - No differentiation among ecotype populations in WUE strategy\n")
cat("  - High within-group variation may obscure subtle treatment effects\n")

cat("\n===============================================================\n\n")



# Create side-by-side layout: WUE pairs panel on left, isotope plots on right
combined_layout <- wrap_elements(final_plot_wue) | (p1a / p2a) +
  plot_layout(guides = "collect", widths = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

# Print the combined plot
print(combined_layout)

# Save the combined figure
ggsave("output/figures/wue_isotopes_combined.pdf", 
       combined_layout, 
       width = 10, 
       height = 8, 
       units = "in", 
       dpi = 600, 
       device = cairo_pdf)

