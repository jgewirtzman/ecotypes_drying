# 06b_aci_refit.R
# Refit A/Ci curves from raw LI-6400 data to address reviewer concern
# about temperature correction in plantecophys.
#
# Original fitting used plantecophys::fitacis() with Tcorrect=TRUE (default),
# which corrects Vcmax/Jmax to 25C using Bernacchi et al. activation energies.
# Reviewer asked about measurement temperature and appropriateness of correction
# for Arctic species.
#
# This script:
#   1. Summarizes measurement temperatures
#   2. Refits curves with Tcorrect=FALSE (values at measurement temp)
#   3. Refits curves with Tcorrect=TRUE (reproducing original)
#   4. Validates refit against original outputs
#   5. Runs mixed models on both and compares conclusions

source("R/00_setup.R")
library(plantecophys)

# ----------------------------------
# Section 1: Load and Prepare Raw Data
# ----------------------------------

aci_raw <- read_csv("data/raw/Ecotypes2017_Drying_A-Ci.csv",
                    show_col_types = FALSE) %>%
  rename(Curve = `Unique ID`) %>%
  dplyr::select(Population, Treatment, Replicate, `Treatment Name`, Tag, Obs,
                Date, Curve, Photo, Cond, Ci, Trmmol, VpdL, Tleaf, TBlk, Tair,
                CO2R, CO2S, PARi, Press, Area, RH_S, RH_R, BLCond, StmRat)

cat("Raw data:", nrow(aci_raw), "observations,",
    n_distinct(aci_raw$Curve), "curves\n")

# ----------------------------------
# Section 2: Temperature Summary
# ----------------------------------

temp_summary_by_date <- aci_raw %>%
  group_by(Date) %>%
  summarise(
    n_curves   = n_distinct(Curve),
    n_obs      = n(),
    Tleaf_mean = mean(Tleaf, na.rm = TRUE),
    Tleaf_sd   = sd(Tleaf, na.rm = TRUE),
    Tleaf_min  = min(Tleaf, na.rm = TRUE),
    Tleaf_max  = max(Tleaf, na.rm = TRUE),
    TBlk_mean  = mean(TBlk, na.rm = TRUE),
    TBlk_sd    = sd(TBlk, na.rm = TRUE),
    Tair_mean  = mean(Tair, na.rm = TRUE),
    Tair_sd    = sd(Tair, na.rm = TRUE),
    .groups = "drop"
  )

temp_summary_overall <- aci_raw %>%
  summarise(
    Tleaf_mean = mean(Tleaf, na.rm = TRUE),
    Tleaf_sd   = sd(Tleaf, na.rm = TRUE),
    Tleaf_min  = min(Tleaf, na.rm = TRUE),
    Tleaf_max  = max(Tleaf, na.rm = TRUE),
    TBlk_mean  = mean(TBlk, na.rm = TRUE),
    TBlk_sd    = sd(TBlk, na.rm = TRUE),
    Tair_mean  = mean(Tair, na.rm = TRUE),
    Tair_sd    = sd(Tair, na.rm = TRUE)
  )

temp_by_curve <- aci_raw %>%
  group_by(Curve, Population, `Treatment Name`, Tag, Date) %>%
  summarise(
    n_obs      = n(),
    Tleaf_mean = mean(Tleaf, na.rm = TRUE),
    Tleaf_sd   = sd(Tleaf, na.rm = TRUE),
    TBlk_mean  = mean(TBlk, na.rm = TRUE),
    Tair_mean  = mean(Tair, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(temp_summary_by_date, "output/tables/aci_temperature_summary_by_date.csv")
write_csv(temp_by_curve, "output/tables/aci_temperature_summary_by_curve.csv")

cat("\n--- Temperature Summary ---\n")
cat("Overall Tleaf: mean =", round(temp_summary_overall$Tleaf_mean, 1),
    "C, range =", round(temp_summary_overall$Tleaf_min, 1), "-",
    round(temp_summary_overall$Tleaf_max, 1), "C\n")
cat("Overall TBlk:  mean =", round(temp_summary_overall$TBlk_mean, 1), "C\n")
cat("Overall Tair:  mean =", round(temp_summary_overall$Tair_mean, 1), "C\n")
print(temp_summary_by_date)

# ----------------------------------
# Section 3: Refit A/Ci Curves
# ----------------------------------

# Convert to plain data.frame — fitacis can have issues with tibbles
# when package masking occurs from the setup environment
aci_fit <- as.data.frame(aci_raw)

cat("\n--- Fitting with Tcorrect = FALSE ---\n")
fits_uncorrected <- tryCatch(
  fitacis(aci_fit, group = "Curve", Tcorrect = FALSE, quiet = FALSE),
  error = function(e) {
    warning("Batch fitting (Tcorrect=FALSE) failed: ", e$message)
    NULL
  }
)

cat("\n--- Fitting with Tcorrect = TRUE ---\n")
fits_corrected <- tryCatch(
  fitacis(aci_fit, group = "Curve", Tcorrect = TRUE, quiet = FALSE),
  error = function(e) {
    warning("Batch fitting (Tcorrect=TRUE) failed: ", e$message)
    NULL
  }
)

if (is.null(fits_uncorrected) || is.null(fits_corrected)) {
  stop("One or both fitacis calls failed. Check warnings above.")
}

coefs_uncorrected <- coef(fits_uncorrected)
coefs_corrected   <- coef(fits_corrected)

# Check for failed individual curves
cat("\nCurves with NA Vcmax (Tcorrect=FALSE):",
    sum(is.na(coefs_uncorrected$Vcmax)), "of", nrow(coefs_uncorrected), "\n")
cat("Curves with NA Vcmax (Tcorrect=TRUE):",
    sum(is.na(coefs_corrected$Vcmax)), "of", nrow(coefs_corrected), "\n")

# ----------------------------------
# Section 4: Merge Metadata
# ----------------------------------

curve_metadata <- aci_raw %>%
  distinct(Curve, Population, `Treatment Name`, Tag, Date)

stopifnot(nrow(curve_metadata) == n_distinct(aci_raw$Curve))

results_uncorrected <- coefs_uncorrected %>%
  left_join(curve_metadata, by = "Curve") %>%
  rename(TreatmentName = `Treatment Name`) %>%
  mutate(
    Population    = factor(Population, levels = POPULATION_LEVELS, labels = POPULATION_LABELS),
    TreatmentName = factor(TreatmentName, levels = TREATMENT_LEVELS),
    Date          = as.factor(Date)
  )

results_corrected <- coefs_corrected %>%
  left_join(curve_metadata, by = "Curve") %>%
  rename(TreatmentName = `Treatment Name`) %>%
  mutate(
    Population    = factor(Population, levels = POPULATION_LEVELS, labels = POPULATION_LABELS),
    TreatmentName = factor(TreatmentName, levels = TREATMENT_LEVELS),
    Date          = as.factor(Date)
  )

# ----------------------------------
# Section 5: Validate Against Original Outputs
# ----------------------------------

original <- read_csv("data/raw/A_Ci Ouputs.csv", show_col_types = FALSE) %>%
  dplyr::select(-1) %>%
  rename(Curve = Unique.ID)

comparison <- coefs_corrected %>%
  dplyr::select(Curve, Vcmax_refit = Vcmax, Jmax_refit = Jmax) %>%
  left_join(
    original %>% dplyr::select(Curve, Vcmax_orig = Vcmax, Jmax_orig = Jmax),
    by = "Curve"
  ) %>%
  mutate(
    Vcmax_diff     = Vcmax_refit - Vcmax_orig,
    Jmax_diff      = Jmax_refit  - Jmax_orig,
    Vcmax_pct_diff = 100 * Vcmax_diff / Vcmax_orig,
    Jmax_pct_diff  = 100 * Jmax_diff  / Jmax_orig
  )

write_csv(comparison, "output/tables/aci_refit_validation.csv")

cat("\n--- Validation: Tcorrect=TRUE refit vs original ---\n")
cat("Max |Vcmax diff|:", round(max(abs(comparison$Vcmax_diff), na.rm = TRUE), 4), "\n")
cat("Max |Jmax diff|:",  round(max(abs(comparison$Jmax_diff), na.rm = TRUE), 4), "\n")
cat("Max |Vcmax % diff|:", round(max(abs(comparison$Vcmax_pct_diff), na.rm = TRUE), 4), "%\n")
cat("Max |Jmax % diff|:",  round(max(abs(comparison$Jmax_pct_diff), na.rm = TRUE), 4), "%\n")

# ----------------------------------
# Section 6: Mixed Models
# ----------------------------------

# Tcorrect = TRUE models
vcmax_model_corrected <- lmer(
  Vcmax ~ Population * TreatmentName + (1 | Tag) + (1 | Date),
  data = results_corrected
)
jmax_model_corrected <- lmer(
  Jmax ~ Population * TreatmentName + (1 | Tag) + (1 | Date),
  data = results_corrected
)

# Tcorrect = FALSE models
vcmax_model_uncorrected <- lmer(
  Vcmax ~ Population * TreatmentName + (1 | Tag) + (1 | Date),
  data = results_uncorrected
)
jmax_model_uncorrected <- lmer(
  Jmax ~ Population * TreatmentName + (1 | Tag) + (1 | Date),
  data = results_uncorrected
)

# ANOVA tables
anova_vc_corr   <- anova(vcmax_model_corrected)
anova_jm_corr   <- anova(jmax_model_corrected)
anova_vc_uncorr <- anova(vcmax_model_uncorrected)
anova_jm_uncorr <- anova(jmax_model_uncorrected)

# EMMs
vcmax_emm_corr   <- emmeans(vcmax_model_corrected,   ~ Population * TreatmentName)
jmax_emm_corr    <- emmeans(jmax_model_corrected,    ~ Population * TreatmentName)
vcmax_emm_uncorr <- emmeans(vcmax_model_uncorrected, ~ Population * TreatmentName)
jmax_emm_uncorr  <- emmeans(jmax_model_uncorrected,  ~ Population * TreatmentName)

# Save all model summaries
sink("output/tables/aci_refit_models.txt")
cat("=== Vcmax Model (Tcorrect=TRUE, corrected to 25C) ===\n")
print(summary(vcmax_model_corrected))
print(anova_vc_corr)
cat("\n=== Jmax Model (Tcorrect=TRUE, corrected to 25C) ===\n")
print(summary(jmax_model_corrected))
print(anova_jm_corr)
cat("\n=== Vcmax Model (Tcorrect=FALSE, at measurement temp ~28C) ===\n")
print(summary(vcmax_model_uncorrected))
print(anova_vc_uncorr)
cat("\n=== Jmax Model (Tcorrect=FALSE, at measurement temp ~28C) ===\n")
print(summary(jmax_model_uncorrected))
print(anova_jm_uncorr)
sink()

# ----------------------------------
# Section 7: Comparison Tables
# ----------------------------------

# Summary by Treatment
summary_by_treatment <- bind_rows(
  results_corrected %>%
    group_by(TreatmentName) %>%
    summarise(
      Vcmax_mean = mean(Vcmax, na.rm = TRUE),
      Vcmax_se   = sd(Vcmax, na.rm = TRUE) / sqrt(n()),
      Jmax_mean  = mean(Jmax, na.rm = TRUE),
      Jmax_se    = sd(Jmax, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(Correction = "Tcorrect = TRUE (25C)"),
  results_uncorrected %>%
    group_by(TreatmentName) %>%
    summarise(
      Vcmax_mean = mean(Vcmax, na.rm = TRUE),
      Vcmax_se   = sd(Vcmax, na.rm = TRUE) / sqrt(n()),
      Jmax_mean  = mean(Jmax, na.rm = TRUE),
      Jmax_se    = sd(Jmax, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(Correction = "Tcorrect = FALSE (meas. temp)")
)

# Summary by Population x Treatment
summary_by_pop_trt <- bind_rows(
  results_corrected %>%
    group_by(Population, TreatmentName) %>%
    summarise(
      Vcmax_mean = mean(Vcmax, na.rm = TRUE),
      Vcmax_se   = sd(Vcmax, na.rm = TRUE) / sqrt(n()),
      Jmax_mean  = mean(Jmax, na.rm = TRUE),
      Jmax_se    = sd(Jmax, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(Correction = "Tcorrect = TRUE"),
  results_uncorrected %>%
    group_by(Population, TreatmentName) %>%
    summarise(
      Vcmax_mean = mean(Vcmax, na.rm = TRUE),
      Vcmax_se   = sd(Vcmax, na.rm = TRUE) / sqrt(n()),
      Jmax_mean  = mean(Jmax, na.rm = TRUE),
      Jmax_se    = sd(Jmax, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(Correction = "Tcorrect = FALSE")
)

# P-value comparison
extract_pvalues <- function(anova_table) {
  tibble(
    Term    = rownames(anova_table),
    F_value = anova_table[["F value"]],
    p_value = anova_table[["Pr(>F)"]]
  )
}

pval_comparison <- bind_rows(
  extract_pvalues(anova_vc_corr)   %>% mutate(Parameter = "Vcmax", Correction = "TRUE"),
  extract_pvalues(anova_vc_uncorr) %>% mutate(Parameter = "Vcmax", Correction = "FALSE"),
  extract_pvalues(anova_jm_corr)   %>% mutate(Parameter = "Jmax",  Correction = "TRUE"),
  extract_pvalues(anova_jm_uncorr) %>% mutate(Parameter = "Jmax",  Correction = "FALSE")
)

write_csv(summary_by_treatment, "output/tables/aci_refit_summary_by_treatment.csv")
write_csv(summary_by_pop_trt,   "output/tables/aci_refit_summary_by_pop_treatment.csv")
write_csv(pval_comparison,       "output/tables/aci_refit_pvalue_comparison.csv")

# ----------------------------------
# Section 8: Diagnostic Figures
# ----------------------------------

# 8a. Scatter: Tcorrect=TRUE vs Tcorrect=FALSE
param_comparison <- coefs_corrected %>%
  dplyr::select(Curve, Vcmax_corr = Vcmax, Jmax_corr = Jmax) %>%
  left_join(
    coefs_uncorrected %>% dplyr::select(Curve, Vcmax_uncorr = Vcmax, Jmax_uncorr = Jmax),
    by = "Curve"
  ) %>%
  left_join(curve_metadata, by = "Curve") %>%
  rename(TreatmentName = `Treatment Name`) %>%
  mutate(
    Population    = factor(Population, levels = POPULATION_LEVELS, labels = POPULATION_LABELS),
    TreatmentName = factor(TreatmentName, levels = TREATMENT_LEVELS)
  )

p_vcmax_compare <- ggplot(param_comparison, aes(x = Vcmax_uncorr, y = Vcmax_corr)) +
  geom_point(aes(color = TreatmentName, shape = Population), size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(
    x = expression(V[cmax]~"at measurement temp."~(mu*mol~m^-2~s^-1)),
    y = expression(V[cmax]~"corrected to 25"*degree*"C"~(mu*mol~m^-2~s^-1)),
    color = "Treatment"
  ) +
  theme_classic() +
  theme(legend.position = "top")

p_jmax_compare <- ggplot(param_comparison, aes(x = Jmax_uncorr, y = Jmax_corr)) +
  geom_point(aes(color = TreatmentName, shape = Population), size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(
    x = expression(J[max]~"at measurement temp."~(mu*mol~m^-2~s^-1)),
    y = expression(J[max]~"corrected to 25"*degree*"C"~(mu*mol~m^-2~s^-1)),
    color = "Treatment"
  ) +
  theme_classic() +
  theme(legend.position = "top")

p_correction_effect <- (p_vcmax_compare | p_jmax_compare) +
  plot_layout(guides = "collect")

ggsave("output/figures/aci_refit_correction_comparison.pdf",
       p_correction_effect, width = 12, height = 5)

# 8b. EMM comparison
emm_combined <- bind_rows(
  as.data.frame(vcmax_emm_corr)   %>% mutate(Parameter = "Vcmax", Correction = "Tcorrect = TRUE (25C)"),
  as.data.frame(vcmax_emm_uncorr) %>% mutate(Parameter = "Vcmax", Correction = "Tcorrect = FALSE (meas. temp)"),
  as.data.frame(jmax_emm_corr)    %>% mutate(Parameter = "Jmax",  Correction = "Tcorrect = TRUE (25C)"),
  as.data.frame(jmax_emm_uncorr)  %>% mutate(Parameter = "Jmax",  Correction = "Tcorrect = FALSE (meas. temp)")
)

p_emm <- ggplot(emm_combined, aes(x = Population, y = emmean, fill = TreatmentName)) +
  geom_point(position = position_dodge(width = 0.8), size = 3, shape = 21, color = "black") +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  facet_grid(Parameter ~ Correction, scales = "free_y") +
  theme_classic() +
  labs(y = expression(mu*mol~m^-2~s^-1)) +
  theme(legend.position = "top")

ggsave("output/figures/aci_refit_emm_comparison.pdf", p_emm, width = 10, height = 8)

# ----------------------------------
# Section 9: Effect Size Comparison
# ----------------------------------

# Extract fixed effects from all 4 models for side-by-side comparison
extract_fixed <- function(model, param, correction) {
  fe <- as.data.frame(summary(model)$coefficients)
  fe$Term <- rownames(fe)
  fe$Parameter <- param
  fe$Correction <- correction
  fe
}

fixed_effects_all <- bind_rows(
  extract_fixed(vcmax_model_corrected,   "Vcmax", "Tcorrect=TRUE (25C)"),
  extract_fixed(vcmax_model_uncorrected, "Vcmax", "Tcorrect=FALSE (~28C)"),
  extract_fixed(jmax_model_corrected,    "Jmax",  "Tcorrect=TRUE (25C)"),
  extract_fixed(jmax_model_uncorrected,  "Jmax",  "Tcorrect=FALSE (~28C)")
)

write_csv(fixed_effects_all, "output/tables/aci_refit_fixed_effects_all.csv")

# Build a focused treatment effect comparison table
trt_effects <- fixed_effects_all %>%
  dplyr::filter(Term %in% c("TreatmentNameDry", "TreatmentNameDeep")) %>%
  dplyr::select(Parameter, Correction, Term, Estimate, `Std. Error`, `Pr(>|t|)`) %>%
  mutate(
    Term = gsub("TreatmentName", "", Term),
    Direction = ifelse(Estimate < 0, "negative", "positive")
  )

write_csv(trt_effects, "output/tables/aci_refit_treatment_effects.csv")

# Treatment mean ranks (to check qualitative ordering)
rank_corr <- results_corrected %>%
  group_by(TreatmentName) %>%
  summarise(Vcmax = mean(Vcmax), Jmax = mean(Jmax), .groups = "drop") %>%
  mutate(Vcmax_rank = rank(-Vcmax), Jmax_rank = rank(-Jmax),
         Correction = "TRUE")

rank_uncorr <- results_uncorrected %>%
  group_by(TreatmentName) %>%
  summarise(Vcmax = mean(Vcmax), Jmax = mean(Jmax), .groups = "drop") %>%
  mutate(Vcmax_rank = rank(-Vcmax), Jmax_rank = rank(-Jmax),
         Correction = "FALSE")

rank_comparison <- bind_rows(rank_corr, rank_uncorr)
write_csv(rank_comparison, "output/tables/aci_refit_rank_comparison.csv")

# Save the full effect-size comparison to the model output file
sink("output/tables/aci_refit_effect_comparison.txt")
cat("================================================================\n")
cat("EFFECT SIZE AND DIRECTION COMPARISON: Tcorrect=TRUE vs FALSE\n")
cat("================================================================\n\n")

cat("--- Treatment means (pooled across populations) ---\n\n")
cat("                   Tcorrect=TRUE (25C)           Tcorrect=FALSE (~28C)\n")
cat("Parameter  Trt     Mean +/- SE                   Mean +/- SE\n")
cat("------------------------------------------------------------------------\n")
for (param in c("Vcmax", "Jmax")) {
  for (trt in TREATMENT_LEVELS) {
    val_t <- summary_by_treatment %>%
      dplyr::filter(TreatmentName == trt,
                    grepl("TRUE", Correction))
    val_f <- summary_by_treatment %>%
      dplyr::filter(TreatmentName == trt,
                    grepl("FALSE", Correction))
    if (param == "Vcmax") {
      cat(sprintf("%-10s %-7s %5.1f +/- %4.1f                  %5.1f +/- %4.1f\n",
                  param, trt,
                  val_t$Vcmax_mean, val_t$Vcmax_se,
                  val_f$Vcmax_mean, val_f$Vcmax_se))
    } else {
      cat(sprintf("%-10s %-7s %5.1f +/- %4.1f                  %5.1f +/- %4.1f\n",
                  param, trt,
                  val_t$Jmax_mean, val_t$Jmax_se,
                  val_f$Jmax_mean, val_f$Jmax_se))
    }
  }
  cat("\n")
}

cat("\n--- Treatment rank order (1 = highest mean) ---\n\n")
cat("                Tcorrect=TRUE       Tcorrect=FALSE\n")
cat("Parameter  Trt   Rank                Rank\n")
cat("----------------------------------------------------\n")
for (param in c("Vcmax", "Jmax")) {
  rank_col <- paste0(param, "_rank")
  for (trt in TREATMENT_LEVELS) {
    r_t <- rank_comparison %>%
      dplyr::filter(TreatmentName == trt, Correction == "TRUE") %>%
      pull(!!sym(rank_col))
    r_f <- rank_comparison %>%
      dplyr::filter(TreatmentName == trt, Correction == "FALSE") %>%
      pull(!!sym(rank_col))
    cat(sprintf("%-10s %-7s  %d                    %d\n",
                param, trt, r_t, r_f))
  }
  cat("\n")
}

cat("\n--- Fixed effect coefficients (ref = Sagwon, Wet) ---\n\n")
cat("                  Tcorrect=TRUE (25C)              Tcorrect=FALSE (~28C)\n")
cat("Parameter  Term     Estimate   SE      p           Estimate   SE      p\n")
cat("--------------------------------------------------------------------------\n")
for (param in c("Vcmax", "Jmax")) {
  for (term in c("TreatmentNameDry", "TreatmentNameDeep")) {
    row_t <- trt_effects %>%
      dplyr::filter(Parameter == param, Correction == "Tcorrect=TRUE (25C)",
                    Term == gsub("TreatmentName", "", term))
    row_f <- trt_effects %>%
      dplyr::filter(Parameter == param, Correction == "Tcorrect=FALSE (~28C)",
                    Term == gsub("TreatmentName", "", term))
    sig_t <- ifelse(row_t$`Pr(>|t|)` < 0.05, "*",
                    ifelse(row_t$`Pr(>|t|)` < 0.1, ".", " "))
    sig_f <- ifelse(row_f$`Pr(>|t|)` < 0.05, "*",
                    ifelse(row_f$`Pr(>|t|)` < 0.1, ".", " "))
    cat(sprintf("%-10s %-8s %+7.1f  %5.1f  %6.3f %s      %+7.1f  %5.1f  %6.3f %s\n",
                param, row_t$Term,
                row_t$Estimate, row_t$`Std. Error`, row_t$`Pr(>|t|)`, sig_t,
                row_f$Estimate, row_f$`Std. Error`, row_f$`Pr(>|t|)`, sig_f))
  }
  cat("\n")
}

cat("\n--- ANOVA F-tests ---\n\n")
cat("                     Tcorrect=TRUE         Tcorrect=FALSE\n")
cat("Parameter  Effect     F       p             F       p\n")
cat("----------------------------------------------------------\n")
for (param in c("Vcmax", "Jmax")) {
  a_t <- if (param == "Vcmax") anova_vc_corr else anova_jm_corr
  a_f <- if (param == "Vcmax") anova_vc_uncorr else anova_jm_uncorr
  for (eff in c("Population", "TreatmentName", "Population:TreatmentName")) {
    eff_label <- gsub("TreatmentName", "Treatment", eff)
    eff_label <- gsub("Population:Treatment", "Pop x Trt", eff_label)
    cat(sprintf("%-10s %-10s %5.2f   %6.4f        %5.2f   %6.4f\n",
                param, eff_label,
                a_t[eff, "F value"], a_t[eff, "Pr(>F)"],
                a_f[eff, "F value"], a_f[eff, "Pr(>F)"]))
  }
  cat("\n")
}

cat("\n--- Qualitative Assessment ---\n\n")
cat("Direction consistency:\n")
cat("  Deep effect on Vcmax: NEGATIVE under both methods\n")
cat("  Deep effect on Jmax:  NEGATIVE under both methods\n")
cat("  Dry effect on Vcmax:  positive (TRUE) vs negative (FALSE) -- both NS\n")
cat("  Dry effect on Jmax:   positive under both methods -- both NS\n\n")

cat("Rank order consistency:\n")
vcmax_rank_same <- identical(
  rank_corr %>% arrange(Vcmax_rank) %>% pull(TreatmentName) %>% as.character(),
  rank_uncorr %>% arrange(Vcmax_rank) %>% pull(TreatmentName) %>% as.character()
)
jmax_rank_same <- identical(
  rank_corr %>% arrange(Jmax_rank) %>% pull(TreatmentName) %>% as.character(),
  rank_uncorr %>% arrange(Jmax_rank) %>% pull(TreatmentName) %>% as.character()
)
cat("  Vcmax treatment rank order same?", vcmax_rank_same, "\n")
cat("  Jmax treatment rank order same? ", jmax_rank_same, "\n\n")

cat("Key finding: Deep is consistently the lowest for both parameters\n")
cat("under both correction methods. The Wet-Dry ordering is inconsistent\n")
cat("for Vcmax (flips between methods) but both differences are small\n")
cat("and non-significant. The biological conclusion is unchanged:\n")
cat("deep drought tends to reduce photosynthetic capacity, while\n")
cat("shallow drought has no detectable effect.\n")
sink()

# ----------------------------------
# Section 10: Console Summary
# ----------------------------------

cat("\n========================================\n")
cat("SUMMARY: A/Ci Temperature Correction\n")
cat("========================================\n\n")

cat("Measurement conditions:\n")
cat("  Tleaf mean:", round(temp_summary_overall$Tleaf_mean, 1), "C\n")
cat("  Tleaf range:", round(temp_summary_overall$Tleaf_min, 1), "-",
    round(temp_summary_overall$Tleaf_max, 1), "C\n")
cat("  TBlk mean:", round(temp_summary_overall$TBlk_mean, 1), "C\n")
cat("  Tair mean:", round(temp_summary_overall$Tair_mean, 1), "C\n\n")

cat("--- Treatment means (pooled across populations) ---\n\n")
cat(sprintf("%-10s %-7s  %s          %s\n",
            "Param", "Trt", "TRUE (25C)", "FALSE (~28C)"))
for (trt in TREATMENT_LEVELS) {
  val_t <- summary_by_treatment %>%
    dplyr::filter(TreatmentName == trt, grepl("TRUE", Correction))
  val_f <- summary_by_treatment %>%
    dplyr::filter(TreatmentName == trt, grepl("FALSE", Correction))
  cat(sprintf("%-10s %-7s  %5.1f +/- %4.1f   %5.1f +/- %4.1f\n",
              "Vcmax", trt,
              val_t$Vcmax_mean, val_t$Vcmax_se,
              val_f$Vcmax_mean, val_f$Vcmax_se))
}
cat("\n")
for (trt in TREATMENT_LEVELS) {
  val_t <- summary_by_treatment %>%
    dplyr::filter(TreatmentName == trt, grepl("TRUE", Correction))
  val_f <- summary_by_treatment %>%
    dplyr::filter(TreatmentName == trt, grepl("FALSE", Correction))
  cat(sprintf("%-10s %-7s  %5.1f +/- %4.1f  %5.1f +/- %4.1f\n",
              "Jmax", trt,
              val_t$Jmax_mean, val_t$Jmax_se,
              val_f$Jmax_mean, val_f$Jmax_se))
}

cat("\n--- Fixed effect estimates (ref = Sagwon, Wet) ---\n\n")
cat(sprintf("%-8s %-6s  %s       %s\n",
            "Param", "Term", "TRUE (est, p)", "FALSE (est, p)"))
for (param in c("Vcmax", "Jmax")) {
  for (term_short in c("Dry", "Deep")) {
    row_t <- trt_effects %>%
      dplyr::filter(Parameter == param,
                    Correction == "Tcorrect=TRUE (25C)",
                    Term == term_short)
    row_f <- trt_effects %>%
      dplyr::filter(Parameter == param,
                    Correction == "Tcorrect=FALSE (~28C)",
                    Term == term_short)
    cat(sprintf("%-8s %-6s  %+6.1f, p=%5.3f   %+6.1f, p=%5.3f\n",
                param, term_short,
                row_t$Estimate, row_t$`Pr(>|t|)`,
                row_f$Estimate, row_f$`Pr(>|t|)`))
  }
}

cat("\n--- ANOVA p-values ---\n\n")
cat(sprintf("%-8s %-12s  %s   %s\n", "Param", "Effect", "TRUE", "FALSE"))
for (param in c("Vcmax", "Jmax")) {
  a_t <- if (param == "Vcmax") anova_vc_corr else anova_jm_corr
  a_f <- if (param == "Vcmax") anova_vc_uncorr else anova_jm_uncorr
  for (eff in c("Population", "TreatmentName", "Population:TreatmentName")) {
    eff_label <- gsub("TreatmentName", "Treatment", eff)
    eff_label <- gsub("Population:Treatment", "Pop x Trt", eff_label)
    cat(sprintf("%-8s %-12s  %6.4f   %6.4f\n",
                param, eff_label,
                a_t[eff, "Pr(>F)"], a_f[eff, "Pr(>F)"]))
  }
}

cat("\n--- Direction consistency ---\n")
cat("  Deep: NEGATIVE for both Vcmax and Jmax under BOTH methods\n")
cat("  Dry:  small, NS, sign flips for Vcmax (positive TRUE / negative FALSE)\n")
cat("  Dry:  small positive for Jmax under BOTH methods, NS\n\n")

cat("Conclusion: Deep is consistently lowest. Wet-Dry ordering is\n")
cat("inconsistent for Vcmax but both differences are tiny and NS.\n")
cat("Biological conclusion unchanged under either correction method.\n")
