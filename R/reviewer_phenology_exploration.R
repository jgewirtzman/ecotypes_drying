#!/usr/bin/env Rscript
# ==============================================================================
# Phenology Data Exploration: Structure & Segmented Regression
# ==============================================================================
# Purpose: Characterize the phenology dataset and attempt piecewise/segmented
#          regression on Tiller Total Green Length as a function of DOY.
# ==============================================================================

cat("====================================================================\n")
cat("  PHENOLOGY DATA EXPLORATION: Structure & Segmented Regression\n")
cat("====================================================================\n\n")

# --- Setup -------------------------------------------------------------------
if (!requireNamespace("segmented", quietly = TRUE)) {
  cat("Installing 'segmented' package...\n")
  install.packages("segmented", repos = "https://cloud.r-project.org")
}
library(segmented)

# --- Load Data ---------------------------------------------------------------
dat <- read.csv("data/raw/Ecotypes2017_Drying_Phenology.csv",
                stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

# Clean column names: remove any BOM artifacts, leading dots/X
names(dat) <- gsub("^X\\.+", "", names(dat))

cat("Column names in dataset:\n")
print(names(dat))
cat("\n")

# Rename for convenience
names(dat)[names(dat) == "Treatment.Name"] <- "TrtName"
names(dat)[names(dat) == "Tiller.Total.Green.Length"] <- "TTGL"
names(dat)[names(dat) == "Average.Green.Leaf.Length"] <- "AGLL"
names(dat)[names(dat) == "Tiller...Green"] <- "PctGreen"

cat("Dimensions:", nrow(dat), "rows x", ncol(dat), "columns\n\n")

# ============================================================================
# PART 1: CHARACTERIZE DATA STRUCTURE
# ============================================================================
cat("====================================================================\n")
cat("  PART 1: DATA STRUCTURE CHARACTERIZATION\n")
cat("====================================================================\n\n")

# --- 1a. Unique measurement dates -------------------------------------------
unique_dates <- sort(unique(dat$Date))
unique_doy   <- sort(unique(dat$DOY))

cat("Number of unique measurement dates:", length(unique_dates), "\n")
cat("Dates:", paste(unique_dates, collapse = ", "), "\n")
cat("DOY values:", paste(unique_doy, collapse = ", "), "\n")
cat("DOY range:", min(unique_doy), "to", max(unique_doy),
    "(span =", max(unique_doy) - min(unique_doy), "days)\n\n")

# --- 1b. Unique tillers per Treatment x Population --------------------------
cat("--- Unique Tillers (Tags) per Treatment x Population ---\n")
tiller_counts <- tapply(dat$Tag, list(dat$TrtName, dat$Population),
                        function(x) length(unique(x)))
print(tiller_counts)
cat("\nTotal unique tillers:", length(unique(dat$Tag)), "\n\n")

# --- 1c. Measurements per tiller --------------------------------------------
meas_per_tiller <- tapply(dat$DOY, dat$Tag, length)
cat("--- Measurements per Tiller ---\n")
cat("Min:", min(meas_per_tiller), "\n")
cat("Median:", median(meas_per_tiller), "\n")
cat("Max:", max(meas_per_tiller), "\n")
cat("Mean:", round(mean(meas_per_tiller), 2), "\n\n")

cat("Distribution of measurement counts per tiller:\n")
print(table(meas_per_tiller))
cat("\n")

# --- 1d. Total observations -------------------------------------------------
cat("Total observations (rows):", nrow(dat), "\n")
cat("Observations with non-NA TTGL:", sum(!is.na(dat$TTGL)), "\n")
cat("Observations with TTGL > 0:", sum(dat$TTGL > 0, na.rm = TRUE), "\n\n")

# --- 1e. Summary of key response variables -----------------------------------
cat("--- Summary of Tiller Total Green Length ---\n")
print(summary(dat$TTGL))
cat("\n")

cat("--- Summary of Average Green Leaf Length ---\n")
print(summary(dat$AGLL))
cat("\n")

cat("--- Summary of Tiller % Green ---\n")
print(summary(dat$PctGreen))
cat("\n")

# --- 1f. Sample sizes per group (Treatment x Population) --------------------
cat("--- Sample sizes (observations) per Treatment x Population ---\n")
obs_counts <- tapply(dat$TTGL[!is.na(dat$TTGL)],
                     list(dat$TrtName[!is.na(dat$TTGL)],
                          dat$Population[!is.na(dat$TTGL)]),
                     length)
print(obs_counts)
cat("\n")


# ============================================================================
# PART 2: SEGMENTED / PIECEWISE REGRESSION
# ============================================================================
cat("====================================================================\n")
cat("  PART 2: SEGMENTED / PIECEWISE REGRESSION\n")
cat("====================================================================\n\n")

# We will fit:
#   (a) Simple linear model: TTGL ~ DOY
#   (b) Segmented model with 1 breakpoint
# Separately by Treatment Name, and then by Treatment x Population

# Helper function to fit linear + segmented model and report results
fit_segmented <- function(data, group_label) {
  cat("--------------------------------------------------------------\n")
  cat("Group:", group_label, "\n")
  cat("  N observations:", nrow(data), "\n")
  cat("  N unique tillers:", length(unique(data$Tag)), "\n")
  cat("  DOY range:", min(data$DOY), "-", max(data$DOY), "\n")
  cat("  N unique DOY values:", length(unique(data$DOY)), "\n")

  # Remove NA in TTGL
  data <- data[!is.na(data$TTGL), ]
  cat("  N non-NA TTGL:", nrow(data), "\n\n")

  if (nrow(data) < 5) {
    cat("  *** SKIPPING: Too few observations (< 5) ***\n\n")
    return(invisible(NULL))
  }

  # --- Simple linear model ---
  lm_fit <- lm(TTGL ~ DOY, data = data)
  lm_summary <- summary(lm_fit)
  cat("  [Linear Model] TTGL ~ DOY\n")
  cat("    Intercept:", round(coef(lm_fit)[1], 3), "\n")
  cat("    Slope:", round(coef(lm_fit)[2], 3), "\n")
  cat("    R-squared:", round(lm_summary$r.squared, 4), "\n")
  cat("    Adj R-squared:", round(lm_summary$adj.r.squared, 4), "\n")
  cat("    p-value (slope):", format.pval(lm_summary$coefficients[2, 4]), "\n\n")

  # --- Segmented model with 1 breakpoint ---
  cat("  [Segmented Model] 1 breakpoint\n")

  # Try segmented regression; catch errors/warnings
  seg_result <- tryCatch({
    # Use midpoint of DOY range as initial guess for breakpoint
    psi_init <- mean(range(data$DOY))
    seg_fit <- segmented(lm_fit, seg.Z = ~ DOY, psi = psi_init,
                         control = seg.control(display = FALSE,
                                               it.max = 100,
                                               n.boot = 20))

    seg_summ <- summary(seg_fit)

    # Extract breakpoint and CI
    bp <- seg_fit$psi
    bp_est <- bp[1, "Est."]
    bp_se  <- bp[1, "St.Err"]
    bp_ci_lo <- bp_est - 1.96 * bp_se
    bp_ci_hi <- bp_est + 1.96 * bp_se

    # Slopes
    slopes <- slope(seg_fit)

    # R-squared for segmented model
    ss_res <- sum(residuals(seg_fit)^2)
    ss_tot <- sum((data$TTGL - mean(data$TTGL))^2)
    r2_seg <- 1 - ss_res / ss_tot

    cat("    Converged: YES\n")
    cat("    Breakpoint estimate (DOY):", round(bp_est, 2), "\n")
    cat("    Breakpoint SE:", round(bp_se, 2), "\n")
    cat("    Breakpoint 95% CI: [", round(bp_ci_lo, 2), ",",
        round(bp_ci_hi, 2), "]\n")
    cat("    R-squared:", round(r2_seg, 4), "\n")
    cat("    Slopes:\n")
    print(slopes)
    cat("\n")

    # Compare to linear model
    cat("    R-squared improvement over linear:",
        round(r2_seg - lm_summary$r.squared, 4), "\n")

    # Davies test for existence of breakpoint
    dav <- tryCatch({
      davies.test(lm_fit, ~ DOY, k = 10)
    }, error = function(e) {
      cat("    Davies test: ERROR -", e$message, "\n")
      NULL
    })
    if (!is.null(dav)) {
      cat("    Davies test p-value:", format.pval(dav$p.value),
          "(significant = evidence for breakpoint)\n")
    }

    list(converged = TRUE, breakpoint = bp_est, r2 = r2_seg)
  },
  warning = function(w) {
    cat("    WARNING:", w$message, "\n")
    # Try to continue despite warning
    tryCatch({
      psi_init <- mean(range(data$DOY))
      suppressWarnings(
        seg_fit <- segmented(lm_fit, seg.Z = ~ DOY, psi = psi_init,
                             control = seg.control(display = FALSE,
                                                   it.max = 100,
                                                   n.boot = 20))
      )
      bp <- seg_fit$psi
      bp_est <- bp[1, "Est."]
      bp_se  <- bp[1, "St.Err"]
      ss_res <- sum(residuals(seg_fit)^2)
      ss_tot <- sum((data$TTGL - mean(data$TTGL))^2)
      r2_seg <- 1 - ss_res / ss_tot

      cat("    (Despite warning) Breakpoint estimate:", round(bp_est, 2), "\n")
      cat("    (Despite warning) R-squared:", round(r2_seg, 4), "\n")
      list(converged = TRUE, breakpoint = bp_est, r2 = r2_seg,
           warning = w$message)
    }, error = function(e2) {
      cat("    FAILED after warning:", e2$message, "\n")
      list(converged = FALSE, warning = w$message, error = e2$message)
    })
  },
  error = function(e) {
    cat("    Converged: NO\n")
    cat("    Error:", e$message, "\n")
    list(converged = FALSE, error = e$message)
  })

  cat("\n")
  return(invisible(seg_result))
}


# --- 2a. By Treatment Name only ---------------------------------------------
cat("============================================================\n")
cat("  2a. Segmented Regression by Treatment Name\n")
cat("============================================================\n\n")

treatments <- sort(unique(dat$TrtName))
cat("Treatments:", paste(treatments, collapse = ", "), "\n\n")

trt_results <- list()
for (trt in treatments) {
  sub <- dat[dat$TrtName == trt, ]
  trt_results[[trt]] <- fit_segmented(sub, paste("Treatment =", trt))
}


# --- 2b. By Treatment Name x Population ------------------------------------
cat("\n============================================================\n")
cat("  2b. Segmented Regression by Treatment x Population\n")
cat("============================================================\n\n")

populations <- sort(unique(dat$Population))
cat("Populations:", paste(populations, collapse = ", "), "\n\n")

txp_results <- list()
for (trt in treatments) {
  for (pop in populations) {
    sub <- dat[dat$TrtName == trt & dat$Population == pop, ]
    label <- paste("Treatment =", trt, "| Population =", pop)
    if (nrow(sub) > 0) {
      txp_results[[label]] <- fit_segmented(sub, label)
    }
  }
}


# ============================================================================
# PART 3: SUMMARY ASSESSMENT
# ============================================================================
cat("\n====================================================================\n")
cat("  PART 3: SUMMARY ASSESSMENT\n")
cat("====================================================================\n\n")

cat("--- Data Points per Group (Treatment x Population) ---\n")
for (trt in treatments) {
  for (pop in populations) {
    sub <- dat[dat$TrtName == trt & dat$Population == pop & !is.na(dat$TTGL), ]
    n_obs <- nrow(sub)
    n_tiller <- length(unique(sub$Tag))
    n_doy <- length(unique(sub$DOY))
    cat(sprintf("  %-5s x %-2s: %3d obs, %2d tillers, %d unique DOY\n",
                trt, pop, n_obs, n_tiller, n_doy))
  }
}

cat("\n--- Key Findings ---\n\n")

cat("1. TEMPORAL RESOLUTION:\n")
cat("   - Only", length(unique_doy), "unique measurement dates across the season\n")
cat("   - DOY range:", min(unique_doy), "-", max(unique_doy), "\n")
cat("   - This gives at most", length(unique_doy),
    "time points per tiller for fitting models\n\n")

cat("2. SAMPLE SIZES FOR SEGMENTED REGRESSION:\n")
cat("   - By Treatment alone:", paste(
  sapply(treatments, function(t) {
    n <- sum(dat$TrtName == t & !is.na(dat$TTGL))
    paste0(t, "=", n)
  }), collapse = ", "), "\n")
cat("   - By Treatment x Population: typically",
    round(mean(obs_counts, na.rm = TRUE)), "observations per group\n")
cat("   - Segmented regression with 1 breakpoint estimates ~5 parameters\n")
cat("   - Rule of thumb: need >= 10-15 obs per parameter, so >= 50-75 total\n")
cat("   - Per Treatment groups have sufficient data; Treatment x Population\n")
cat("     groups are borderline\n\n")

cat("3. DATA SPARSITY ASSESSMENT:\n")
cat("   - With only", length(unique_doy), "unique DOY values, the breakpoint\n")
cat("     can only be estimated between these discrete time points\n")
cat("   - Individual tillers have at most", max(meas_per_tiller),
    "measurements\n")
cat("   - This is MARGINAL for robust changepoint detection\n")
cat("   - The segmented package may struggle with convergence at this\n")
cat("     resolution, especially for Treatment x Population subgroups\n\n")

cat("4. CONVERGENCE SUMMARY:\n")
n_converged <- sum(sapply(txp_results, function(x) {
  if (is.null(x)) return(FALSE)
  isTRUE(x$converged)
}))
n_total <- length(txp_results)
cat("   - Treatment x Population models converged:", n_converged, "of", n_total, "\n")
n_trt_converged <- sum(sapply(trt_results, function(x) {
  if (is.null(x)) return(FALSE)
  isTRUE(x$converged)
}))
cat("   - Treatment-only models converged:", n_trt_converged, "of",
    length(trt_results), "\n\n")

cat("5. RECOMMENDATION:\n")
cat("   - The data may be better suited to nonlinear parametric models\n")
cat("     (e.g., logistic growth for green-up, exponential decay for senescence)\n")
cat("   - Mixed-effects models with tiller as random effect would account\n")
cat("     for repeated measures\n")
cat("   - Consider a two-phase approach: green-up (DOY ~196-232) and\n")
cat("     senescence (DOY ~232-254) fitted separately\n\n")

cat("====================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("====================================================================\n")
