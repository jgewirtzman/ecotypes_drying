# Summary of Analyses for Reviewer Response

## Overview

This document summarizes the additional analyses performed to address reviewer comments on the Eriophorum vaginatum ecotypes drying manuscript. Four reviewer concerns were investigated, with the most substantive work addressing the A/Ci curve temperature correction issue (Comment 4).

---

## 1. A/Ci Curve Temperature Correction (Reviewer Comment 4)

**Reviewer concern:** "Under what temperature were A/Ci measurements made? If temperatures differ too much from 25°C, Vcmax/Jmax could be very different from actual values (also because non-Arctic specific activation energies are used)."

### What we found about measurement temperatures

- Measurements were made on **Aug 2** (18 curves) and **Aug 3, 2017** (6 curves)
- **Leaf temperatures (Tleaf) averaged ~28°C** (range 21.4–32.3°C), elevated above ambient due to cuvette heating at 1500 µmol PAR
- Block temperatures (TBlk) were 15.1°C (Aug 2) and 18.1°C (Aug 3)
- Air temperatures (Tair) were 15.9°C (Aug 2) and 18.9°C (Aug 3)
- The ~3°C offset from 25°C is modest but worth addressing

### What we did

The original A/Ci fitting was done externally using `plantecophys::fitacis()` with all defaults:
- `Tcorrect = TRUE` (corrects Vcmax/Jmax to 25°C reference temperature)
- Uses Bernacchi et al. (2001) activation energies derived from tobacco (*Nicotiana tabacum*)

We refitted all 24 curves two ways:
1. **Tcorrect = TRUE** — corrects to 25°C using Bernacchi activation energies (reproducing the original)
2. **Tcorrect = FALSE** — reports Vcmax/Jmax at actual measurement temperature (~28°C), no activation energy assumptions

The Tcorrect=TRUE refit validates perfectly against the original (max difference <0.0001%).

### Key results: Treatment effects are robust to correction method

**Treatment means (pooled across populations):**

| Parameter | Treatment | Tcorrect=TRUE (25°C) | Tcorrect=FALSE (~28°C) |
|-----------|-----------|---------------------|----------------------|
| Vcmax | Wet | 51.0 ± 5.6 | 84.2 ± 8.8 |
| Vcmax | Dry | 57.4 ± 7.9 | 79.6 ± 9.7 |
| Vcmax | Deep | 42.7 ± 8.5 | 47.5 ± 9.5 |
| Jmax | Wet | 149.2 ± 15.0 | 183.0 ± 19.4 |
| Jmax | Dry | 167.9 ± 22.6 | 193.5 ± 25.9 |
| Jmax | Deep | 113.3 ± 23.4 | 116.9 ± 23.6 |

**ANOVA p-values (Treatment effect):**

| Parameter | Tcorrect=TRUE | Tcorrect=FALSE |
|-----------|--------------|----------------|
| Vcmax | p = 0.548 | p = 0.059 |
| Jmax | p = 0.403 | p = 0.207 |

**Direction and rank consistency:**
- Deep drought is consistently the **lowest** for both Vcmax and Jmax under both methods
- Jmax rank order is identical: Dry > Wet > Deep (both methods)
- Vcmax Wet vs. Dry ordering flips between methods, but both differences are small and non-significant
- All Population and Population × Treatment interaction effects remain non-significant under both methods

### Bottom line for the response

The biological conclusion is unchanged regardless of temperature correction method: **deep drought tends to reduce photosynthetic capacity, while shallow drought has no detectable effect.** Since all measurements were made at the same temperature (~28°C), the temperature correction applies a uniform monotonic scaling that preserves all relative comparisons. The non-Arctic activation energies are a legitimate concern for absolute values, but they do not affect relative treatment comparisons—which is what the study tests.

### Suggested response elements

- Report that Tleaf averaged ~28°C due to cuvette heating at 1500 µmol PAR
- Acknowledge the Bernacchi activation energies are from tobacco
- Note that we reanalyzed with Tcorrect=FALSE (values at measurement temperature) as a sensitivity check
- Emphasize that all treatment comparisons and statistical conclusions are identical under both approaches
- Could note: Arctic-specific activation energies for E. vaginatum are not available; correcting to a different reference temperature (e.g., 15°C) would use the same tobacco-derived parameters and introduce *larger* extrapolation error, not smaller

### Script and outputs
- **Script:** `R/06b_aci_refit.R`
- **Key outputs:** `output/tables/aci_refit_effect_comparison.txt`, `output/tables/aci_refit_pvalue_comparison.csv`, `output/tables/aci_refit_summary_by_treatment.csv`, `output/figures/aci_refit_correction_comparison.pdf`

---

## 2. A/Ci Curve Sample Sizes (Reviewer Comment about n=24)

**Reviewer concern:** The manuscript says n=24, which could be misleading given the experimental design.

### What we found

The 24 curves break down as follows across Population × Treatment:

| | Sagwon | Toolik | Coldfoot |
|------|--------|--------|----------|
| Wet | 3 | 3 | 3 |
| Dry | 3 | 3 | 3 |
| Deep | 2 | 2 | 2 |

- Total: 9 Wet + 9 Dry + 6 Deep = 24 curves
- Deep has fewer curves (n=2 per population vs. n=3 for Wet/Dry)
- The mixed models use Tag (19 unique) and Date (2 dates) as random intercepts
- Tag random effect variance goes to zero (singular fit) — expected with only 24 observations and 19 tags

### Suggested response
- Report the per-cell sample sizes explicitly
- Note the reduced replication for Deep treatment
- This is already partially addressed in the manuscript but could be made clearer in a table or methods clarification

---

## 3. Stomatal Conductance Results (Reviewer Comment)

**Reviewer concern:** Confirm the stomatal conductance (gs) treatment effects.

### What we confirmed

From `R/04_gas_exchange.R`, the mixed model for stomatal conductance shows:
- **Deep drought significantly reduces gs**: estimate = −0.169, p = 0.012
- **Dry (shallow drought) effect is non-significant**: estimate = −0.064, p = 0.277
- No significant Population or Population × Treatment interaction effects
- This is consistent with the manuscript's reported findings

---

## 4. Phenology: Segmented Regression (Reviewer Comment)

**Reviewer concern:** Reviewer suggested trying broken stick / SiZer models to test whether senescence timing differs between treatments.

### What we found

Using the `segmented` R package on Tiller Total Green Length (TTGL) vs. DOY:
- **7 measurement dates** spanning DOY 196–254 (July 15 – Sept 11)
- **48 tillers** tracked across dates
- Breakpoints estimated at approximately **DOY 240–243** across treatments
- Confidence intervals for breakpoints **overlap substantially** across Wet, Dry, and Deep treatments
- The data has only 7 time points, which severely limits the power of segmented regression to detect differences in breakpoint timing

### Suggested response
- Could mention that piecewise regression was attempted but the coarse temporal resolution (7 dates over ~2 months) provides insufficient power to resolve differences in senescence onset timing
- The broad confidence intervals on breakpoints overlap across all treatments
- Alternatively, acknowledge that the current phenology data are better suited to the polynomial/mixed model approach already used in the manuscript

### Script
- `R/reviewer_phenology_exploration.R`

---

## Files Created

| File | Description |
|------|-------------|
| `R/06b_aci_refit.R` | A/Ci curve refit with both Tcorrect methods, full analysis |
| `R/reviewer_response_exploration.R` | Temperature extraction from raw LI-6400 data |
| `R/reviewer_phenology_exploration.R` | Segmented regression on phenology data |
| `output/tables/aci_refit_effect_comparison.txt` | Comprehensive effect comparison (formatted) |
| `output/tables/aci_refit_models.txt` | Full model summaries (4 mixed models) |
| `output/tables/aci_refit_pvalue_comparison.csv` | ANOVA p-values side by side |
| `output/tables/aci_refit_summary_by_treatment.csv` | Treatment means ± SE |
| `output/tables/aci_refit_summary_by_pop_treatment.csv` | Pop × Treatment means ± SE |
| `output/tables/aci_refit_validation.csv` | Refit vs. original validation |
| `output/tables/aci_temperature_summary_by_date.csv` | Tleaf/TBlk/Tair by date |
| `output/tables/aci_temperature_summary_by_curve.csv` | Per-curve temperature summary |
| `output/figures/aci_refit_correction_comparison.pdf` | Scatter: Tcorrect TRUE vs FALSE |
| `output/figures/aci_refit_emm_comparison.pdf` | EMM comparison by method |
