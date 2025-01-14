# 10_model_diagnostics.R
# Model diagnostic checks for all fitted models
source("R/00_setup.R")

# Load fitted models from previous scripts
# Note: requires running scripts 01-07 first
source("R/02_phenology.R")     # contains phenology models
source("R/03_canopy.R")        # contains ndvi/lai models

# Function to run cross-validation for any mixed model
cv_mixed_model <- function(formula, data, response_var, folds = 10) {
  set.seed(123)
  data$folds <- sample(rep(1:folds, length.out = nrow(data)))
  mse_list <- vector("numeric", length = folds)
  
  for (i in 1:folds) {
    train_data <- data[data$folds != i, ]
    test_data <- data[data$folds == i, ]
    model <- lmer(formula, data = train_data)
    predictions <- predict(model, newdata = test_data, allow.new.levels = TRUE)
    mse_list[i] <- mean((test_data[[response_var]] - predictions)^2)
  }
  
  mean(mse_list)
}

# Function to create diagnostic plots for a model
create_diagnostic_plots <- function(model, title) {
  par(mfrow = c(2, 2))
  
  # Residuals vs Fitted
  plot(fitted(model), resid(model),
       xlab = "Fitted values", ylab = "Residuals",
       main = paste(title, "- Residuals vs Fitted"))
  abline(h = 0, lty = 2)
  
  # Q-Q plot
  qqnorm(resid(model), main = paste(title, "- Normal Q-Q"))
  qqline(resid(model))
  
  # Scale-Location plot
  plot(fitted(model), sqrt(abs(resid(model))),
       xlab = "Fitted values", ylab = "√|Residuals|",
       main = paste(title, "- Scale-Location"))
  
  # Residuals vs Time (if DOY exists)
  if("DOY" %in% names(model@frame)) {
    plot(model@frame$DOY, resid(model),
         xlab = "Day of Year", ylab = "Residuals",
         main = paste(title, "- Residuals vs Time"))
    abline(h = 0, lty = 2)
  }
  
  par(mfrow = c(1, 1))
}

# Run diagnostics and save results
sink("output/tables/model_diagnostics.txt")

# Phenology models
cat("\n=== Phenology Models ===\n")
cat("\nTiller Total Green Length Model:\n")
print(summary(tot_gl_quadratic))
cat("\nAIC =", AIC(tot_gl_quadratic))
cat("\nBIC =", BIC(tot_gl_quadratic))
cat("\nRandom effects variance:", capture.output(VarCorr(tot_gl_quadratic)))

# NDVI/LAI models
cat("\n\n=== NDVI/LAI Models ===\n")
cat("\nNDVI Quadratic Model:\n")
print(summary(ndvi_quadratic))
cat("\nAIC =", AIC(ndvi_quadratic))
cat("\nBIC =", BIC(ndvi_quadratic))
cat("\nRandom effects variance:", capture.output(VarCorr(ndvi_quadratic)))

cat("\nLAI Quadratic Model:\n")
print(summary(lai_quadratic))
cat("\nAIC =", AIC(lai_quadratic))
cat("\nBIC =", BIC(lai_quadratic))
cat("\nRandom effects variance:", capture.output(VarCorr(lai_quadratic)))

sink()

# Create diagnostic plots
pdf("output/figures/model_diagnostics.pdf", width = 10, height = 8)
create_diagnostic_plots(tot_gl_quadratic, "Total Green Length")
create_diagnostic_plots(ndvi_quadratic, "NDVI")
create_diagnostic_plots(lai_quadratic, "LAI")
dev.off()

# Calculate and save cross-validation results
cv_results <- data.frame(
  Model = c("Total Green Length", "NDVI", "LAI"),
  CV_MSE = c(
    cv_mixed_model(tot_gl_quadratic@call$formula, pheno1, "TillerTotalGreenLength"),
    cv_mixed_model(ndvi_quadratic@call$formula, ndvi_data, "NDVI"),
    cv_mixed_model(lai_quadratic@call$formula, lai_data, "LAI")
  )
)

write_csv(cv_results, "output/tables/cross_validation_results.csv")

# Create model comparison table
model_comparison <- data.frame(
  Model = c("Total Green Length", "NDVI", "LAI"),
  AIC = c(AIC(tot_gl_quadratic), AIC(ndvi_quadratic), AIC(lai_quadratic)),
  BIC = c(BIC(tot_gl_quadratic), BIC(ndvi_quadratic), BIC(lai_quadratic)),
  R2_marginal = c(
    r.squaredGLMM(tot_gl_quadratic)[1],
    r.squaredGLMM(ndvi_quadratic)[1],
    r.squaredGLMM(lai_quadratic)[1]
  ),
  R2_conditional = c(
    r.squaredGLMM(tot_gl_quadratic)[2],
    r.squaredGLMM(ndvi_quadratic)[2],
    r.squaredGLMM(lai_quadratic)[2]
  )
)

write_csv(model_comparison, "output/tables/model_comparison.csv")

# Create and save residual correlation plot
residual_correlations <- data.frame(
  Time = c(pheno1$DOY, ndvi_data$DOY, lai_data$DOY),
  Residuals = c(resid(tot_gl_quadratic), resid(ndvi_quadratic), resid(lai_quadratic)),
  Model = factor(rep(c("Total Green Length", "NDVI", "LAI"), 
                     c(nrow(pheno1), nrow(ndvi_data), nrow(lai_data))))
)

residual_plot <- ggplot(residual_correlations, aes(x = Time, y = Residuals, color = Model)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic() +
  labs(x = "Day of Year", y = "Model Residuals") +
  theme(legend.position = "top")

ggsave("output/figures/residual_correlations.pdf", residual_plot, width = 10, height = 6)