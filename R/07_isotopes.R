# 07_isotopes.R
source("R/00_setup.R")

# Read isotope data
d13c <- read_csv("data/raw/Ecotypes2017_Drying_d13C.csv") %>%
  mutate(
    Population = factor(Population, 
                        levels = POPULATION_LEVELS, 
                        labels = POPULATION_LABELS),
    `Treatment Name` = factor(`Treatment Name`, 
                              levels = TREATMENT_LEVELS)
  )

# Function to calculate difference between two populations
diff_boot <- function(data, indices, pop1, pop2) {
  d <- data[indices, ]
  mean1 <- mean(d$d13C[d$Population == pop1])
  mean2 <- mean(d$d13C[d$Population == pop2])
  return(mean1 - mean2)
}


# Define your custom bootstrap function (mean calculation)
bootfun <- function(data, indices) {
  d <- data[indices]
  return(mean(d))
}

# Analysis by Population only
results_pop <- d13c %>%
  group_by(Population) %>%
  dplyr::summarise(
    mean = mean(d13C, na.rm = TRUE),
    boot = list(boot(d13C, bootfun, R = 1000)),
    ci_lower = quantile(boot[[1]]$t, 0.025),
    ci_upper = quantile(boot[[1]]$t, 0.975),
    .groups = 'drop'
  )

# Analysis by Treatment only
results_treat <- d13c %>%
  group_by(`Treatment Name`) %>%
  summarise(
    mean = mean(d13C, na.rm = TRUE),
    boot = list(boot(d13C, bootfun, R = 1000)),
    ci_lower = quantile(boot[[1]]$t, 0.025),
    ci_upper = quantile(boot[[1]]$t, 0.975),
    .groups = 'drop'
  )

# Extract bootstrap samples into long format dataframes
boot_samples_pop <- map2_dfr(results_pop$boot, results_pop$Population, ~{
  tibble(
    d13C = as.vector(.x$t),
    Population = .y
  )
})

boot_samples_treat <- map2_dfr(results_treat$boot, results_treat$`Treatment Name`, ~{
  tibble(
    d13C = as.vector(.x$t),
    `Treatment Name` = .y
  )
})

# Create Population Plot
p1a <- ggplot(boot_samples_pop, aes(x = d13C, 
                                    y = fct_rev(`Population`)))+
  geom_density_ridges_gradient(aes(fill = after_stat(x)), 
                               scale = 1,
                               alpha = 0.6, 
                               bandwidth = 0.7,
                               rel_min_height = 0.01,
                               height = 0.3) +
  geom_point(data = d13c, 
             aes(x = d13C, 
                 color = Population), 
             size = 3, 
             alpha = 0.6,
             shape = 21,
             stroke = 1) +
  theme_classic() +
  labs(x = expression(delta^13*"C (\u2030)"),
       y = "Population", 
       fill = "iWUE", 
       color = "Population") +
  scale_fill_viridis_c(
    direction = 1,
    breaks = seq(-28, -25, 1),
    limits = c(-29, -24),
    labels = c("Lower", "", "", "Higher"),
    oob = scales::squish,  # Squishes values outside limits
    rescaler = function(x, from) scales::rescale_mid(x, from = from, to = c(0, 1), mid = -26.5)
  )+
  scale_color_manual(values = custom_colors_pop) +
  theme(legend.position = "right") +
  # Add guide order to ensure Population comes first
  guides(color = guide_legend(order = 1),
         fill = guide_colorbar(order = 2))

# Create Treatment Plot
p2a <- ggplot(boot_samples_treat, aes(x = d13C, 
                                      y = fct_rev(`Treatment Name`)))+
  geom_density_ridges_gradient(aes(fill = after_stat(x)), 
                               scale = 1,
                               alpha = 0.6, 
                               bandwidth = 0.7,
                               rel_min_height = 0.01,
                               height = 0.3) +
  geom_point(data = d13c, 
             aes(x = d13C, 
                 color = `Treatment Name`), 
             size = 3, 
             alpha = 0.6,
             shape = 21,
             stroke = 1) +
  theme_classic() +
  labs(x = expression(delta^13*"C (\u2030)"),
       y = "Treatment", 
       fill = "iWUE", 
       color = "Treatment") +
  scale_fill_viridis_c(
    direction = 1,
    breaks = seq(-28, -25, 1),
    limits = c(-29, -24),
    labels = c("Lower", "", "", "Higher"),
    oob = scales::squish,  # Squishes values outside limits
    rescaler = function(x, from) scales::rescale_mid(x, from = from, to = c(0, 1), mid = -26.5)
  )+
  scale_color_manual(values = TREATMENT_COLORS) +
  theme(legend.position = "right") +
  # Add guide order to ensure Treatment comes last
  guides(color = guide_legend(order = 3),
         fill = guide_colorbar(order = 2))

# Combine plots using patchwork
wplot <- (p1a + theme(legend.position = "right")) / 
  (p2a + theme(legend.position = "right")) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left"
  )

# Print the combined plot
print(wplot)
ggsave("output/figures/isotopes.pdf", wplot, width = 6, height = 6, units="in", dpi=600)
ggsave("output/figures/isotopes.pdf", wplot, width = 6, height = 6, units = "in", dpi = 600, device = cairo_pdf)



# Calculate and save summary statistics
d13c_summary <- d13c %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_d13c = mean(d13C, na.rm = TRUE),
    se_d13c = sd(d13C, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write_csv(d13c_summary, "output/tables/d13c_summary.csv")

# Calculate and save summary statistics
d13c_summary <- d13c %>%
  group_by(Population, `Treatment Name`) %>%
  summarize(
    mean_d13c = mean(d13C, na.rm = TRUE),
    se_d13c = sd(d13C, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write_csv(d13c_summary, "output/tables/d13c_summary.csv")
