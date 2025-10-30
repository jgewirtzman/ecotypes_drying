# Load libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

# ===== LOAD DATA =====

rain <- read.csv('/Users/jongewirtzman/Downloads/edc_data_download (3)/1-hour_data.csv')
par <- read.csv('/Users/jongewirtzman/Downloads/edc_data_download (5)/1-hour_data.csv')

# ===== RAIN DATA ANALYSIS =====

# Calculate percentage of hours with rain (excluding NAs)
pct_rain <- sum(rain$rain > 0, na.rm = TRUE) / sum(!is.na(rain$rain)) * 100
cat("Percentage of hours with rain:", pct_rain, "%\n")

# Filter for study period (July 17 - Sept 11)
rain_subset <- rain[rain$date >= "2017-07-17" & rain$date <= "2017-09-11", ]
pct_rain_period <- sum(rain_subset$rain > 0, na.rm = TRUE) / sum(!is.na(rain_subset$rain)) * 100
cat("Percentage of hours with rain (July 17 - Sept 11):", pct_rain_period, "%\n")

# Compare rain and rain2 measurements
comparison <- rain[!is.na(rain$rain) & !is.na(rain$rain2), ]
cat("\nCorrelation between rain and rain2:", cor(comparison$rain, comparison$rain2), "\n")

# Linear regression model
model <- lm(rain2 ~ rain, data = comparison)
summary(model)

# Fill missing rain2 values using model predictions
rain$rain2_filled <- rain$rain2
missing_rain2 <- is.na(rain$rain2) & !is.na(rain$rain)
predicted_values <- predict(model, newdata = rain[missing_rain2, ])
rain$rain2_filled[missing_rain2] <- ifelse(
  rain$rain[missing_rain2] > 0,
  pmax(0, predicted_values),
  0
)

# ===== DAILY RAINFALL AGGREGATION =====

daily_rain <- rain %>%
  group_by(date) %>%
  summarise(
    daily_rain1 = sum(rain, na.rm = TRUE),
    daily_rain2 = sum(rain2_filled, na.rm = TRUE)
  )

daily_rain_subset <- daily_rain[daily_rain$date >= "2017-07-17" & 
                                  daily_rain$date <= "2017-09-11", ]

# Count days with >1mm and >5mm rain
days_over_1mm <- sum(daily_rain_subset$daily_rain2 > 1)
days_over_5mm <- sum(daily_rain_subset$daily_rain2 > 5)
total_days <- nrow(daily_rain_subset)

cat("\nPeriod: July 17 - Sept 11, 2017\n")
cat("Total days:", total_days, "\n")
cat("Days with >1mm rain:", days_over_1mm, "(", round(days_over_1mm/total_days*100, 1), "%)\n")
cat("Days with >5mm rain:", days_over_5mm, "(", round(days_over_5mm/total_days*100, 1), "%)\n")

# ===== 2-WEEK PERIOD ANALYSIS =====

daily_rain_subset$date <- as.Date(daily_rain_subset$date)
daily_rain_subset$period <- NA
start_date <- as.Date("2017-07-17")
period_breaks <- seq(start_date, as.Date("2017-09-11"), by = "14 days")
period_breaks <- c(period_breaks, as.Date("2017-09-12"))

for (i in 1:(length(period_breaks) - 1)) {
  mask <- daily_rain_subset$date >= period_breaks[i] & 
    daily_rain_subset$date < period_breaks[i + 1]
  daily_rain_subset$period[mask] <- i
}

period_summary <- daily_rain_subset %>%
  group_by(period) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    total_days = n(),
    days_over_5mm = sum(daily_rain2 > 5),
    pct_over_5mm = round(sum(daily_rain2 > 5) / n() * 100, 1)
  )

print(period_summary)

# ===== PAR DATA PROCESSING =====

par$date <- as.Date(par$date)
daily_par <- par %>%
  group_by(date) %>%
  summarise(daily_par = sum(quantum * 3600 / 1000000, na.rm = TRUE))

# ===== SOIL MOISTURE CALIBRATION AND PROCESSING =====

# Read and process calibration data
curve <- read_csv("data/raw/gwc_calibration.csv")

# Extract empty pot weights
empty_pot_data <- curve[grep("Empty Pot", curve$Date), ]
empty_pot_weights <- setNames(empty_pot_data$`Mass (g)`, empty_pot_data$`Pot #`)

# Process calibration data
curve <- curve %>%
  filter(!grepl("Empty Pot", Date)) %>%
  mutate(`Wet Soil Mass` = `Mass (g)` - empty_pot_weights[`Pot #`]) %>%
  group_by(`Pot #`) %>%
  mutate(
    `Dry Soil Mass` = min(`Wet Soil Mass`),
    `Water Content` = `Wet Soil Mass` - `Dry Soil Mass`,
    `Gravimetric Water Content` = `Water Content` / `Dry Soil Mass`
  ) %>%
  ungroup() %>%
  na.omit()

# Fit calibration model
calibration_model <- lm(`Gravimetric Water Content` ~ poly(`PER (mean)`, 3), data = curve)

# Read and process soil moisture data
soil_gwc <- read_csv("data/raw/Ecotypes2017_Drying_Soil Moisture.csv")
names(soil_gwc)[names(soil_gwc) == "PER"] <- "PER (mean)"

# Predict GWC and prepare data
soil_gwc$GWC <- predict(calibration_model, newdata = soil_gwc)
soil_gwc$date <- as.Date(paste0("2017-", soil_gwc$DOY), format = "%Y-%j")

# Calculate daily mean GWC by treatment
daily_gwc_treatment <- soil_gwc %>%
  group_by(date, `Treatment Name`) %>%
  summarize(mean_gwc = mean(GWC * 100, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(`Treatment Name`))

# ===== PLOTTING DATA PREPARATION =====

# Prepare plotting data
daily_rain_plot <- daily_rain %>%
  filter(date >= "2017-07-01" & date <= "2017-09-11") %>%
  mutate(date = as.Date(date))

daily_par_plot <- daily_par %>%
  filter(date >= "2017-07-01" & date <= "2017-09-11")

# Define exclusion periods
exclusion_periods_df <- data.frame(
  xmin = as.Date(c("2017-07-01", "2017-07-28", "2017-08-27")),
  xmax = as.Date(c("2017-07-17", "2017-08-11", "2017-09-11"))
)

# Mark days with >5mm rain during exclusion periods
daily_rain_plot$is_exclusion_rain <- FALSE
for (i in 1:nrow(exclusion_periods_df)) {
  daily_rain_plot$is_exclusion_rain <- daily_rain_plot$is_exclusion_rain | 
    (daily_rain_plot$date >= exclusion_periods_df$xmin[i] & 
       daily_rain_plot$date <= exclusion_periods_df$xmax[i] & 
       daily_rain_plot$daily_rain1 > 5)
}

# Calculate scaling factor for dual y-axis
rainfall_max <- max(daily_rain_plot$daily_rain1, na.rm = TRUE)
gwc_max <- max(daily_gwc_treatment$mean_gwc, na.rm = TRUE)
scale_factor <- rainfall_max / gwc_max

# ===== CREATE PLOTS =====

# Combined rainfall and soil moisture plot
combined_rain_soil <- ggplot() +
  geom_rect(data = exclusion_periods_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "gray80", alpha = 0.3) +
  geom_col(data = daily_rain_plot,
           aes(x = date, y = daily_rain1, fill = is_exclusion_rain),
           alpha = 0.7, width = 1) +
  scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "grey50"), guide = "none") +
  geom_line(data = daily_gwc_treatment,
            aes(x = date, y = mean_gwc * scale_factor, color = `Treatment Name`),
            size = 1.2) +
  geom_point(data = daily_gwc_treatment,
             aes(x = date, y = mean_gwc * scale_factor, color = `Treatment Name`),
             size = 2) +
  scale_x_date(limits = c(as.Date("2017-07-01"), as.Date("2017-09-11")),
               date_breaks = "1 week", date_labels = "%b %d") +
  scale_y_continuous(
    name = "Daily Rainfall (mm)",
    sec.axis = sec_axis(~./scale_factor, name = "Soil Moisture (GWC %)")
  ) +
  scale_color_manual(
    values = c("Wet" = "#1f77b4", "Dry" = "#8c564b", "Deep" = "#d62728"),
    name = "Treatment"
  ) +
  labs(
    # title = "Daily Rainfall, Soil Moisture, and PAR: July 1 - Sept 11, 2017",
    # subtitle = "Grey bars = days with >5mm rain during exclusion periods",
    x = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    #axis.title.y.right = element_text(color = "darkgreen"),
    #axis.text.y.right = element_text(color = "darkgreen"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

# PAR plot
# Create rectangles for covered days only
covered_days_df <- daily_rain_plot %>%
  filter(is_exclusion_rain) %>%
  dplyr::select(date)

par_plot <- ggplot() +
  geom_rect(data = covered_days_df,
            aes(xmin = date - 0.5, xmax = date + 0.5, ymin = -Inf, ymax = Inf),
            fill = "gray80", alpha = 0.3) +
  geom_line(data = daily_par_plot, aes(x = date, y = daily_par),
            color = "orange", size = 1) +
  geom_point(data = daily_par_plot, aes(x = date, y = daily_par),
             color = "orange", size = 1.5) +
  scale_x_date(limits = c(as.Date("2017-07-01"), as.Date("2017-09-11")),
               date_breaks = "1 week", date_labels = "%b %d") +
  labs(y = "Daily PAR (mol/m²/day)", x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

# Timeline plot
timeline_data <- data.frame(
  date = seq(as.Date("2017-07-01"), as.Date("2017-09-12"), by = "day")
)
timeline_data$covered <- timeline_data$date %in% 
  daily_rain_plot$date[daily_rain_plot$is_exclusion_rain]

timeline_plot <- ggplot(timeline_data, aes(x = date, y = 0.5)) +
  geom_tile(aes(fill = covered), height = 0.8, width = 1) +
  scale_fill_manual(
    values = c("FALSE" = "lightgray", "TRUE" = "darkgray"),
    labels = c("Uncovered", "Covered"),
    name = ""
  ) +
  scale_x_date(limits = c(as.Date("2017-07-01"), as.Date("2017-09-12")),
               date_breaks = "1 week", date_labels = "%b %d") +
  scale_y_continuous(limits = c(0, 1), breaks = 0.5, labels = "Tarp\nStatus") +
  labs(x = "Date", y = "") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

# Combine all plots
combined_plots <- combined_rain_soil / par_plot / timeline_plot + 
  plot_layout(heights = c(4, 2, 1))

print(combined_plots)

# ===== PAR ANALYSIS FOR COVERED DAYS =====

covered_dates <- daily_rain_plot$date[daily_rain_plot$is_exclusion_rain]

par_summary <- daily_par_plot %>%
  summarise(
    total_par = sum(daily_par, na.rm = TRUE),
    par_covered = sum(daily_par[date %in% covered_dates], na.rm = TRUE),
    pct_covered = par_covered / total_par * 100
  )

cat("\n=== PAR Analysis for Covered Days ===\n")
cat(sprintf("Total PAR (Jul 1 - Sep 11): %.2f mol/m²\n", par_summary$total_par))
cat(sprintf("PAR on covered days: %.2f mol/m² (%.1f%%)\n", 
            par_summary$par_covered, par_summary$pct_covered))
cat(sprintf("Covered days: %d out of 73 (%.1f%%)\n", 
            length(covered_dates), length(covered_dates)/73*100))





# ===== CUMULATIVE PAR BY TREATMENT =====

# Create full date range
all_dates <- data.frame(
  date = seq(as.Date("2017-07-01"), as.Date("2017-09-11"), by = "day")
)

# Merge with PAR data
all_dates <- left_join(all_dates, daily_par_plot, by = "date")
all_dates$daily_par[is.na(all_dates$daily_par)] <- 0

# Merge with covered status
all_dates <- left_join(all_dates, 
                       daily_rain_plot %>% dplyr::select(date, is_exclusion_rain),
                       by = "date")
all_dates$is_exclusion_rain[is.na(all_dates$is_exclusion_rain)] <- FALSE

# Calculate PAR received by each treatment
# Wet always gets full PAR; Dry and Deep get 60% when covered
all_dates$Wet_par <- all_dates$daily_par  
all_dates$Dry_par <- ifelse(all_dates$is_exclusion_rain, 
                            all_dates$daily_par * 0.6, 
                            all_dates$daily_par)
all_dates$Deep_par <- ifelse(all_dates$is_exclusion_rain, 
                             all_dates$daily_par * 0.6, 
                             all_dates$daily_par)

# Calculate cumulative PAR
all_dates$Wet_cumpar <- cumsum(all_dates$Wet_par)
all_dates$Dry_cumpar <- cumsum(all_dates$Dry_par)
all_dates$Deep_cumpar <- cumsum(all_dates$Deep_par)

# Reshape for plotting
cumpar_long <- all_dates %>%
  dplyr::select(date, Wet_cumpar, Dry_cumpar, Deep_cumpar) %>%
  tidyr::pivot_longer(cols = c(Wet_cumpar, Dry_cumpar, Deep_cumpar),
                      names_to = "Treatment",
                      values_to = "cumulative_par") %>%
  mutate(Treatment = gsub("_cumpar", "", Treatment))
cumpar_long <- cumpar_long %>%
  mutate(Treatment = ifelse(Treatment %in% c("Dry", "Deep"), "Dry/Deep", Treatment)) %>%
  distinct()

# Create cumulative PAR plot
cumpar_plot <- ggplot(cumpar_long, aes(x = date, y = cumulative_par, color = Treatment)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_date(limits = c(as.Date("2017-07-01"), as.Date("2017-09-11")),
               date_breaks = "1 week", date_labels = "%b %d") +
  scale_color_manual(
    values = c("Wet" = "#1f77b4", "Dry/Deep" = "#d62728"),
    name = "Treatment"
  ) +
  labs(
    y = "Cumulative PAR (mol/m²)"
  ) +
  theme_classic() +
  labs(x = "Date")+
  theme(
    legend.position = "right",
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank()
  )

print(cumpar_plot)

# Print summary statistics
cat("\n=== Cumulative PAR by Treatment (Jul 1 - Sep 11) ===\n")
for (trt in c("Wet", "Dry", "Deep")) {
  total <- tail(all_dates[[paste0(trt, "_cumpar")]], 1)
  covered_days <- sum(all_dates$is_exclusion_rain)
  par_lost <- sum(all_dates$daily_par[all_dates$is_exclusion_rain] * 0.4)  # 40% reduction
  
  if (trt == "Wet") {
    cat(sprintf("%s: %.2f mol/m² (no coverage)\n", trt, total))
  } else {
    cat(sprintf("%s: %.2f mol/m² (lost %.2f mol/m² over %d covered days, 60%% transmission)\n", 
                trt, total, par_lost, covered_days))
  }
}

# Update combined plots to include cumulative PAR
combined_plots_with_cumpar <- combined_rain_soil / par_plot / timeline_plot / cumpar_plot + 
  plot_layout(heights = c(3, 2, 1, 2))

print(combined_plots_with_cumpar)



# ===== CUMULATIVE LAI BY TREATMENT =====

# Read NDVI data and calculate LAI
ndvi_timeseries <- read_csv("data/raw/Ecotypes2017_Drying_NDVI.csv") %>%
  select_if(~ !all(is.na(.))) %>%
  filter(rowSums(is.na(.)) < ncol(.)) %>%
  mutate(
    LAI = 0.03 * exp(7.65 * NDVI),
    date = as.Date(paste0("2017-", DOY), format = "%Y-%j")
  )

# Calculate daily mean LAI by treatment
daily_lai_treatment <- ndvi_timeseries %>%
  group_by(date, `Treatment Name`) %>%
  summarize(mean_lai = mean(LAI, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(`Treatment Name`))

# Calculate cumulative LAI (sum of daily LAI values)
cumulative_lai_treatment <- daily_lai_treatment %>%
  arrange(`Treatment Name`, date) %>%
  group_by(`Treatment Name`) %>%
  mutate(cumulative_lai = cumsum(mean_lai)) %>%
  ungroup() %>%
  mutate(`Treatment Name` = factor(`Treatment Name`, levels = c("Wet", "Dry", "Deep")))

# Create cumulative LAI plot (WITH date on x-axis)
cumlai_plot <- ggplot(cumulative_lai_treatment, 
                      aes(x = date, y = cumulative_lai, color = `Treatment Name`)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_date(limits = c(as.Date("2017-07-01"), as.Date("2017-09-11")),
               date_breaks = "1 week", date_labels = "%b %d") +
  scale_color_manual(
    values = c("Wet" = "#1f77b4", "Dry" = "#8c564b", "Deep" = "#d62728"),
    name = "Treatment"
  ) +
  labs(
    y = "Cumulative LAI",
    x = "Date"
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  )

print(cumlai_plot)

# Print summary statistics
cat("\n=== Cumulative LAI by Treatment (Jul 1 - Sep 11) ===\n")
lai_summary <- cumulative_lai_treatment %>%
  group_by(`Treatment Name`) %>%
  summarize(
    final_cumulative_lai = last(cumulative_lai),
    .groups = "drop"
  )
print(lai_summary)

# Update cumpar_plot to remove x-axis labels
cumpar_plot <- ggplot(cumpar_long, aes(x = date, y = cumulative_par, color = Treatment)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_date(limits = c(as.Date("2017-07-01"), as.Date("2017-09-11")),
               date_breaks = "1 week", date_labels = "%b %d") +
  scale_color_manual(
    values = c("Wet" = "#1f77b4", "Dry/Deep" = "#d62728"),
    name = "Treatment"
  ) +
  labs(
    y = "Cumulative PAR (mol/m²)",
    x = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

# Update combined plots to include both cumulative PAR and LAI
combined_plots_full <- combined_rain_soil / par_plot / timeline_plot / 
  cumpar_plot / cumlai_plot + 
  plot_layout(heights = c(3, 2, 1, 2, 2))

print(combined_plots_full)
