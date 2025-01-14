library(ggplot2)
library(sf)
library(cowplot)
library(elevatr)
library(ggspatial)
library(tigris)
library(dplyr)
library(ggnewscale)  # Added this library

# Custom colors for site points
custom_colors <- c("Coldfoot" = "#3E6347",  
                   "Toolik" = "#90C49B", 
                   "Sagwon" = "#79A9C8")

# Define the locations
locations <- data.frame(
  name = c("Toolik", "Sagwon", "Coldfoot"),
  lat = c(68.6269, 69.4245, 67.2528),
  lon = c(-149.5975, -148.6946, -150.1847)
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Create a polygon for our area of interest
aoi <- st_polygon(list(rbind(
  c(-152, 67),
  c(-152, 70),
  c(-147, 70),
  c(-147, 67),
  c(-152, 67)
))) %>%
  st_sfc() %>%
  st_sf() %>%
  st_set_crs(4326)

# Get elevation data
elevation <- get_elev_raster(aoi, z = 6, src = "aws")

# Convert elevation raster to data frame for ggplot
elevation_df <- raster::as.data.frame(elevation, xy = TRUE)
names(elevation_df) <- c("x", "y", "elevation")

# Get Alaska boundary data
alaska_sf <- states(cb = TRUE) %>%
  filter(STUSPS == "AK")
alaska_counties <- counties(state = "AK", cb = TRUE)

# Create main map
main_map <- ggplot() +
  # First layer for elevation
  geom_raster(data = elevation_df, 
              aes(x = x, y = y, fill = elevation)) +
  scale_fill_gradientn(
    colors = c("#2F4F3F", "#627F63", "#9FA78C", "#D9CC9B", "#B29866", "#8C6D45", "#FFFFFF"),
    name = "Elevation (m)",
    breaks = c(0, 500, 1000, 1500, 2000),
    limits = c(0, 2000)
  ) +
  # Add new scale for points
  new_scale_fill() +
  geom_sf(data = locations, aes(fill = name), size = 4, shape = 21, color = "black", stroke = 1.5) +
  scale_fill_manual(
    values = custom_colors,
    guide = "none"
  ) +
  geom_sf_text(data = locations, aes(label = name),
               nudge_x = 0.3, nudge_y = 0.2, size = 5, fontface = "bold", color = "white", stroke = 0.3) +
  coord_sf(xlim = c(-152, -147), 
           ylim = c(67, 70)) +
  annotation_scale(
    location = "bl",
    bar_cols = c("black", "white"),
    text_family = "sans",
    width_hint = 0.25
  ) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_minimal(),
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in")
  ) +
  scale_x_continuous(
    breaks = c(-152, -151, -150, -149, -148, -147),
    labels = c("152°W", "151°W", "150°W", "149°W", "148°W", "147°W")
  ) +
  scale_y_continuous(
    breaks = seq(67, 70, by = 0.5),
    labels = function(x) paste0(x, "°N")
  ) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 8),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Create inset map using Albers projection for Alaska
alaska_albers <- st_transform(alaska_sf, 
                              crs = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154")
alaska_counties_albers <- st_transform(alaska_counties,
                                       crs = st_crs(alaska_albers))
study_box <- st_as_sfc(st_bbox(c(xmin = -152, xmax = -147, 
                                 ymin = 67, ymax = 70))) %>%
  st_set_crs(4326) %>%
  st_transform(st_crs(alaska_albers))

inset <- ggplot() +
  geom_sf(data = alaska_counties_albers, fill = "white", color = "black", size = 0.1) +
  geom_sf(data = alaska_albers, fill = NA, color = "black", size = 0.2) +
  geom_sf(data = study_box, fill = NA, color = "#B24444", linewidth = 1) +
  theme_void()

# Combine the maps
final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset, x = 0.65, y = 0.65, width = 0.3, height = 0.3)

# Display the map
final_map
