# Load libraries
# install.packages(c("stats19", "dplyr", "sf", "ggplot2", "viridis", "readr"))

library(stats19)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(readr)

# Set data preparation parameters
years_to_download <- 2020:2024
grid_cellsize <- 2500
max_candidate_sites <- 80
set.seed(42)

# Define Birmingham study area
bbox_lonlat <- st_as_sfc(st_bbox(c(xmin = -2.20, xmax = -1.60, ymin = 52.30, ymax = 52.65), crs = 4326))

# Download collision data for one specified year
download_stats19_year <- function(year) {
  cat("Reading collision data for", year, "\n")
  get_stats19(year = year, type = "collision", output_format = "sf") %>% mutate(year_downloaded = year)
}

# Download and combine all years
collisions_all <- bind_rows(lapply(years_to_download, download_stats19_year))
total_collisions <- nrow(collisions_all)

cat("Total collisions across all years:", total_collisions, "\n")

# Filter collisions to Birmingham area
bbox_proj <- st_transform(bbox_lonlat, st_crs(collisions_all))
collisions_all <- st_filter(collisions_all, bbox_proj)

filtered_collisions <- nrow(collisions_all)

cat("Collisions after bounding box filter:", filtered_collisions, "\n")

if (filtered_collisions == 0) stop("No collisions found inside the selected bounding box.")

# Convert data to metre-based CRS
collisions_all <- st_transform(collisions_all, 27700)

# Find severity column
severity_column <- if ("collision_severity" %in% names(collisions_all)) {
  "collision_severity"
} else if ("enhanced_severity_collision" %in% names(collisions_all)) {
  "enhanced_severity_collision"
} else {
  stop("No recognised severity column found.")
}

# Add severity weights
collisions_all <- collisions_all %>%
  mutate(
    severity = .data[[severity_column]],
    sev_weight = case_when(
      severity %in% c("Slight", "slight", 3, "3") ~ 1,
      severity %in% c("Serious", "serious", 2, "2") ~ 2,
      severity %in% c("Fatal", "fatal", 1, "1") ~ 3,
      TRUE ~ 1
    )
  )

# Create demand grid
grid_geom <- st_make_grid(collisions_all, cellsize = grid_cellsize)

grid_sf <- st_sf(
  cell_id = seq_along(grid_geom),
  geometry = grid_geom,
  crs = st_crs(collisions_all)
)

# Assign collisions to grid cells
collisions_joined <- st_join(collisions_all, grid_sf, left = FALSE)

if (nrow(collisions_joined) == 0) stop("No collisions were assigned to grid cells.")

# Aggregate deterministic demand
demand_points_det <- collisions_joined %>%
  group_by(cell_id) %>%
  summarise(accidents = n(), weight = sum(sev_weight), .groups = "drop")

# Create demand centroids
demand_centroids_det <- st_centroid(demand_points_det)

cat("Deterministic demand cells:", nrow(demand_points_det), "\n")

# Aggregate yearly demand
demand_by_year <- collisions_joined %>%
  st_drop_geometry() %>%
  group_by(year_downloaded, cell_id) %>%
  summarise(accidents = n(), weight = sum(sev_weight), .groups = "drop")

# Create scenario matrix
cell_ids <- demand_points_det$cell_id
scenario_years <- sort(unique(demand_by_year$year_downloaded))

scenario_matrix <- matrix(0, nrow = length(cell_ids), ncol = length(scenario_years), dimnames = list(cell_ids, scenario_years))

# Fill yearly demand scenarios
for (s in seq_along(scenario_years)) {
  yr <- scenario_years[s]
  temp <- demand_by_year %>% filter(year_downloaded == yr)
  scenario_matrix[match(temp$cell_id, cell_ids), s] <- temp$weight
}

# Keep cells with demand
keep_rows <- rowSums(scenario_matrix) > 0

scenario_matrix <- scenario_matrix[keep_rows, , drop = FALSE]
demand_points_det <- demand_points_det[keep_rows, ]
demand_centroids_det <- st_centroid(demand_points_det)

# Store deterministic weights
demand_weights_det <- demand_points_det$weight

# Check row alignment
stopifnot(all(rownames(scenario_matrix) == as.character(demand_points_det$cell_id)))

# Store stochastic scenarios
scenario_weights_list <- lapply(seq_len(ncol(scenario_matrix)), function(k) scenario_matrix[, k])

# Set equal scenario probabilities
scenario_probs <- rep(1 / length(scenario_weights_list), length(scenario_weights_list))

cat("Demand cells after removing zero rows:", nrow(demand_points_det), "\n")
cat("Number of stochastic scenarios:", length(scenario_weights_list), "\n")

# Build candidate site pool
candidate_pool <- demand_centroids_det %>%
  mutate(candidate_pool_id = row_number(), sample_prob = weight / sum(weight))

# Sample candidate sites
if (nrow(candidate_pool) > max_candidate_sites) {
  candidate_sites <- candidate_pool %>% slice_sample(n = max_candidate_sites, weight_by = sample_prob)
} else {
  candidate_sites <- candidate_pool
}

# Add candidate site IDs
candidate_sites <- candidate_sites %>%
  select(-sample_prob) %>%
  mutate(site_id = row_number())

cat("Candidate sites:", nrow(candidate_sites), "\n")

# Extract coordinates
demand_coords <- st_coordinates(demand_centroids_det)
site_coords <- st_coordinates(candidate_sites)

# Calculate Euclidean distances
distance_matrix <- sqrt(
  outer(demand_coords[, 1], site_coords[, 1], "-")^2 +
    outer(demand_coords[, 2], site_coords[, 2], "-")^2
)

cat("Distance matrix size:", paste(dim(distance_matrix), collapse = " x "), "\n")

# Validate input dimensions
stopifnot(nrow(distance_matrix) == length(demand_weights_det))
stopifnot(nrow(distance_matrix) == nrow(scenario_matrix))
stopifnot(ncol(distance_matrix) == nrow(candidate_sites))
stopifnot(abs(sum(scenario_probs) - 1) < 1e-8)

# Save model inputs
saveRDS(demand_points_det, "demand_points_det.rds")
saveRDS(candidate_sites, "candidate_sites.rds")
saveRDS(demand_centroids_det, "demand_centroids_det.rds")
saveRDS(distance_matrix, "distance_matrix.rds")
saveRDS(demand_weights_det, "demand_weights_det.rds")
saveRDS(scenario_matrix, "scenario_matrix.rds")
saveRDS(scenario_weights_list, "scenario_weights_list.rds")
saveRDS(scenario_probs, "scenario_probs.rds")
saveRDS(grid_sf, "grid_sf.rds")

# Save report CSV files
write_csv(st_drop_geometry(demand_points_det), "demand_points_det.csv")
write_csv(st_drop_geometry(candidate_sites), "candidate_sites.csv")
write_csv(as.data.frame(st_coordinates(demand_centroids_det)), "demand_centroids_det_coords.csv")
write_csv(as.data.frame(st_coordinates(candidate_sites)), "candidate_sites_coords.csv")
write_csv(as.data.frame(distance_matrix), "distance_matrix.csv")
write_csv(data.frame(weight = demand_weights_det), "demand_weights_det.csv")
write_csv(as.data.frame(scenario_matrix), "scenario_matrix.csv")
write_csv(data.frame(probability = scenario_probs), "scenario_probs.csv")

# Create input summary
input_summary <- data.frame(
  item = c(
    "Years used",
    "Total collisions before filtering",
    "Collisions after bounding box filter",
    "Demand points",
    "Candidate sites",
    "Stochastic scenarios",
    "Grid cell size",
    "Candidate site selection",
    "Distance matrix size"
  ),
  value = c(
    paste(min(years_to_download), max(years_to_download), sep = "-"),
    total_collisions,
    filtered_collisions,
    nrow(demand_points_det),
    nrow(candidate_sites),
    length(scenario_weights_list),
    paste0(grid_cellsize, "m x ", grid_cellsize, "m"),
    "Weighted sampling from demand centroids",
    paste(dim(distance_matrix), collapse = " x ")
  )
)

# Save input summary
write_csv(input_summary, "input_summary.csv")

# Create Birmingham boundary for plotting
bbox_map <- st_sf(
  geometry = st_transform(bbox_lonlat, st_crs(demand_points_det))
)

# Get bounding box limits
bbox_limits <- st_bbox(bbox_map)

# Plot demand and candidate sites
plot_demand_sites <- ggplot() +
  geom_sf(
    data = bbox_map,
    fill = NA,
    colour = "grey40",
    linewidth = 0.6
  ) +
  geom_sf(
    data = demand_points_det,
    aes(colour = weight),
    size = 0.9,
    alpha = 0.25
  ) +
  scale_colour_viridis_c(
    option = "plasma",
    trans = "sqrt",
    name = "Demand weight"
  ) +
  geom_sf(
    data = candidate_sites,
    aes(shape = "Candidate ambulance site"),
    size = 3,
    fill = "#56B4E9",
    colour = "black",
    stroke = 1.2
  ) +
  scale_shape_manual(
    name = "Site type",
    values = c("Candidate ambulance site" = 21)
  ) +
  coord_sf(
    xlim = c(bbox_limits["xmin"], bbox_limits["xmax"]),
    ylim = c(bbox_limits["ymin"], bbox_limits["ymax"]),
    expand = FALSE
  ) +
  labs(
    title = "RTC Demand and Candidate Ambulance Sites"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(plot_demand_sites)

# Save full map
ggsave(
  "rtc_demand_candidate_sites.png",
  plot_demand_sites,
  width = 8,
  height = 6,
  dpi = 300
)
print(plot_demand_sites)

cat("Data preparation complete.\n")

