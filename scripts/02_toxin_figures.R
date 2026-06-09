library(tidyverse)
library(lubridate)
library(scales)

#Which toxins are most common?
#Which congeners dominate?
  #Which are rare?

## Detection frequecy by toxin 


# 04_toxin_figures.R



# Load clean master toxin dataset --------------------------------------

master_tox <- readRDS(
  "~/Desktop/Project/Brooks_lake_2025/data_clean/toxins/master_tox.rds"
)

# Output folder for figures --------------------------------------------

fig_dir <- "~/Desktop/Project/Brooks_lake_2025/figures/toxins"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# Define toxin columns -------------------------------------------------

toxin_cols <- c(
  "rr", "yr", "lr", "la", "dm_lr", "ly",
  "nod", "lf", "wr", "atx", "hatx", "cyl"
)

microcystin_cols <- c(
  "rr", "yr", "lr", "la", "dm_lr", "ly", "lf", "wr"
)


# Clean plotting dataset -----------------------------------------------

exclude_sites <- c(
  "SPATT_BLANK",
  "GRAB_BLANK",
  "BKS_BL_BON",
  "BKS_BL_CREEK",
  "BYS_BR_SH_01",
  "BYS_PC_SH_02",
  "BYS_FR_SH_03",
  "BYS_CR_SH_04",
  "BYS_BY_BU_SS",
  "BYS_BY_BU_DD",
  "BKS_RN_INFLOW",
  "BKS_UB_NCOVE",
  "Boysen Dup" 

)

tox_plot <- master_tox %>%
  mutate(
    method = str_to_lower(method),
    sample_type = str_to_lower(sample_type),
    site_type = str_to_lower(site_type),
    lake = str_to_lower(lake),
    sample_status = str_to_lower(sample_status),
    sample_date = as.Date(sample_date)
  ) %>%
  filter(
    !site_id %in% exclude_sites,
    !str_detect(sample_type, "duplicate|blank"),
    !str_detect(site_id, "BLANK"),
    sample_status != "not_analyzed"
  )


saveRDS(
  tox_plot,
  "data_clean/toxins/tox_plot.rds"
)
# Split grab and SPATT -------------------------------------------------

grab_tox <- tox_plot %>%
  filter(method == "grab")

spatt_tox <- tox_plot %>%
  filter(method == "spatt")

saveRDS(
  spatt_tox,
  "data_clean/toxins/spatt_tox.rds"
)

saveRDS(
  grab_tox,
  "data_clean/toxins/grab_tox.rds"
)


# Figure 1: Congener detection by lake ---------------------------------

congener_detection_lake <- grab_tox %>%
  select(lake, all_of(toxin_cols)) %>%
  pivot_longer(
    cols = all_of(toxin_cols),
    names_to = "toxin",
    values_to = "concentration"
  ) %>%
  group_by(lake, toxin) %>%
  summarize(
    n_samples = sum(!is.na(concentration)),
    detections = sum(concentration > 0, na.rm = TRUE),
    detection_frequency = detections / n_samples,
    .groups = "drop"
  )

fig1 <- ggplot(
  congener_detection_lake,
  aes(
    x = toxin,
    y = lake,
    fill = detection_frequency
  )
) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    labels = percent_format(),
    na.value = "grey90"
  ) +
  labs(
    title = "Grab Sample Toxin Detection Frequency by Lake",
    x = "Toxin / congener",
    y = "Lake",
    fill = "Detection frequency"
  ) +
  theme_minimal()

fig1

ggsave(
  file.path(fig_dir, "fig1_congener_detection_by_lake.png"),
  fig1,
  width = 9,
  height = 5,
  dpi = 300
)


# Figure 2: Surface buoy dynamics over time ----------------------------

surface_buoy_grab <- grab_tox %>%
  filter(site_type == "buoy_surface")

fig2 <- ggplot(
  surface_buoy_grab,
  aes(
    x = sample_date,
    y = total_mc
  )
) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_line(aes(group = site_id), alpha = 0.7) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "Surface Buoy Grab Sample Total Microcystins Through Time",
    x = "Sample date",
    y = "Total microcystins, log1p scale"
  ) +
  theme_minimal()

fig2

ggsave(
  file.path(fig_dir, "fig2_surface_buoy_total_mc_time_series.png"),
  fig2,
  width = 10,
  height = 6,
  dpi = 300
)


# Figure 3: Surface vs depth buoy dynamics over time -------------------

buoy_grab <- grab_tox %>%
  filter(site_type %in% c("buoy_surface", "buoy_depth"))

fig3 <- ggplot(
  buoy_grab,
  aes(
    x = sample_date,
    y = total_mc,
    color = site_type,
    group = site_type
  )
) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "Grab Sample Surface vs Depth Buoy Total Microcystins",
    x = "Sample date",
    y = "Total microcystins, log1p scale",
    color = "Site type"
  ) +
  theme_minimal()

fig3

ggsave(
  file.path(fig_dir, "fig3_surface_vs_depth_buoy_total_mc.png"),
  fig3,
  width = 10,
  height = 6,
  dpi = 300
)


# Figure 4: Shore vs buoy by lake over time ----------------------------

shore_buoy_grab <- grab_tox %>%
  filter(site_type %in% c("shore", "buoy_surface", "buoy_depth"))

fig4 <- ggplot(
  shore_buoy_grab,
  aes(
    x = sample_date,
    y = total_mc,
    color = site_type
  )
) +
  geom_point(size = 2.2, alpha = 0.75) +
  geom_line(
    aes(group = interaction(site_id, site_type)),
    alpha = 0.45
  ) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "Grab Sample Shore, Surface Buoy, and Depth Buoy Total Microcystins",
    x = "Sample date",
    y = "Total microcystins, log1p scale",
    color = "Site type"
  ) +
  theme_minimal()

fig4

ggsave(
  file.path(fig_dir, "fig4_shore_surface_depth_total_mc_time_series.png"),
  fig4,
  width = 11,
  height = 7,
  dpi = 300
)

## specifically for BROOKS ONLY
# Brooks-only grab data ------------------------------------------------

brooks_grab <- grab_tox %>%
  filter(lake == "brooks")

# Summarize Brooks shore sites by date
brooks_shore_mean <- brooks_grab %>%
  filter(site_type == "shore") %>%
  group_by(sample_date) %>%
  summarize(
    total_mc = mean(total_mc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(site_group = "shore_mean")

# Keep buoy surface and buoy depth as observed values
brooks_buoy <- brooks_grab %>%
  filter(site_type %in% c("buoy_surface", "buoy_depth")) %>%
  select(sample_date, total_mc, site_group = site_type)

# Combine shore mean with buoy data
brooks_shore_buoy_mean <- bind_rows(
  brooks_shore_mean,
  brooks_buoy
)

# Plot
fig_brooks_mean <- ggplot(
  brooks_shore_buoy_mean,
  aes(
    x = sample_date,
    y = total_mc,
    color = site_group,
    group = site_group
  )
) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_line(alpha = 0.8) +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "Brooks Lake Grab Samples: Mean Shore vs Buoy Sites",
    x = "Sample date",
    y = "Total microcystins, log1p scale",
    color = "Site group"
  ) +
  theme_minimal()

fig_brooks_mean





# Add shoreline region metadata
brooks_site_metadata <- tribble(
  ~site_id, ~shore_region,
  "BKS_BL_SH_01", "south",
  "BKS_BL_SH_02", "south",
  "BKS_BL_SH_03", "south",
  "BKS_BL_SH_04", "west",
  "BKS_BL_SH_05", "north",
  "BKS_BL_SH_06", "north",
  "BKS_BL_SH_07", "east"
)

# Brooks shoreline data only
brooks_shore_position <- brooks_grab %>%
  filter(site_type == "shore") %>%
  left_join(
    brooks_site_metadata,
    by = "site_id"
  )

# Plot shoreline sites by region
fig_brooks_position <- ggplot(
  brooks_shore_position,
  aes(
    x = sample_date,
    y = total_mc,
    color = shore_region,
    group = site_id
  )
) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_line(alpha = 0.65) +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "Brooks Lake Shoreline Grab Samples by Shore Region",
    x = "Sample date",
    y = "Total microcystins, log1p scale",
    color = "Shore region"
  ) +
  theme_minimal()

fig_brooks_position

ggsave(
  file.path(fig_dir, "brooks_mean_shore_vs_buoy.png"),
  fig_brooks_mean,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "brooks_shoreline_by_region.png"),
  fig_brooks_position,
  width = 9,
  height = 5,
  dpi = 300
)


# Figure 5: SPATT vs grab detection frequency --------------------------

method_detection <- tox_plot %>%
  filter(method %in% c("grab", "spatt")) %>%
  select(method, all_of(toxin_cols)) %>%
  pivot_longer(
    cols = all_of(toxin_cols),
    names_to = "toxin",
    values_to = "concentration"
  ) %>%
  group_by(method, toxin) %>%
  summarize(
    n_samples = sum(!is.na(concentration)),
    detections = sum(concentration > 0, na.rm = TRUE),
    detection_frequency = detections / n_samples,
    .groups = "drop"
  )

fig5 <- ggplot(
  method_detection,
  aes(
    x = toxin,
    y = detection_frequency,
    fill = method
  )
) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Toxin Detection Frequency by Method",
    x = "Toxin / congener",
    y = "Detection frequency",
    fill = "Method"
  ) +
  theme_minimal()

fig5

ggsave(
  file.path(fig_dir, "fig5_detection_frequency_grab_vs_spatt.png"),
  fig5,
  width = 9,
  height = 5,
  dpi = 300
)


# Figure 6: SPATT detection through time -------------------------------

spatt_detection_time <- spatt_tox %>%
  mutate(
    total_toxin_detected = if_any(
      all_of(toxin_cols),
      ~ !is.na(.x) & .x > 0
    )
  )

fig6 <- ggplot(
  spatt_detection_time,
  aes(
    x = sample_date,
    y = site_id,
    fill = total_toxin_detected
  )
) +
  geom_tile(color = "white") +
  facet_wrap(~ lake, scales = "free_y") +
  labs(
    title = "SPATT Toxin Detection Through Time",
    x = "Sample date",
    y = "Site ID",
    fill = "Any toxin detected?"
  ) +
  theme_minimal()

fig6

ggsave(
  file.path(fig_dir, "fig6_spatt_detection_through_time.png"),
  fig6,
  width = 11,
  height = 7,
  dpi = 300
)