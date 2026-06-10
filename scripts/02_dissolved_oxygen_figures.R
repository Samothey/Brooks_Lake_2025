# ======================================================================
# 03_all_lakes_profile_physical_story.R
#
# Purpose:
#   Use DEQ vertical profile data to compare physical structure across lakes:
#   A. Temperature profile snapshots
#   B. DO profile snapshots
#   C. Low-oxygen layer thickness
#   D. Surface-bottom temperature difference
#
# Data:
#   lake_profile_compiled.xlsx
# ======================================================================

library(tidyverse)
library(lubridate)
library(readxl)
library(scico)
library(scales)

# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis/all_lakes_all_lakes_profiles"
data_out_dir <- "data_clean/analysis/all_lakes/all_lakes_profiles"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load and clean profile data
# ======================================================================

lake_profile_compiled <- read_excel(
  "data_clean/deq/lake_profile_compiled.xlsx"
)

profiles <- lake_profile_compiled %>%
  mutate(
    lake = str_to_title(lake),
    date = as.Date(date),
    depth = as.numeric(depth),
    temp = as.numeric(temp),
    do_mg_l = as.numeric(do_mg_l)
  ) %>%
  filter(
    !is.na(lake),
    !is.na(date),
    !is.na(depth)
  ) %>%
  arrange(lake, date, depth)

# Optional ordering
profiles <- profiles %>%
  mutate(
    lake = factor(
      lake,
      levels = c(
        "Brooks",
        "Lower Jade",
        "Rainbow",
        "Upper Brooks"
      )
    )
  )

# ======================================================================
# 3. Temperature profile snapshots
# ======================================================================

p_temp_profiles <- ggplot(
  profiles %>% filter(!is.na(temp)),
  aes(
    x = date,
    y = depth,
    fill = temp
  )
) +
  geom_tile(
    width = 2,
    height = 0.8
  ) +
  scale_y_reverse() +
  scale_fill_scico(
    palette = "lajolla",
    limits = c(4, 20),
    oob = squish
  ) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %d"
  ) +
  facet_wrap(
    ~ lake,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "Depth (m)",
    fill = "Temp (°C)",
    title = "Temperature profile snapshots across lakes"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_temp_profiles

ggsave(
  file.path(fig_dir, "all_lakes_temperature_profile_snapshots.png"),
  p_temp_profiles,
  width = 12,
  height = 8,
  dpi = 300
)

# ======================================================================
# 4. DO profile snapshots
# ======================================================================

p_do_profiles <- ggplot(
  profiles %>% filter(!is.na(do_mg_l)),
  aes(
    x = date,
    y = depth,
    fill = do_mg_l
  )
) +
  geom_tile(
    width = 2,
    height = 0.8
  ) +
  scale_y_reverse() +
  scale_fill_scico(
    palette = "lajolla",
    direction = -1,
    limits = c(0, 16),
    oob = squish
  ) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %d"
  ) +
  facet_wrap(
    ~ lake,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "Depth (m)",
    fill = "DO (mg/L)",
    title = "Dissolved oxygen profile snapshots across lakes"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_do_profiles

ggsave(
  file.path(fig_dir, "all_lakes_do_profile_snapshots.png"),
  p_do_profiles,
  width = 12,
  height = 8,
  dpi = 300
)

# ======================================================================
# 5. Interpolate DO profiles and calculate low-oxygen thickness
# ======================================================================

profile_do_interp <- profiles %>%
  filter(!is.na(do_mg_l)) %>%
  group_by(lake, date) %>%
  nest() %>%
  mutate(
    max_depth = map_dbl(data, ~ max(.x$depth, na.rm = TRUE)),
    raster = map2(data, max_depth, function(df, max_d) {
      depth_grid <- seq(0, max_d, by = 0.1)
      
      tibble(
        depth = depth_grid,
        do_mg_l = approx(
          x = df$depth,
          y = df$do_mg_l,
          xout = depth_grid,
          rule = 2
        )$y
      )
    })
  ) %>%
  select(lake, date, raster) %>%
  unnest(raster)

do_layer_summary <- profile_do_interp %>%
  group_by(lake, date) %>%
  summarise(
    max_depth = max(depth, na.rm = TRUE),
    anoxic_thickness_m = sum(do_mg_l <= 1, na.rm = TRUE) * 0.1,
    hypoxic_thickness_m = sum(do_mg_l <= 2, na.rm = TRUE) * 0.1,
    anoxic_fraction = anoxic_thickness_m / max_depth,
    hypoxic_fraction = hypoxic_thickness_m / max_depth,
    min_do_mgl = min(do_mg_l, na.rm = TRUE),
    bottom_do_mgl = do_mg_l[which.max(depth)],
    .groups = "drop"
  )

saveRDS(
  do_layer_summary,
  file.path(data_out_dir, "all_lakes_do_layer_summary_2025.rds")
)

# ======================================================================
# 6. Low-oxygen thickness through time
# ======================================================================

do_layer_long <- do_layer_summary %>%
  select(
    lake,
    date,
    anoxic_thickness_m,
    hypoxic_thickness_m
  ) %>%
  pivot_longer(
    cols = c(anoxic_thickness_m, hypoxic_thickness_m),
    names_to = "metric",
    values_to = "thickness_m"
  ) %>%
  mutate(
    metric = recode(
      metric,
      anoxic_thickness_m = "Anoxic layer ≤ 1 mg/L",
      hypoxic_thickness_m = "Hypoxic layer ≤ 2 mg/L"
    )
  )

p_low_oxygen_thickness <- ggplot(
  do_layer_long,
  aes(
    x = date,
    y = thickness_m,
    color = metric
  )
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(
    ~ lake,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "Layer thickness (m)",
    color = NULL,
    title = "Thickness of low-oxygen water through time"
  ) +
  theme_bw()

p_low_oxygen_thickness

ggsave(
  file.path(fig_dir, "all_lakes_low_oxygen_layer_thickness.png"),
  p_low_oxygen_thickness,
  width = 12,
  height = 7,
  dpi = 300
)

# ======================================================================
# 7. Percent of water column low oxygen
# ======================================================================

do_fraction_long <- do_layer_summary %>%
  select(
    lake,
    date,
    anoxic_fraction,
    hypoxic_fraction
  ) %>%
  pivot_longer(
    cols = c(anoxic_fraction, hypoxic_fraction),
    names_to = "metric",
    values_to = "fraction"
  ) %>%
  mutate(
    metric = recode(
      metric,
      anoxic_fraction = "Anoxic fraction ≤ 1 mg/L",
      hypoxic_fraction = "Hypoxic fraction ≤ 2 mg/L"
    )
  )

p_low_oxygen_fraction <- ggplot(
  do_fraction_long,
  aes(
    x = date,
    y = fraction,
    color = metric
  )
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ lake) +
  labs(
    x = NULL,
    y = "Percent of water column",
    color = NULL,
    title = "Fraction of water column with low oxygen"
  ) +
  theme_bw()

p_low_oxygen_fraction

ggsave(
  file.path(fig_dir, "all_lakes_low_oxygen_fraction.png"),
  p_low_oxygen_fraction,
  width = 12,
  height = 7,
  dpi = 300
)

# ======================================================================
# 8. Surface-bottom temperature difference
# ======================================================================

profile_temp_summary <- profiles %>%
  filter(!is.na(temp)) %>%
  group_by(lake, date) %>%
  summarise(
    surface_depth = min(depth, na.rm = TRUE),
    bottom_depth = max(depth, na.rm = TRUE),
    surface_temp = temp[which.min(depth)],
    bottom_temp = temp[which.max(depth)],
    temp_diff = surface_temp - bottom_temp,
    .groups = "drop"
  )

p_temp_diff <- ggplot(
  profile_temp_summary,
  aes(
    x = date,
    y = temp_diff,
    color = lake
  )
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  labs(
    x = NULL,
    y = "Surface-bottom temperature difference (°C)",
    color = NULL,
    title = "Surface-bottom temperature difference across lakes"
  ) +
  theme_bw()

p_temp_diff

ggsave(
  file.path(fig_dir, "all_lakes_surface_bottom_temperature_difference.png"),
  p_temp_diff,
  width = 10,
  height = 5,
  dpi = 300
)

# ======================================================================
# 9. Save summaries
# ======================================================================

saveRDS(
  profile_temp_summary,
  file.path(data_out_dir, "all_lakes_profile_temp_summary_2025.rds")
)

saveRDS(
  profile_do_interp,
  file.path(data_out_dir, "all_lakes_profile_do_interpolated_2025.rds")
)




