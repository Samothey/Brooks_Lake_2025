# ======================================================================
# Rainbow hybrid buoy + DEQ profile temperature heatmap
#
# Purpose:
#   Create an exploratory full-season Rainbow temperature heatmap using:
#   1. Actual buoy thermistor-string data while complete profiles exist
#   2. DEQ vertical profile interpolation after buoy thermistor failure
#
# Important:
#   - Buoy data are used only for hours/days where all 7 temperature
#     depths are available.
#   - DEQ profile interpolation is used for the rest of the season.
#   - This is exploratory and should be described as a hybrid visualization.
# ======================================================================

library(tidyverse)
library(lubridate)
library(janitor)
library(readxl)
library(scico)
library(scales)

# ======================================================================
# 1. Paths
# ======================================================================

fig_dir <- "figures/integrated_analysis/rainbow"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load Rainbow buoy data
# ======================================================================

Rainbow_WQ_2025_cleaned <- read_csv(
  "data_clean/buoy/Rainbow_WQ_2025_cleaned.csv",
  show_col_types = FALSE
)

rainbow_buoy <- Rainbow_WQ_2025_cleaned %>%
  clean_names() %>%
  mutate(
    datetime = datetime_mst
  )

# ======================================================================
# 3. Build hourly buoy temperature table
# ======================================================================

rainbow_temp_hourly <- rainbow_buoy %>%
  select(
    datetime,
    w_temp_c,
    w_temp1m,
    w_temp2m,
    w_temp3m,
    w_temp4m,
    w_temp5_5m,
    w_temp7m
  ) %>%
  pivot_longer(
    cols = starts_with("w_temp"),
    names_to = "depth_raw",
    values_to = "temp_c"
  ) %>%
  mutate(
    depth = case_when(
      depth_raw == "w_temp_c" ~ 0.75,
      depth_raw == "w_temp1m" ~ 1,
      depth_raw == "w_temp2m" ~ 2,
      depth_raw == "w_temp3m" ~ 3,
      depth_raw == "w_temp4m" ~ 4,
      depth_raw == "w_temp5_5m" ~ 5.5,
      depth_raw == "w_temp7m" ~ 7,
      TRUE ~ NA_real_
    ),
    datetime_hour = floor_date(datetime, unit = "hour")
  ) %>%
  filter(
    !is.na(datetime_hour),
    !is.na(depth),
    !is.na(temp_c)
  ) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 4. Keep only complete buoy hours
# ======================================================================

complete_buoy_hours <- rainbow_temp_hourly %>%
  count(datetime_hour, name = "n_depths") %>%
  filter(n_depths == 7)

rainbow_temp_hourly_complete <- rainbow_temp_hourly %>%
  semi_join(
    complete_buoy_hours,
    by = "datetime_hour"
  )

# ======================================================================
# 5. Convert complete buoy hours to daily mean profiles
# ======================================================================

rainbow_buoy_daily <- rainbow_temp_hourly_complete %>%
  mutate(
    date = as.Date(datetime_hour)
  ) %>%
  group_by(date, depth) %>%
  summarise(
    temp = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 6. Load Rainbow DEQ profile data
# ======================================================================

lake_profile_compiled <- read_excel(
  "data_clean/deq/lake_profile_compiled.xlsx"
)

rainbow_profiles <- lake_profile_compiled %>%
  mutate(
    date = as.Date(date),
    lake = str_to_lower(lake),
    depth = as.numeric(depth),
    temp = as.numeric(temp)
  ) %>%
  filter(
    lake %in% c("rainbow", "rainbow lake"),
    !is.na(depth),
    !is.na(temp)
  ) %>%
  arrange(date, depth)

# ======================================================================
# 7. Interpolation grids
# ======================================================================

depth_grid <- seq(0, 8, by = 0.2)

# Use the full seasonal range covered by DEQ profiles
date_grid <- seq(
  min(rainbow_profiles$date, na.rm = TRUE),
  max(rainbow_profiles$date, na.rm = TRUE),
  by = "day"
)

# ======================================================================
# 8. Interpolate buoy daily profiles vertically
# ======================================================================

rainbow_buoy_daily_interp <- rainbow_buoy_daily %>%
  group_by(date) %>%
  nest() %>%
  mutate(
    data = map(data, ~ arrange(.x, depth)),
    raster = map(data, function(df) {
      tibble(
        depth = depth_grid,
        temp = approx(
          x = df$depth,
          y = df$temp,
          xout = depth_grid,
          rule = 2
        )$y,
        source = "Buoy"
      )
    })
  ) %>%
  select(date, raster) %>%
  unnest(raster)

# ======================================================================
# 9. Interpolate DEQ profiles vertically at profile dates
# ======================================================================

rainbow_profile_depth_interp <- rainbow_profiles %>%
  group_by(date) %>%
  nest() %>%
  mutate(
    data = map(data, ~ arrange(.x, depth)),
    raster = map(data, function(df) {
      tibble(
        depth = depth_grid,
        temp = approx(
          x = df$depth,
          y = df$temp,
          xout = depth_grid,
          rule = 2
        )$y
      )
    })
  ) %>%
  select(date, raster) %>%
  unnest(raster)

# ======================================================================
# 10. Interpolate DEQ profiles through time
# ======================================================================

rainbow_profile_interp_full <- rainbow_profile_depth_interp %>%
  group_by(depth) %>%
  nest() %>%
  mutate(
    raster = map(data, function(df) {
      tibble(
        date = date_grid,
        temp = approx(
          x = as.numeric(df$date),
          y = df$temp,
          xout = as.numeric(date_grid),
          rule = 2
        )$y,
        source = "Profile interpolation"
      )
    })
  ) %>%
  select(depth, raster) %>%
  unnest(raster)

# ======================================================================
# 11. Build hybrid dataset
# ======================================================================
# Use buoy data through the final day with complete buoy profiles.
# Use profile interpolation after that date.

last_buoy_date <- max(rainbow_buoy_daily_interp$date, na.rm = TRUE)

rainbow_hybrid_temp <- bind_rows(
  rainbow_buoy_daily_interp %>%
    filter(date <= last_buoy_date),
  
  rainbow_profile_interp_full %>%
    filter(date > last_buoy_date)
)

# ======================================================================
# 12. QA checks
# ======================================================================

last_buoy_date

rainbow_hybrid_temp %>%
  count(source)

rainbow_hybrid_temp %>%
  summarise(
    first_date = min(date),
    last_date = max(date),
    min_temp = min(temp, na.rm = TRUE),
    max_temp = max(temp, na.rm = TRUE)
  )

# ======================================================================
# 13. Plot hybrid heatmap
# ======================================================================

p_rainbow_hybrid_temp <- ggplot(
  rainbow_hybrid_temp,
  aes(
    x = date,
    y = depth,
    fill = temp
  )
) +
  geom_tile(
    width = 1,
    height = 0.2
  ) +
  geom_vline(
    xintercept = last_buoy_date,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_y_reverse(
    breaks = seq(0, 8, by = 1)
  ) +
  scale_fill_scico(
    palette = "lajolla",
    limits = c(4, 19),
    oob = squish
  ) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %d"
  ) +
  labs(
    x = NULL,
    y = "Depth (m)",
    fill = "Temp (°C)",
    title = "Rainbow temperature: buoy record and profile interpolation"
  ) +
  annotate(
    "text",
    x = min(rainbow_hybrid_temp$date) + 7,
    y = 0.3,
    label = "Buoy",
    size = 3
  ) +
  annotate(
    "text",
    x = last_buoy_date + 35,
    y = 0.3,
    label = "Profile interpolation",
    size = 3
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_rainbow_hybrid_temp

ggsave(
  file.path(fig_dir, "rainbow_hybrid_buoy_profile_temperature_heatmap.png"),
  p_rainbow_hybrid_temp,
  width = 10,
  height = 5,
  dpi = 300
)