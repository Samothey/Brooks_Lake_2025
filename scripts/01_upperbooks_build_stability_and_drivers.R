# ======================================================================
# 01_upperbrooks_build_daily_drivers.R
#
# Purpose:
#   Build daily Upper Brooks Lake physical driver dataset:
#   - daily shared weather drivers
#   - Schmidt stability hourly and daily
#   - relative stratification state
#   - deep dissolved oxygen
#   - daily drivers table
# ======================================================================

library(tidyverse)
library(lubridate)
library(zoo)
library(rLakeAnalyzer)
library(janitor)

# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis/upperbrooks"
data_out_dir <- "data_clean/analysis/upperbrooks"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load cleaned datasets
# ======================================================================

brooks_weather_noaa_daily_2025 <- readRDS(
  "data_clean/weather/brooks_weather_noaa_daily_2025.rds"
)

upper_brooks_2025_clean <- readRDS(
  "data_clean/buoy/upper_brooks_2025_clean.rds"
)

# ======================================================================
# 3. Prepare Upper Brooks buoy data
# ======================================================================

upper_brooks <- upper_brooks_2025_clean %>%
  clean_names() %>%
  mutate(
    datetime = as.POSIXct(datetime_mst, tz = "America/Denver")
  )

# ======================================================================
# 4. Temperature profiles to hourly means
# ======================================================================

wtemp_ub <- upper_brooks %>%
  select(datetime, temp_1m, temp_2m, temp_3m, temp_4m) %>%
  pivot_longer(
    cols = starts_with("temp_"),
    names_to = "depth",
    values_to = "temp_c"
  ) %>%
  mutate(
    depth = str_remove(depth, "temp_"),
    depth = str_remove(depth, "m"),
    depth = as.numeric(depth)
  ) %>%
  filter(!is.na(temp_c))

wtemp_hourly_ub <- wtemp_ub %>%
  mutate(datetime_hour = round_date(datetime, unit = "hour")) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 5. Interpolate temperature profiles for heatmap
# ======================================================================

interp_hourly_ub <- wtemp_hourly_ub %>%
  group_by(datetime_hour) %>%
  nest() %>%
  mutate(n_depths = map_dbl(data, nrow)) %>%
  filter(n_depths > 1) %>%
  mutate(data = map(data, ~ arrange(.x, depth))) %>%
  mutate(
    interp_fun = map(data, ~ approxfun(.x$depth, .x$temp_c)),
    raster = map2(data, interp_fun, function(df, func) {
      tibble(
        depth = seq(min(df$depth), max(df$depth), by = 0.2),
        temp_c = func(depth)
      )
    })
  ) %>%
  select(datetime_hour, raster) %>%
  unnest(raster) %>%
  filter(!is.na(temp_c))

# ======================================================================
# 6. Calculate hourly Schmidt stability
# ======================================================================

upper_brooks_temp_wide <- wtemp_hourly_ub %>%
  mutate(
    depth_name = paste0("wtr_", depth)
  ) %>%
  select(
    datetime = datetime_hour,
    depth_name,
    temp_c
  ) %>%
  pivot_wider(
    names_from = depth_name,
    values_from = temp_c
  )

temp_cols_ub <- c(
  "wtr_1",
  "wtr_2",
  "wtr_3",
  "wtr_4"
)

temp_depths_ub <- c(1, 2, 3, 4)

upper_brooks_temp_wide_full <- upper_brooks_temp_wide %>%
  select(datetime, all_of(temp_cols_ub)) %>%
  filter(if_all(all_of(temp_cols_ub), ~ !is.na(.x)))

# Upper Brooks bathymetry
# Surface area = 24 acres
# Max depth = 5 m

acre_to_m2 <- 4046.856

bthD_ub <- c(0, 1, 2, 3, 4, 5)

bthA_ub <- c(
  24,
  22,
  18,
  12,
  6,
  0.5
) * acre_to_m2

schmidt_hourly_ub <- upper_brooks_temp_wide_full %>%
  rowwise() %>%
  mutate(
    schmidt_stability = schmidt.stability(
      wtr = c_across(all_of(temp_cols_ub)),
      depths = temp_depths_ub,
      bthD = bthD_ub,
      bthA = bthA_ub
    )
  ) %>%
  ungroup() %>%
  select(
    datetime_hour = datetime,
    schmidt_stability
  ) %>%
  arrange(datetime_hour) %>%
  mutate(
    stability_24hr = as.numeric(
      rollmean(
        schmidt_stability,
        k = 24,
        fill = NA,
        align = "center"
      )
    )
  )

# ======================================================================
# 7. Relative stability states
# ======================================================================
# Upper Brooks is shallow, so we classify stability relative to its own
# seasonal distribution rather than using the Brooks / Lower Jade absolute
# thresholds.

stability_q_ub <- quantile(
  schmidt_hourly_ub$stability_24hr,
  probs = c(0.33, 0.66),
  na.rm = TRUE
)

schmidt_hourly_ub <- schmidt_hourly_ub %>%
  mutate(
    strat_state_hourly = case_when(
      stability_24hr <= stability_q_ub[[1]] ~ "Low stability",
      stability_24hr <= stability_q_ub[[2]] ~ "Moderate stability",
      stability_24hr > stability_q_ub[[2]] ~ "High stability",
      TRUE ~ NA_character_
    ),
    strat_state_hourly = factor(
      strat_state_hourly,
      levels = c("Low stability", "Moderate stability", "High stability")
    )
  )

# ======================================================================
# 8. Daily Schmidt stability
# ======================================================================

schmidt_daily_ub <- schmidt_hourly_ub %>%
  mutate(date = as.Date(datetime_hour)) %>%
  filter(!is.na(stability_24hr)) %>%
  group_by(date) %>%
  summarise(
    stability_daily = mean(stability_24hr, na.rm = TRUE),
    stability_max = max(stability_24hr, na.rm = TRUE),
    stability_min = min(stability_24hr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    strat_state = case_when(
      stability_daily <= stability_q_ub[[1]] ~ "Low stability",
      stability_daily <= stability_q_ub[[2]] ~ "Moderate stability",
      stability_daily > stability_q_ub[[2]] ~ "High stability",
      TRUE ~ NA_character_
    ),
    strat_state = factor(
      strat_state,
      levels = c("Low stability", "Moderate stability", "High stability")
    )
  )

# ======================================================================
# 9. Daily dissolved oxygen
# ======================================================================
# Update these column names if your Upper Brooks DO columns differ.

do_ub <- upper_brooks %>%
  select(datetime, do_mgl_1m, do_mgl_4m) %>%
  pivot_longer(
    cols = c(do_mgl_1m, do_mgl_4m),
    names_to = "depth",
    values_to = "do_mgl"
  ) %>%
  mutate(
    depth = case_when(
      depth == "do_mgl_1m" ~ "Surface",
      depth == "do_mgl_4m" ~ "4 m",
      TRUE ~ depth
    )
  ) %>%
  filter(!is.na(do_mgl))

do_hourly_ub <- do_ub %>%
  mutate(datetime_hour = round_date(datetime, unit = "hour")) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

deep_do_daily_ub <- do_hourly_ub %>%
  filter(depth == "4 m") %>%
  mutate(date = as.Date(datetime_hour)) %>%
  group_by(date) %>%
  summarise(
    deep_do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 10. Build Upper Brooks daily drivers table
# ======================================================================

upperbrooks_daily_drivers <- brooks_weather_noaa_daily_2025 %>%
  mutate(date = as.Date(date)) %>%
  left_join(
    schmidt_daily_ub,
    by = "date"
  ) %>%
  left_join(
    deep_do_daily_ub,
    by = "date"
  )

# ======================================================================
# 11. Save outputs
# ======================================================================

saveRDS(
  schmidt_hourly_ub,
  file.path(data_out_dir, "schmidt_hourly_upperbrooks_2025.rds")
)

saveRDS(
  schmidt_daily_ub,
  file.path(data_out_dir, "schmidt_daily_upperbrooks_2025.rds")
)

saveRDS(
  deep_do_daily_ub,
  file.path(data_out_dir, "deep_do_daily_upperbrooks_2025.rds")
)

saveRDS(
  upperbrooks_daily_drivers,
  file.path(data_out_dir, "upperbrooks_daily_drivers_2025.rds")
)

# ======================================================================
# 12. Quick checks
# ======================================================================

glimpse(upperbrooks_daily_drivers)

upperbrooks_daily_drivers %>%
  summarise(
    n_days = n(),
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE),
    n_stability_days = sum(!is.na(stability_daily)),
    n_deep_do_days = sum(!is.na(deep_do_mgl))
  )