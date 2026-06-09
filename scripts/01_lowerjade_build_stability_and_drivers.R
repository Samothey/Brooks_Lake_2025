# ======================================================================
# 01_lowerjade_build_daily_drivers.R
#
# Purpose:
#   Build daily Lower Jade Lake physical driver dataset:
#   - Schmidt stability hourly and daily
#   - stratification state
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

fig_dir <- "figures/integrated_analysis"
data_out_dir <- "data_clean/analysis"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load Lower Jade buoy dataset + weather stations data
# ======================================================================

LowerJade_buoy_2025 <- read_csv(
  "data_clean/buoy/LowerJade_buoy_2025.csv",
  show_col_types = FALSE
)

brooks_weather_noaa_daily_2025 <- readRDS(
  "data_clean/weather/brooks_weather_noaa_daily_2025.rds"
)
# ======================================================================
# 3. Prepare Lower Jade buoy data
# ======================================================================

lowerjade <- LowerJade_buoy_2025 %>%
  clean_names() %>%
  mutate(
    datetime = as.POSIXct(datetime, tz = "America/Denver")
  )

# ======================================================================
# 4. Temperature profiles to hourly means
# ======================================================================

wtemp_lj <- lowerjade %>%
  select(datetime, matches("^temp_[0-9]+m$")) %>%
  pivot_longer(
    cols = matches("^temp_[0-9]+m$"),
    names_to = "depth",
    values_to = "temp_c"
  ) %>%
  mutate(
    depth = str_remove(depth, "temp_"),
    depth = str_remove(depth, "m"),
    depth = as.numeric(depth)
  ) %>%
  filter(!is.na(temp_c))

wtemp_hourly_lj <- wtemp_lj %>%
  mutate(datetime_hour = round_date(datetime, unit = "hour")) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 5. Calculate hourly Schmidt stability
# ======================================================================

lowerjade_temp_wide <- wtemp_hourly_lj %>%
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

temp_cols_lj <- c(
  "wtr_1",
  "wtr_2",
  "wtr_4",
  "wtr_5",
  "wtr_6",
  "wtr_7",
  "wtr_8",
  "wtr_9",
  "wtr_10",
  "wtr_11",
  "wtr_12",
  "wtr_13",
  "wtr_14"
)

temp_depths_lj <- c(
  1, 2, 4, 5, 6, 7, 8,
  9, 10, 11, 12, 13, 14
)

lowerjade_temp_wide_full <- lowerjade_temp_wide %>%
  select(datetime, all_of(temp_cols_lj)) %>%
  filter(if_all(all_of(temp_cols_lj), ~ !is.na(.x)))

# Lower Jade bathymetry
# Surface area = 16 acres
# Max depth = 15.2 m

acre_to_m2 <- 4046.856

bthD_lj <- c(
  0, 1, 2, 4, 5, 6, 7, 8,
  9, 10, 11, 12, 13, 14, 15.2
)

bthA_lj <- c(
  16,
  15.7,
  15.2,
  13.6,
  12.5,
  11.2,
  9.9,
  8.6,
  7.4,
  6.1,
  4.8,
  3.5,
  2.2,
  1.2,
  0.1
) * acre_to_m2

schmidt_hourly_lj <- lowerjade_temp_wide_full %>%
  rowwise() %>%
  mutate(
    schmidt_stability = schmidt.stability(
      wtr = c_across(all_of(temp_cols_lj)),
      depths = temp_depths_lj,
      bthD = bthD_lj,
      bthA = bthA_lj
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
    ),
    strat_state_hourly = case_when(
      stability_24hr <= 25 ~ "Mixed/Weak",
      stability_24hr <= 50 ~ "Moderate",
      stability_24hr > 50 ~ "Strong",
      TRUE ~ NA_character_
    ),
    strat_state_hourly = factor(
      strat_state_hourly,
      levels = c("Mixed/Weak", "Moderate", "Strong")
    )
  )

# ======================================================================
# 6. Daily Schmidt stability
# ======================================================================

schmidt_daily_lj <- schmidt_hourly_lj %>%
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
      stability_daily <= 25 ~ "Mixed/Weak",
      stability_daily <= 50 ~ "Moderate",
      stability_daily > 50 ~ "Strong",
      TRUE ~ NA_character_
    ),
    strat_state = factor(
      strat_state,
      levels = c("Mixed/Weak", "Moderate", "Strong")
    )
  )

# ======================================================================
# 7. Daily dissolved oxygen
# ======================================================================

do_lj <- lowerjade %>%
  select(datetime, do_mgl_1m, do_mgl_14m) %>%
  pivot_longer(
    cols = c(do_mgl_1m, do_mgl_14m),
    names_to = "depth",
    values_to = "do_mgl"
  ) %>%
  mutate(
    depth = case_when(
      depth == "do_mgl_1m" ~ "Surface",
      depth == "do_mgl_14m" ~ "14 m",
      TRUE ~ depth
    )
  ) %>%
  filter(!is.na(do_mgl))

do_hourly_lj <- do_lj %>%
  mutate(datetime_hour = round_date(datetime, unit = "hour")) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

deep_do_daily_lj <- do_hourly_lj %>%
  filter(depth == "14 m") %>%
  mutate(date = as.Date(datetime_hour)) %>%
  group_by(date) %>%
  summarise(
    deep_do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 8. Build Lower Jade daily drivers table
# ======================================================================
# Note: weather station data are from the Brooks station / NOAA daily file
# and are used as shared watershed-scale weather drivers.

lowerjade_daily_drivers <- brooks_weather_noaa_daily_2025 %>%
  mutate(date = as.Date(date)) %>%
  left_join(
    schmidt_daily_lj %>%
      mutate(date = as.Date(date)),
    by = "date"
  ) %>%
  left_join(
    deep_do_daily_lj %>%
      mutate(date = as.Date(date)),
    by = "date"
  )

# ======================================================================
# 9. Save outputs
# ======================================================================

saveRDS(
  schmidt_hourly_lj,
  file.path(data_out_dir, "schmidt_hourly_lowerjade_2025.rds")
)

saveRDS(
  schmidt_daily_lj,
  file.path(data_out_dir, "schmidt_daily_lowerjade_2025.rds")
)

saveRDS(
  deep_do_daily_lj,
  file.path(data_out_dir, "deep_do_daily_lowerjade_2025.rds")
)

saveRDS(
  lowerjade_daily_drivers,
  file.path(data_out_dir, "lowerjade_daily_drivers_2025.rds")
)

# ======================================================================
# 10. Quick checks
# ======================================================================

glimpse(lowerjade_daily_drivers)

lowerjade_daily_drivers %>%
  summarise(
    n_days = n(),
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE),
    n_stability_days = sum(!is.na(stability_daily)),
    n_deep_do_days = sum(!is.na(deep_do_mgl))
  )