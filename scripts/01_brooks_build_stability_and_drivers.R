# ======================================================================
# 01_brooks_build_daily_drivers.R
#
# Purpose:
#   Build daily Brooks Lake physical driver dataset:
#   - daily weather
#   - Schmidt stability (hourly adnd aily)
#   - stratification state
#   - deep dissolved oxygen
#.  - build daily drivers table
# ======================================================================

library(tidyverse)
library(lubridate)
library(zoo)
library(rLakeAnalyzer)

# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis"
data_out_dir <- "data_clean/analysis"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)


# ======================================================================
# 2. Load cleaned datasets
# ======================================================================

brooks_weather_noaa_daily_2025 <- readRDS(
  "data_clean/weather/brooks_weather_daily_2025.rds"
)

deq_nutrients_clean_2025 <- readRDS(
  "data_clean/deq/deq_nutrients_clean_2025.rds"
)

brooks_wq_2025_cleaned <- read_csv(
  "data_clean/buoy/brooks_wq_2025_cleaned.csv",
  show_col_types = FALSE
)
# ======================================================================
# 3. Prepare Brooks buoy temperature profile data
# ======================================================================

brooks_buoy <- brooks_wq_2025_cleaned %>%
  select(-any_of("...1")) %>%
  rename(WTemp0.75m = WTempC) %>%
  mutate(
    datetime = DATETIME_MST,
    datetime = ymd_hms(datetime, tz = "America/Denver")
  )

wtemp <- brooks_buoy %>%
  select(datetime, matches("^WTemp")) %>%
  pivot_longer(
    cols = matches("^WTemp"),
    names_to = "depth",
    values_to = "temp_c"
  ) %>%
  mutate(
    depth = str_remove(depth, "WTemp"),
    depth = str_remove(depth, "m"),
    depth = as.numeric(depth)
  ) %>%
  filter(!is.na(temp_c))

wtemp_hourly <- wtemp %>%
  mutate(datetime_hour = round_date(datetime, unit = "hour")) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 4. Interpolate temperature profiles for heatmap
# ======================================================================

interp_hourly <- wtemp_hourly %>%
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
# 5. Calculate hourly Schmidt stability
# ======================================================================

brooks_temp_wide <- wtemp_hourly %>%
  mutate(
    depth_name = paste0("wtr_", depth)
  ) %>%
  select(datetime = datetime_hour, depth_name, temp_c) %>%
  pivot_wider(
    names_from = depth_name,
    values_from = temp_c
  )

temp_cols <- c(
  "wtr_0.75",
  "wtr_1",
  "wtr_4",
  "wtr_7",
  "wtr_9",
  "wtr_10",
  "wtr_13",
  "wtr_15"
)

temp_depths <- c(0.75, 1, 4, 7, 9, 10, 13, 15)

brooks_temp_wide_full <- brooks_temp_wide %>%
  select(datetime, all_of(temp_cols)) %>%
  filter(if_all(all_of(temp_cols), ~ !is.na(.x)))

# Estimated Brooks bathymetry
bthD <- c(0.75, 1, 4, 7, 9, 10, 13, 15)

bthA <- c(
  866030,
  820000,
  650000,
  430000,
  300000,
  230000,
  90000,
  10000
)

schmidt_hourly <- brooks_temp_wide_full %>%
  rowwise() %>%
  mutate(
    schmidt_stability = schmidt.stability(
      wtr = c_across(all_of(temp_cols)),
      depths = temp_depths,
      bthD = bthD,
      bthA = bthA
    )
  ) %>%
  ungroup() %>%
  select(datetime_hour = datetime, schmidt_stability) %>%
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

schmidt_daily <- schmidt_hourly %>%
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

do_brooks <- brooks_buoy %>%
  select(datetime, odomgL, odomgL_15m) %>%
  pivot_longer(
    cols = c(odomgL, odomgL_15m),
    names_to = "depth",
    values_to = "do_mgl"
  ) %>%
  mutate(
    depth = case_when(
      depth == "odomgL" ~ "Surface",
      depth == "odomgL_15m" ~ "15 m",
      TRUE ~ depth
    )
  ) %>%
  filter(!is.na(do_mgl))

do_hourly <- do_brooks %>%
  mutate(datetime_hour = round_date(datetime, unit = "hour")) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

deep_do_daily <- do_hourly %>%
  filter(depth == "15 m") %>%
  mutate(date = as.Date(datetime_hour)) %>%
  group_by(date) %>%
  summarise(
    deep_do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------------------------------------------------
# Join daily drivers
# ----------------------------------------------------------------------
## you may have to manually load the brooks weather noaa daily 2025 because it never generates correclty ? 

brooks_daily_drivers <- brooks_weather_noaa_daily_2025 %>%
mutate(date = as.Date(date)) %>%
left_join(schmidt_daily, by = "date") %>%
left_join(deep_do_daily, by = "date")


# ----------------------------------------------------------------------
# 9. Save outputs
# ----------------------------------------------------------------------


saveRDS(
  schmidt_hourly,
  file.path(data_out_dir, "schmidt_hourly_2025.rds")
)

saveRDS(
  schmidt_daily,
  file.path(data_out_dir, "schmidt_daily_2025.rds")
)

saveRDS(
  deep_do_daily,
  file.path(data_out_dir, "deep_do_daily_2025.rds")
)

saveRDS(
  brooks_daily_drivers,
  file.path(data_out_dir, "brooks_daily_drivers_2025.rds")
)