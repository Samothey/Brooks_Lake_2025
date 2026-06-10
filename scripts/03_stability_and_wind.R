library(tidyverse)
library(lubridate)
library(janitor)
library(scales)

# ======================================================================
# 1. Load data
# ======================================================================

brooks_weather_10min_clean_2025 <- readRDS(
  "data_clean/weather/brooks_weather_10min_clean_2025.rds"
)

schmidt_hourly_2025 <- readRDS(
  "data_clean/analysis/brooks/schmidt_hourly_2025.rds"
)

# ======================================================================
# 2. Fix weather timestamps
# ======================================================================
# datetime is already parsed, but it represents logger time GMT-07:00.
# force_tz() assigns the correct fixed logger timezone without changing
# the clock time. with_tz() converts it to America/Denver.

weather_10min <- brooks_weather_10min_clean_2025 %>%
  clean_names() %>%
  mutate(
    datetime_logger = force_tz(datetime, tzone = "Etc/GMT+7"),
    datetime_mountain = with_tz(datetime_logger, tzone = "America/Denver"),
    datetime_hour = floor_date(datetime_mountain, unit = "hour")
  )

# ======================================================================
# 3. Convert 10-min weather to hourly
# ======================================================================

weather_hourly <- weather_10min %>%
  group_by(datetime_hour) %>%
  summarise(
    air_temp_mean_c = mean(air_temp_c, na.rm = TRUE),
    air_temp_max_c = max(air_temp_c, na.rm = TRUE),
    wind_speed_mean_ms = mean(wind_speed_ms, na.rm = TRUE),
    wind_speed_max_ms = max(wind_speed_ms, na.rm = TRUE),
    gust_speed_max_ms = max(gust_speed_ms, na.rm = TRUE),
    relative_humidity_mean_pct = mean(relative_humidity_pct, na.rm = TRUE),
    par_mean_umol_m2_s = mean(par_umol_m2_s, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 4. Prepare Brooks hourly stability
# ======================================================================

brooks_stability_hourly <- schmidt_hourly_2025 %>%
  mutate(
    datetime_hour = with_tz(datetime_hour, tzone = "America/Denver")
  ) %>%
  select(
    datetime_hour,
    schmidt_stability,
    stability_24hr
  )

# ======================================================================
# 5. Join hourly stability and weather
# ======================================================================

brooks_stability_weather_hourly <- brooks_stability_hourly %>%
  left_join(
    weather_hourly,
    by = "datetime_hour"
  )

# ======================================================================
# 6. Check join
# ======================================================================

brooks_stability_weather_hourly %>%
  summarise(
    first_time = min(datetime_hour, na.rm = TRUE),
    last_time = max(datetime_hour, na.rm = TRUE),
    n_hours = n(),
    n_weather_matches = sum(!is.na(wind_speed_mean_ms)),
    n_stability = sum(!is.na(stability_24hr))
  )

# ======================================================================
# 7. Plot stability and gusts
# ======================================================================

p_brooks_stability_gust <- ggplot(
  brooks_stability_weather_hourly,
  aes(datetime_hour)
) +
  geom_col(
    aes(
      y = gust_speed_max_ms *
        max(stability_24hr, na.rm = TRUE) /
        max(gust_speed_max_ms, na.rm = TRUE)
    ),
    fill = "#1f78b4",
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_line(
    aes(y = stability_24hr),
    color = "black",
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  labs(
    x = NULL,
    y = "Schmidt stability; gusts scaled",
    title = "Brooks stability and wind gust events"
  ) +
  theme_bw()

p_brooks_stability_gust

# ======================================================================
# 8. Zoom July event
# ======================================================================

p_brooks_july_event <- brooks_stability_weather_hourly %>%
  filter(
    datetime_hour >= ymd_hms("2025-07-10 00:00:00", tz = "America/Denver"),
    datetime_hour <= ymd_hms("2025-07-20 23:00:00", tz = "America/Denver")
  ) %>%
  ggplot(
    aes(datetime_hour)
  ) +
  geom_col(
    aes(
      y = gust_speed_max_ms *
        max(stability_24hr, na.rm = TRUE) /
        max(gust_speed_max_ms, na.rm = TRUE)
    ),
    fill = "#1f78b4",
    alpha = 0.4,
    na.rm = TRUE
  ) +
  geom_line(
    aes(y = stability_24hr),
    color = "black",
    linewidth = 1.2,
    na.rm = TRUE
  ) +
  labs(
    x = NULL,
    y = "Schmidt stability; gusts scaled",
    title = "Brooks July stability drop and wind event"
  ) +
  theme_bw()

p_brooks_july_event


###################################################

# ======================================================================
# 03_brooks_stability_change_storm_events.R
#
# Purpose:
#   Explore whether wind/gust events coincide with drops in Brooks Lake
#   Schmidt stability.
#
# What this does:
#   1. Converts 10-min weather to hourly weather.
#   2. Joins hourly weather to hourly Schmidt stability.
#   3. Calculates 24-hour stability change.
#   4. Flags large stability-loss events.
#   5. Plots stability, gusts, and stability change.
# ======================================================================

library(tidyverse)
library(lubridate)
library(janitor)
library(scales)

# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis/brooks/storm_events"
data_out_dir <- "data_clean/analysis/brooks/storm_events"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load data
# ======================================================================

brooks_weather_10min_clean_2025 <- readRDS(
  "data_clean/weather/brooks_weather_10min_clean_2025.rds"
)

schmidt_hourly_2025 <- readRDS(
  "data_clean/analysis/brooks/schmidt_hourly_2025.rds"
)

# ======================================================================
# 3. Fix weather timestamps
# ======================================================================
# Weather logger reports fixed GMT-07:00.
# The datetime column is already parsed as dttm, so use force_tz().
# Then convert to America/Denver to match buoy/stability timestamps.

weather_10min <- brooks_weather_10min_clean_2025 %>%
  clean_names() %>%
  mutate(
    datetime_logger = force_tz(datetime, tzone = "Etc/GMT+7"),
    datetime_mountain = with_tz(datetime_logger, tzone = "America/Denver"),
    datetime_hour = floor_date(datetime_mountain, unit = "hour")
  )

# ======================================================================
# 4. Convert 10-min weather to hourly
# ======================================================================

weather_hourly <- weather_10min %>%
  group_by(datetime_hour) %>%
  summarise(
    air_temp_mean_c = mean(air_temp_c, na.rm = TRUE),
    air_temp_max_c = max(air_temp_c, na.rm = TRUE),
    wind_speed_mean_ms = mean(wind_speed_ms, na.rm = TRUE),
    wind_speed_max_ms = max(wind_speed_ms, na.rm = TRUE),
    gust_speed_max_ms = max(gust_speed_ms, na.rm = TRUE),
    relative_humidity_mean_pct = mean(relative_humidity_pct, na.rm = TRUE),
    par_mean_umol_m2_s = mean(par_umol_m2_s, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 5. Prepare Brooks hourly stability
# ======================================================================

brooks_stability_hourly <- schmidt_hourly_2025 %>%
  mutate(
    datetime_hour = with_tz(datetime_hour, tzone = "America/Denver")
  ) %>%
  select(
    datetime_hour,
    schmidt_stability,
    stability_24hr
  ) %>%
  arrange(datetime_hour)

# ======================================================================
# 6. Join hourly stability and weather
# ======================================================================

brooks_stability_weather_hourly <- brooks_stability_hourly %>%
  left_join(
    weather_hourly,
    by = "datetime_hour"
  ) %>%
  arrange(datetime_hour)

# ======================================================================
# 7. Calculate stability change
# ======================================================================
# stability_change_24h:
#   positive = lake gained stability compared to 24 hours earlier
#   negative = lake lost stability compared to 24 hours earlier
#
# pct_change_24h:
#   percent change in stability over 24 hours

brooks_stability_weather_hourly <- brooks_stability_weather_hourly %>%
  mutate(
    stability_change_24h =
      stability_24hr - lag(stability_24hr, 24),
    
    pct_change_24h =
      100 * stability_change_24h / lag(stability_24hr, 24)
  )

# ======================================================================
# 8. Flag destabilization events
# ======================================================================
# Adjust these thresholds if needed.
# This flags hours where stability dropped by at least 20% over 24 hours.

stability_loss_threshold_pct <- -20

brooks_stability_weather_hourly <- brooks_stability_weather_hourly %>%
  mutate(
    stability_event = case_when(
      pct_change_24h <= stability_loss_threshold_pct ~ "Large stability loss",
      pct_change_24h >= 20 ~ "Large stability gain",
      TRUE ~ "No major change"
    ),
    stability_event = factor(
      stability_event,
      levels = c(
        "Large stability loss",
        "No major change",
        "Large stability gain"
      )
    )
  )

# ======================================================================
# 9. Summarise biggest stability-loss events
# ======================================================================

biggest_stability_losses <- brooks_stability_weather_hourly %>%
  filter(!is.na(pct_change_24h)) %>%
  arrange(pct_change_24h) %>%
  select(
    datetime_hour,
    stability_24hr,
    stability_change_24h,
    pct_change_24h,
    wind_speed_mean_ms,
    wind_speed_max_ms,
    gust_speed_max_ms,
    air_temp_mean_c
  ) %>%
  slice_head(n = 25)

biggest_stability_losses

saveRDS(
  brooks_stability_weather_hourly,
  file.path(data_out_dir, "brooks_stability_weather_hourly_2025.rds")
)

saveRDS(
  biggest_stability_losses,
  file.path(data_out_dir, "brooks_biggest_stability_losses_2025.rds")
)

# ======================================================================
# 10. Plot absolute stability and gusts
# ======================================================================

p_brooks_stability_gust <- ggplot(
  brooks_stability_weather_hourly,
  aes(datetime_hour)
) +
  geom_col(
    aes(
      y = gust_speed_max_ms *
        max(stability_24hr, na.rm = TRUE) /
        max(gust_speed_max_ms, na.rm = TRUE)
    ),
    fill = "#1f78b4",
    alpha = 0.30,
    na.rm = TRUE
  ) +
  geom_line(
    aes(y = stability_24hr),
    color = "black",
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  labs(
    x = NULL,
    y = "Schmidt stability; gusts scaled",
    title = "Brooks stability and wind gust events"
  ) +
  theme_bw()

p_brooks_stability_gust

ggsave(
  file.path(fig_dir, "brooks_stability_and_gusts_hourly.png"),
  p_brooks_stability_gust,
  width = 12,
  height = 5,
  dpi = 300
)

# ======================================================================
# 11. Plot 24-hour stability change
# ======================================================================

p_brooks_stability_change <- ggplot(
  brooks_stability_weather_hourly,
  aes(datetime_hour, stability_change_24h)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray40"
  ) +
  geom_col(
    aes(fill = stability_change_24h < 0),
    alpha = 0.8,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "#1f78b4",
      "FALSE" = "#d95f02"
    ),
    labels = c(
      "TRUE" = "Stability loss",
      "FALSE" = "Stability gain"
    ),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = expression(Delta*" Schmidt stability over 24 h"),
    title = "Brooks 24-hour stability change"
  ) +
  theme_bw()

p_brooks_stability_change

ggsave(
  file.path(fig_dir, "brooks_24hr_stability_change.png"),
  p_brooks_stability_change,
  width = 12,
  height = 5,
  dpi = 300
)

# ======================================================================
# 12. Plot gusts with 24-hour stability change
# ======================================================================

p_brooks_gust_stability_change <- ggplot(
  brooks_stability_weather_hourly,
  aes(datetime_hour)
) +
  geom_col(
    aes(y = stability_change_24h),
    fill = "gray60",
    alpha = 0.7,
    na.rm = TRUE
  ) +
  geom_line(
    aes(
      y = gust_speed_max_ms *
        max(abs(stability_change_24h), na.rm = TRUE) /
        max(gust_speed_max_ms, na.rm = TRUE)
    ),
    color = "#e31a1c",
    linewidth = 1,
    na.rm = TRUE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = expression(Delta*" stability; gusts scaled"),
    title = "Do gust events coincide with stability loss?"
  ) +
  theme_bw()

p_brooks_gust_stability_change

ggsave(
  file.path(fig_dir, "brooks_gusts_vs_24hr_stability_change.png"),
  p_brooks_gust_stability_change,
  width = 12,
  height = 5,
  dpi = 300
)

# ======================================================================
# 13. Zoom: July turnover / stability-loss window
# ======================================================================

july_window <- brooks_stability_weather_hourly %>%
  filter(
    datetime_hour >= ymd_hms("2025-07-10 00:00:00", tz = "America/Denver"),
    datetime_hour <= ymd_hms("2025-07-20 23:00:00", tz = "America/Denver")
  )

p_brooks_july_event <- ggplot(
  july_window,
  aes(datetime_hour)
) +
  geom_col(
    aes(
      y = gust_speed_max_ms *
        max(stability_24hr, na.rm = TRUE) /
        max(gust_speed_max_ms, na.rm = TRUE)
    ),
    fill = "#1f78b4",
    alpha = 0.35,
    na.rm = TRUE
  ) +
  geom_line(
    aes(y = stability_24hr),
    color = "black",
    linewidth = 1.2,
    na.rm = TRUE
  ) +
  labs(
    x = NULL,
    y = "Schmidt stability; gusts scaled",
    title = "Brooks July stability drop and wind event"
  ) +
  theme_bw()

p_brooks_july_event

ggsave(
  file.path(fig_dir, "brooks_july_stability_drop_wind_event.png"),
  p_brooks_july_event,
  width = 10,
  height = 4.5,
  dpi = 300
)

# ======================================================================
# 14. Zoom: July stability change + gusts
# ======================================================================

p_brooks_july_change_event <- ggplot(
  july_window,
  aes(datetime_hour)
) +
  geom_col(
    aes(y = stability_change_24h),
    fill = "gray60",
    alpha = 0.75,
    na.rm = TRUE
  ) +
  geom_line(
    aes(
      y = gust_speed_max_ms *
        max(abs(stability_change_24h), na.rm = TRUE) /
        max(gust_speed_max_ms, na.rm = TRUE)
    ),
    color = "#e31a1c",
    linewidth = 1,
    na.rm = TRUE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = NULL,
    y = expression(Delta*" stability; gusts scaled"),
    title = "Brooks July 24-hour stability change and wind gusts"
  ) +
  theme_bw()

p_brooks_july_change_event

ggsave(
  file.path(fig_dir, "brooks_july_stability_change_and_gusts.png"),
  p_brooks_july_change_event,
  width = 10,
  height = 4.5,
  dpi = 300
)

# ======================================================================
# 15. Scatter: gust speed vs future stability change
# ======================================================================
# Tests whether strong gusts are associated with stability loss
# over the following 24 hours.

brooks_stability_weather_hourly <- brooks_stability_weather_hourly %>%
  arrange(datetime_hour) %>%
  mutate(
    future_stability_change_24h =
      lead(stability_24hr, 24) - stability_24hr,
    future_pct_change_24h =
      100 * future_stability_change_24h / stability_24hr
  )

p_brooks_gust_future_change <- ggplot(
  brooks_stability_weather_hourly,
  aes(
    x = gust_speed_max_ms,
    y = future_stability_change_24h
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray40"
  ) +
  geom_point(
    alpha = 0.45,
    size = 1.7,
    na.rm = TRUE
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    color = "black",
    na.rm = TRUE
  ) +
  labs(
    x = "Hourly max gust speed (m/s)",
    y = expression("Future 24-h " * Delta * " stability"),
    title = "Wind gusts vs subsequent stability change"
  ) +
  theme_bw()

p_brooks_gust_future_change

ggsave(
  file.path(fig_dir, "brooks_gusts_vs_future_24hr_stability_change.png"),
  p_brooks_gust_future_change,
  width = 7,
  height = 5,
  dpi = 300
)


ggplot(
  brooks_stability_weather_hourly,
  aes(datetime_hour)
) +
  geom_line(
    aes(y = stability_change_24h),
    color = "black"
  ) +
  geom_line(
    aes(
      y = scales::rescale(
        gust_speed_max_ms,
        to = range(stability_change_24h, na.rm = TRUE)
      )
    ),
    color = "red"
  ) +
  geom_hline(yintercept = 0, linetype = 2)