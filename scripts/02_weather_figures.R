# ======================================================================
# Brooks Lake weather plotting
#
# Purpose:
#   1. Read cleaned Brooks weather station data
#   2. Read and clean NOAA/SNOTEL daily data
#   3. Join Brooks daily weather with SNOTEL daily weather
#   4. Make exploratory weather figures for bloom interpretation
#
# Notes:
#   - Brooks weather station data include 10-minute and daily summaries.
#   - NOAA/SNOTEL data are daily values.
#   - NOAA "standard" temperature units are Fahrenheit.
#   - NOAA "standard" precipitation units are inches.
# ======================================================================


# ======================================================================
# 0. Libraries
# ======================================================================

library(tidyverse)
library(readr)
library(readxl)
library(lubridate)
library(scales)


# ======================================================================
# 1. File paths and output folders
# ======================================================================

fig_dir <- "figures/weather"

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================================
# 2. Read data
# ======================================================================

brooks_weather_clean <- readRDS(
  "data_clean/weather/brooks_weather_10min_clean_2025.rds"
)

brooks_weather_daily <- readRDS(
  "data_clean/weather/brooks_weather_daily_2025.rds"
)

snotel_raw <- read_excel(
  "data_raw/weather/SNOWTEL_2025.xlsx"
)


# ======================================================================
# 3. Clean NOAA/SNOTEL daily data
# ======================================================================

snotel_daily <- snotel_raw %>%
  transmute(
    date = as.Date(DATE),

    # Precipitation: NOAA standard units are inches.
    snotel_prcp_in = PRCP,
    snotel_prcp_mm = PRCP * 25.4,

    # Temperature: NOAA standard units are Fahrenheit.
    snotel_tavg_f = TAVG,
    snotel_tmax_f = TMAX,
    snotel_tmin_f = TMIN,

    snotel_tavg_c = (TAVG - 32) * 5 / 9,
    snotel_tmax_c = (TMAX - 32) * 5 / 9,
    snotel_tmin_c = (TMIN - 32) * 5 / 9,

    # Snow variables are retained, but appear to be zero for the 2025
    # June-October Brooks study period.
    snotel_snow_depth_in = SNWD,
    snotel_snow_water_equiv_in = WESD
  )


# ======================================================================
# 4. Join Brooks daily weather with NOAA/SNOTEL daily data
# ======================================================================

brooks_weather_noaa_daily <- brooks_weather_daily %>%
  mutate(date = as.Date(date)) %>%
  left_join(
    snotel_daily,
    by = "date"
  ) %>%
  mutate(
    temp_difference_c = air_temp_mean_c - snotel_tavg_c
  )

saveRDS(
  brooks_weather_noaa_daily,
  "data_clean/weather/brooks_weather_noaa_daily_2025.rds"
)


# ----------------------------------------------------------------------
# Join check
# ----------------------------------------------------------------------

join_check <- brooks_weather_noaa_daily %>%
  summarise(
    n_days = n(),
    n_days_with_snotel = sum(!is.na(snotel_tavg_c)),
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE)
  )

print(join_check)


# ======================================================================
# 5. NOAA/SNOTEL + Brooks comparison figures
# ======================================================================

# ----------------------------------------------------------------------
# 5a. Brooks - SNOTEL temperature difference
# ----------------------------------------------------------------------
# Positive values mean Brooks was warmer than Togwotee Pass SNOTEL.
# Negative values mean Brooks was cooler than Togwotee Pass SNOTEL.

p_temp_difference <- ggplot(
  brooks_weather_noaa_daily,
  aes(x = date, y = temp_difference_c)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  labs(
    x = NULL,
    y = "Brooks - SNOTEL temperature difference (°C)"
  ) +
  theme_bw()

p_temp_difference

ggsave(
  filename = file.path(fig_dir, "brooks_minus_snotel_temperature_difference.png"),
  plot = p_temp_difference,
  width = 10,
  height = 4,
  dpi = 300
)


# ----------------------------------------------------------------------
# 5b. Brooks and SNOTEL mean daily temperature overlay
# ----------------------------------------------------------------------

temp_long <- brooks_weather_noaa_daily %>%
  select(
    date,
    air_temp_mean_c,
    snotel_tavg_c
  ) %>%
  pivot_longer(
    cols = -date,
    names_to = "station",
    values_to = "temperature_c"
  ) %>%
  mutate(
    station = recode(
      station,
      air_temp_mean_c = "Brooks weather station",
      snotel_tavg_c = "Togwotee Pass SNOTEL"
    )
  )

p_temp_overlay <- ggplot(
  temp_long,
  aes(
    x = date,
    y = temperature_c,
    color = station
  )
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  labs(
    x = NULL,
    y = "Mean daily air temperature (°C)",
    color = NULL
  ) +
  theme_bw()

p_temp_overlay

ggsave(
  filename = file.path(fig_dir, "brooks_snotel_temperature_overlay.png"),
  plot = p_temp_overlay,
  width = 10,
  height = 4,
  dpi = 300
)


# ----------------------------------------------------------------------
# 5c. Brooks daily temp range with SNOTEL mean overlaid
# ----------------------------------------------------------------------
# Ribbon = Brooks daily min-max air temperature.
# Solid line = Brooks daily mean air temperature.
# Dashed line = SNOTEL daily mean air temperature.

p_brooks_snotel_temp_range <- ggplot(
  brooks_weather_noaa_daily,
  aes(x = date)
) +
  geom_ribbon(
    aes(
      ymin = air_temp_min_c,
      ymax = air_temp_max_c
    ),
    alpha = 0.25
  ) +
  geom_line(
    aes(y = air_temp_mean_c),
    linewidth = 1
  ) +
  geom_line(
    aes(y = snotel_tavg_c),
    linewidth = 1,
    linetype = "dashed",
    na.rm = TRUE
  ) +
  labs(
    x = NULL,
    y = "Daily air temperature (°C)"
  ) +
  theme_bw()

p_brooks_snotel_temp_range

ggsave(
  filename = file.path(fig_dir, "brooks_snotel_daily_temperature_range.png"),
  plot = p_brooks_snotel_temp_range,
  width = 10,
  height = 4,
  dpi = 300
)


# ----------------------------------------------------------------------
# 5d. Weather forcing panels: precip, wind, and overlaid temp
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Temperature data (both stations together)
# ----------------------------------------------------------------------

temp_panel <- brooks_weather_noaa_daily %>%
  select(
    date,
    air_temp_mean_c,
    snotel_tavg_c
  ) %>%
  pivot_longer(
    cols = c(
      air_temp_mean_c,
      snotel_tavg_c
    ),
    names_to = "station",
    values_to = "value"
  ) %>%
  mutate(
    panel = "Air temperature (°C)",
    station = recode(
      station,
      air_temp_mean_c = "Brooks",
      snotel_tavg_c = "SNOTEL"
    )
  )

# ----------------------------------------------------------------------
# Wind panel
# ----------------------------------------------------------------------
wind_panel <- brooks_weather_noaa_daily %>%
  select(
    date,
    wind_speed_mean_ms,
    gust_speed_max_ms
  ) %>%
  pivot_longer(
    cols = c(
      wind_speed_mean_ms,
      gust_speed_max_ms
    ),
    names_to = "wind_metric",
    values_to = "value"
  ) %>%
  mutate(
    panel = "Wind speed (m/s)",
    wind_metric = recode(
      wind_metric,
      wind_speed_mean_ms = "Mean wind speed",
      gust_speed_max_ms = "Maximum gust"
    )
  )

# ----------------------------------------------------------------------
# Precipitation panel
# ----------------------------------------------------------------------

precip_panel <- brooks_weather_noaa_daily %>%
  transmute(
    date,
    panel = "Precipitation (mm)",
    value
    = snotel_prcp_mm
  )

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------

p_weather_summary <- ggplot() +
  
  # precipitation
  geom_col(
    data = precip_panel,
    aes(
      x = date,
      y = value
    )
  ) +
  
  # wind
  geom_line(
    data = wind_panel,
    aes(
      x = date,
      y = value,
      linetype = wind_metric
    ),
    linewidth = 0.8
  ) +
  
  # temperatures
  geom_line(
    data = temp_panel,
    aes(
      x = date,
      y = value,
      color = station
    ),
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ panel,
    ncol = 1,
    scales = "free_y"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    color = NULL,
    linetype = NULL
  ) +
  
  theme_bw()

p_weather_summary

ggsave(
  file.path(
    fig_dir,
    "weather_summary_precip_wind_temp.png"
  ),
  p_weather_summary,
  width = 10,
  height = 8,
  dpi = 300
)



# ======================================================================
# 6. Brooks-only weather station figures
# ======================================================================

# ----------------------------------------------------------------------
# 6a. 10-minute air temperature
# ----------------------------------------------------------------------

p_air_temp_10min <- ggplot(
  brooks_weather_clean,
  aes(x = datetime, y = air_temp_c)
) +
  geom_line(linewidth = 0.3) +
  labs(
    x = NULL,
    y = "Air temperature (°C)"
  ) +
  theme_bw()

p_air_temp_10min

ggsave(
  filename = file.path(fig_dir, "weather_air_temperature_10min.png"),
  plot = p_air_temp_10min,
  width = 10,
  height = 4,
  dpi = 300
)


# ----------------------------------------------------------------------
# 6b. Daily Brooks air temperature
# ----------------------------------------------------------------------

p_air_temp_daily <- ggplot(
  brooks_weather_daily,
  aes(x = date)
) +
  geom_ribbon(
    aes(
      ymin = air_temp_min_c,
      ymax = air_temp_max_c
    ),
    alpha = 0.25
  ) +
  geom_line(
    aes(y = air_temp_mean_c),
    linewidth = 1
  ) +
  labs(
    x = NULL,
    y = "Daily air temperature (°C)"
  ) +
  theme_bw()

p_air_temp_daily

ggsave(
  filename = file.path(fig_dir, "weather_daily_air_temperature.png"),
  plot = p_air_temp_daily,
  width = 10,
  height = 4,
  dpi = 300
)



# ----------------------------------------------------------------------
# 6d. Daily mean and maximum PAR
# ----------------------------------------------------------------------

p_par_daily <- ggplot(
  brooks_weather_daily,
  aes(x = date)
) +
  geom_line(
    aes(y = par_mean_umol_m2_s),
    linewidth = 1
  ) +
  geom_point(
    aes(y = par_max_umol_m2_s),
    size = 1.5,
    alpha = 0.7
  ) +
  labs(
    x = NULL,
    y = expression(paste("PAR (", mu, "mol ", m^{-2}, " ", s^{-1}, ")"))
  ) +
  theme_bw()

p_par_daily

ggsave(
  filename = file.path(fig_dir, "weather_daily_PAR.png"),
  plot = p_par_daily,
  width = 10,
  height = 4,
  dpi = 300
)



# ----------------------------------------------------------------------
# 6f. Daily mean wind speed
# ----------------------------------------------------------------------

p_wind_speed_daily <- ggplot(
  brooks_weather_daily,
  aes(x = date, y = wind_speed_mean_ms)
) +
  geom_line(linewidth = 1) +
  labs(
    x = NULL,
    y = expression(paste("Mean wind speed (m ", s^{-1}, ")"))
  ) +
  theme_bw()

p_wind_speed_daily

ggsave(
  filename = file.path(fig_dir, "weather_daily_wind_speed.png"),
  plot = p_wind_speed_daily,
  width = 10,
  height = 4,
  dpi = 300
)


# ----------------------------------------------------------------------
# 6g. Daily wind summary: mean wind, max wind, and max gust
# ----------------------------------------------------------------------
# soldi line- wind speed mean
# dashed line- wind speed max
# point- gust spee dmax 

p_daily_wind <- ggplot(
  brooks_weather_daily,
  aes(x = date)
) +
  geom_line(
    aes(y = wind_speed_mean_ms),
    linewidth = 1
  ) +
  geom_line(
    aes(y = wind_speed_max_ms),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_point(
    aes(y = gust_speed_max_ms),
    size = 1.6,
    alpha = 0.7
  ) +
  labs(
    x = NULL,
    y = expression(paste("Wind speed (m ", s^{-1}, ")"))
  ) +
  theme_bw()

p_daily_wind

ggsave(
  filename = file.path(fig_dir, "daily_wind_mean_max_gust.png"),
  plot = p_daily_wind,
  width = 10,
  height = 4,
  dpi = 300
)


# ======================================================================
# 7. Wind direction summaries
# ======================================================================

wind_from_summary_calc <- brooks_weather_clean %>%
  filter(!is.na(wind_from)) %>%
  count(wind_from) %>%
  mutate(
    prop = n / sum(n),
    wind_from = factor(
      wind_from,
      levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
    )
  )

p_wind_from_bar <- ggplot(
  wind_from_summary_calc,
  aes(x = wind_from, y = prop)
) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Direction wind blew from",
    y = "Percent of 10-minute observations"
  ) +
  theme_bw()

p_wind_from_bar


# ----------------------------------------------------------------------
# 7c. Wind TO summary recalculated from 10-minute data
# ----------------------------------------------------------------------
# This is useful if you do not want to rely on the saved summary CSV.
# It asks: which shoreline direction was surface material pushed toward
# most often?

wind_to_summary_calc <- brooks_weather_clean %>%
  filter(!is.na(wind_to)) %>%
  count(wind_to) %>%
  mutate(
    prop = n / sum(n),
    wind_to = factor(
      wind_to,
      levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
    )
  )

p_wind_to_bar <- ggplot(
  wind_to_summary_calc,
  aes(x = wind_to, y = prop)
) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Direction wind blew toward",
    y = "Percent of 10-minute observations"
  ) +
  theme_bw()

p_wind_to_bar

ggsave(
  filename = file.path(fig_dir, "wind_toward_direction_frequency.png"),
  plot = p_wind_to_bar,
  width = 7,
  height = 5,
  dpi = 300
)


# ======================================================================
# End of script
# ======================================================================
