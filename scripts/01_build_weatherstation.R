# ======================================================================
# Brooks Lake weather station data cleaning
# Builds cleaned 10-minute and daily weather datasets
# No plotting in this script
# ======================================================================

library(tidyverse)
library(readr)
library(lubridate)

# ----------------------------------------------------------------------
# 1. Read raw weather station file
# ----------------------------------------------------------------------
# The first row is only the title: "Brooks Lake Weather Station"
# The second row contains the real column names.
# skip = 1 removes the title row.

brooks_weather_raw <- read_csv(
  "data_raw/weather/BrooksWeather_2025.csv",
  skip = 1,
  show_col_types = FALSE
)

# ----------------------------------------------------------------------
# 2. Rename columns by position
# ----------------------------------------------------------------------
# The original names have symbols like µ, °, %, ø, and sensor IDs.
# Renaming by position avoids problems with those messy names.

brooks_weather_raw <- brooks_weather_raw %>%
  rename(
    record_id = 1,
    datetime = 2,
    pressure_inhg = 3,
    par_umol_m2_s = 4,
    air_temp_f = 5,
    relative_humidity_pct = 6,
    wind_speed_mph = 7,
    gust_speed_mph = 8,
    wind_direction_deg = 9
  )

# ----------------------------------------------------------------------
# 3. Field season cutoff
# ----------------------------------------------------------------------
# Keep data only through the final useful weather observation.
# Anything after 10/14/25 03:28:31 PM is removed.

field_season_end <- mdy_hms(
  "10/14/25 03:28:31 PM",
  tz = "Etc/GMT+7"
)

# ----------------------------------------------------------------------
# 4. Convert wind degrees to compass labels
# ----------------------------------------------------------------------
# Weather stations report wind direction in degrees.
#
# Meteorological convention:
#   0° or 360° = wind from north
#   90°        = wind from east
#   180°       = wind from south
#   270°       = wind from west
#
# Important:
# Wind direction means where the wind is coming FROM.
#
# Example:
#   wind_direction_deg = 0
#   means wind is coming from the north and blowing toward the south.
#
# This function converts degrees into 8 direction categories.

direction_label <- function(deg) {
  case_when(
    is.na(deg) ~ NA_character_,
    deg >= 337.5 | deg < 22.5 ~ "N",
    deg >= 22.5  & deg < 67.5 ~ "NE",
    deg >= 67.5  & deg < 112.5 ~ "E",
    deg >= 112.5 & deg < 157.5 ~ "SE",
    deg >= 157.5 & deg < 202.5 ~ "S",
    deg >= 202.5 & deg < 247.5 ~ "SW",
    deg >= 247.5 & deg < 292.5 ~ "W",
    deg >= 292.5 & deg < 337.5 ~ "NW"
  )
}

# ----------------------------------------------------------------------
# 5. Circular mean for wind direction
# ----------------------------------------------------------------------
# Wind direction is circular data.
#
# A regular average can be wrong.
#
# Example:
#   mean(c(350, 10)) gives 180
#
# But 350° and 10° are both near north, so the mean should be near 0/360.
#
# This function:
#   1. Converts degrees to radians
#   2. Converts angles into sine and cosine components
#   3. Averages those components
#   4. Converts the result back to degrees

circular_mean_deg <- function(deg) {
  deg <- deg[!is.na(deg)]
  
  if (length(deg) == 0) {
    return(NA_real_)
  }
  
  radians <- deg * pi / 180
  
  mean_sin <- mean(sin(radians))
  mean_cos <- mean(cos(radians))
  
  mean_angle <- atan2(mean_sin, mean_cos) * 180 / pi
  
  (mean_angle + 360) %% 360
}

# ----------------------------------------------------------------------
# 6. Safe mode function
# ----------------------------------------------------------------------
# Finds the most common value in a group.
# Used to find the dominant wind direction each day.
# This version avoids errors if all values are NA.

safe_mode <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  names(which.max(table(x)))
}

# ----------------------------------------------------------------------
# 7. Clean 10-minute weather data
# ----------------------------------------------------------------------

brooks_weather_clean <- brooks_weather_raw %>%
  mutate(
    # Parse timestamp.
    # Example raw timestamp:
    #   "06/29/25 10:18:31 AM"
    #
    # mdy_hms() means:
    #   month/day/year hour:minute:second AM/PM
    #
    # The logger reports GMT-07:00.
    # "Etc/GMT+7" keeps the logger's fixed UTC-7 offset.
    datetime = mdy_hms(datetime, tz = "Etc/GMT+7"),
    
    # Extract date for daily summaries.
    date = as.Date(datetime),
    
    # --------------------------------------------------------------
    # Unit conversions
    # --------------------------------------------------------------
    
    # Air temperature:
    # Original unit = Fahrenheit
    # New unit = Celsius
    air_temp_c = (air_temp_f - 32) * 5 / 9,
    
    # Wind speed:
    # Original unit = miles per hour
    # New unit = meters per second
    # 1 mph = 0.44704 m/s
    wind_speed_ms = wind_speed_mph * 0.44704,
    
    # Gust speed:
    # Same conversion as wind speed.
    gust_speed_ms = gust_speed_mph * 0.44704,
    
    # Pressure:
    # Original unit = inches of mercury
    # New unit = kilopascals
    # 1 in Hg = 3.38639 kPa
    pressure_kpa = pressure_inhg * 3.38639,
    
    # --------------------------------------------------------------
    # Wind direction FROM
    # --------------------------------------------------------------
    # This is the standard weather-station wind direction.
    #
    # Example:
    #   wind_direction_deg = 0
    #   wind_from = "N"
    #
    # This means wind is coming FROM the north.
    wind_from = direction_label(wind_direction_deg),
    
    # --------------------------------------------------------------
    # Wind direction TO
    # --------------------------------------------------------------
    # For bloom accumulation, we also care where wind may push
    # floating cyanobacteria or surface scums.
    #
    # wind_to is the opposite of wind_from.
    #
    # Example:
    #   wind from north = 0°
    #   wind blows toward south = 180°
    #
    # Add 180 degrees to get the opposite direction.
    # %% 360 keeps the value between 0 and 360.
    wind_to_deg = (wind_direction_deg + 180) %% 360,
    
    # Convert wind_to degrees into compass labels.
    wind_to = direction_label(wind_to_deg)
  ) %>%
  filter(
    datetime <= field_season_end
  ) %>%
  select(
    record_id,
    datetime,
    date,
    air_temp_c,
    air_temp_f,
    relative_humidity_pct,
    pressure_kpa,
    pressure_inhg,
    par_umol_m2_s,
    wind_speed_ms,
    wind_speed_mph,
    gust_speed_ms,
    gust_speed_mph,
    wind_direction_deg,
    wind_from,
    wind_to_deg,
    wind_to
  )

# ----------------------------------------------------------------------
# 8. Check cleaned data
# ----------------------------------------------------------------------

# Confirm start and end timestamp.
range(brooks_weather_clean$datetime)

# Check whether data are mostly 10-minute intervals.
brooks_weather_interval_check <- brooks_weather_clean %>%
  arrange(datetime) %>%
  mutate(
    time_step_min = as.numeric(
      difftime(datetime, lag(datetime), units = "mins")
    )
  ) %>%
  count(time_step_min, sort = TRUE)

brooks_weather_interval_check

# ----------------------------------------------------------------------
# 9. Build daily weather summaries
# ----------------------------------------------------------------------
# These daily values are useful for comparing weather to:
#   - phytoplankton sampling dates
#   - toxin sampling dates
#   - nutrients
#   - buoy stability/mixing metrics
#
# For wind direction, this script keeps both:
#
#   wind_from_mean = average direction wind came from
#   wind_to_mean   = average direction wind blew toward
#
# and:
#
#   dominant_wind_from = most common wind-from category that day
#   dominant_wind_to   = most common wind-to category that day
#
# The mean direction and dominant direction are not always the same.

brooks_weather_daily <- brooks_weather_clean %>%
  group_by(date) %>%
  summarise(
    n_obs = n(),
    
    air_temp_mean_c = mean(air_temp_c, na.rm = TRUE),
    air_temp_min_c = min(air_temp_c, na.rm = TRUE),
    air_temp_max_c = max(air_temp_c, na.rm = TRUE),
    
    relative_humidity_mean_pct = mean(relative_humidity_pct, na.rm = TRUE),
    
    pressure_mean_kpa = mean(pressure_kpa, na.rm = TRUE),
    
    par_mean_umol_m2_s = mean(par_umol_m2_s, na.rm = TRUE),
    par_max_umol_m2_s = max(par_umol_m2_s, na.rm = TRUE),
    
    wind_speed_mean_ms = mean(wind_speed_ms, na.rm = TRUE),
    wind_speed_max_ms = max(wind_speed_ms, na.rm = TRUE),
    gust_speed_max_ms = max(gust_speed_ms, na.rm = TRUE),
    
    # Circular average of the direction wind came FROM.
    wind_from_mean_deg = circular_mean_deg(wind_direction_deg),
    
    # Opposite of mean wind-from direction.
    # This estimates the mean direction wind blew TOWARD.
    wind_to_mean_deg = (wind_from_mean_deg + 180) %% 360,
    
    # Most common 10-minute wind direction categories for the day.
    dominant_wind_from = safe_mode(wind_from),
    dominant_wind_to = safe_mode(wind_to),
    
    .groups = "drop"
  ) %>%
  mutate(
    wind_from_mean = direction_label(wind_from_mean_deg),
    wind_to_mean = direction_label(wind_to_mean_deg)
  )

# ----------------------------------------------------------------------
# 10. Build wind direction summaries
# ----------------------------------------------------------------------

# Overall frequency of wind direction FROM each compass category.
# This describes the main meteorological wind pattern.

wind_from_summary <- brooks_weather_clean %>%
  filter(!is.na(wind_from)) %>%
  count(wind_from) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  arrange(desc(prop))

# Overall frequency of wind direction TO each compass category.
# This may be more useful for bloom accumulation because it shows
# which shoreline direction floating material may be pushed toward.

wind_to_summary <- brooks_weather_clean %>%
  filter(!is.na(wind_to)) %>%
  count(wind_to) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  arrange(desc(prop))

# Smaller daily wind-only table for later plotting or joins.

daily_wind_summary <- brooks_weather_daily %>%
  select(
    date,
    wind_speed_mean_ms,
    wind_speed_max_ms,
    gust_speed_max_ms,
    wind_from_mean_deg,
    wind_from_mean,
    wind_to_mean_deg,
    wind_to_mean,
    dominant_wind_from,
    dominant_wind_to
  )

# ----------------------------------------------------------------------
# 11. Save cleaned data
# ----------------------------------------------------------------------

dir.create(
  "data_clean/weather",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  brooks_weather_clean,
  "data_clean/weather/brooks_weather_10min_clean_2025.rds"
)

write_csv(
  brooks_weather_clean,
  "data_clean/weather/brooks_weather_10min_clean_2025.csv"
)

saveRDS(
  brooks_weather_daily,
  "data_clean/weather/brooks_weather_daily_2025.rds"
)

write_csv(
  brooks_weather_daily,
  "data_clean/weather/brooks_weather_daily_2025.csv"
)

write_csv(
  wind_from_summary,
  "data_clean/weather/brooks_weather_wind_from_summary_2025.csv"
)

write_csv(
  wind_to_summary,
  "data_clean/weather/brooks_weather_wind_to_summary_2025.csv"
)

write_csv(
  daily_wind_summary,
  "data_clean/weather/brooks_weather_daily_wind_summary_2025.csv"
)