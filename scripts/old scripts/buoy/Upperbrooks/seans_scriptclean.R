## Data cleaning Upper Brooks 

library(tidyverse)
library(lubridate)
library(purrr)
library(dplyr)
library(LakeMetabolizer)

# Adds in 00:00:00 times that were removed by converting datetime to character (do this before writing to csv)
convert_to_character <- function(df) {
  na_test <- as.POSIXct(df$datetime, format = "%Y-%m-%d %H:%M") # Checks for NA's (will only find them at midnight)
  na_indices <- is.na(na_test) # Gets indices of NA's 
  df$datetime[na_indices] <- paste(df$datetime[na_indices], "00:00:00") # Replaces with 00:00:00
  
  return(df)
}

# Function which takes atmospheric pressure at x elevation and converts it to atmospheric pressure at y elevation (in mb)
# Px_mbar = input atmospheric pressure
# T_C = Air temp at the same time as the above pressure
adjust_pressure_elevation <- function(Px_mbar, elev_from, elev_to, T_C) {
  g <- 9.80665
  M <- 0.0289644
  R <- 8.3144598
  T <- T_C + 273.15  # convert C to K
  
  delta_h <- elev_to - elev_from
  Py_mbar <- Px_mbar * exp(-(g * M * delta_h) / (R * T))
  return(Py_mbar)
}
# Weather ---------------------------------------------------------------------
weather <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BrooksWeather_2025.csv", skip = 1, header = T) %>% 
  rename_with(~ c("datetime", "atm_pressure_mmHg", "par_surface", "temp_air", "RH", "wind_speed", "gust_speed", "wind_dir"), .cols = 2:9) %>%  # Gives better names 
  select(-X.) %>% # Removes uneccessary index column
  mutate(datetime = mdy_hms(datetime)) %>%  # Changes datetime to date time format
  mutate(datetime = round_date(datetime, "10 minutes")) %>% # Makes sure time stamps are at EVEN 10-minute intervals
  mutate(datetime = datetime + hours(1)) # Corrects timezone offset for weather data 

# Dissolved Oxygen -------------------------------------------------------------
do_surface <- read.table("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Upper Brooks/RAW/DO/UB_Surface.TXT", sep=",", header=TRUE, skip = 7) %>% 
  slice(-1) %>% 
  select(Mountain.Standard.Time, Temperature, Dissolved.Oxygen, Dissolved.Oxygen.Saturation) %>% 
  rename(datetime = Mountain.Standard.Time,
         temp_1m = Temperature,
         do_mgl_1m = Dissolved.Oxygen,
         do_sat_1m = Dissolved.Oxygen.Saturation) %>% 
  mutate(datetime = ymd_hms(datetime)) %>%  #
  mutate(datetime = round_date(datetime, "10 minutes"))

do_depth <- read.table("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Upper Brooks/RAW/DO/UB_Depth.TXT", sep=",", header=TRUE, skip = 7) %>% 
  slice(-1) %>% 
  select(Mountain.Standard.Time, Temperature, Dissolved.Oxygen, Dissolved.Oxygen.Saturation) %>% 
  rename(datetime = Mountain.Standard.Time,
         temp_4m = Temperature,
         do_mgl_4m = Dissolved.Oxygen,
         do_sat_4m = Dissolved.Oxygen.Saturation) %>% 
  mutate(datetime = ymd_hms(datetime)) %>%  #
  mutate(datetime = round_date(datetime, "10 minutes"))

# PAR --------------------------------------------------------------------------
par <- read.table("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Upper Brooks/RAW/PAR/UB_PAR.TXT", sep = ",", header = T, skip = 4) %>% 
  slice(-1) %>% 
  select(Mountain.Standard.Time, Temperature, PAR) %>% 
  rename(datetime = Mountain.Standard.Time,
         temp_4m_2 = Temperature,
         par_4m = PAR) %>% 
  mutate(datetime = ymd_hms(datetime)) %>%  #
  mutate(datetime = round_date(datetime, "10 minutes"))

#  Battery Volts (SURFACE MiniDOT)
battery_surface <- read.table(
  "/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Upper Brooks/RAW/DO/UB_Surface.TXT",
  sep = ",", header = TRUE, skip = 7
) %>%
  slice(-1) %>%
  select(Mountain.Standard.Time, Battery) %>%   # <-- correct
  rename(
    datetime = Mountain.Standard.Time,
    battery_1m = Battery
  ) %>%
  mutate(
    datetime = ymd_hms(datetime),
    datetime = round_date(datetime, "10 minutes")
  )

#  Battery Volts (DEPTH MiniDOT)
battery_depth <- read.table(
  "/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Upper Brooks/RAW/DO/UB_Depth.TXT",
  sep = ",", header = TRUE, skip = 7
) %>%
  slice(-1) %>%
  select(Mountain.Standard.Time, Battery) %>%   # <-- FIXED
  rename(
    datetime = Mountain.Standard.Time,
    battery_4m = Battery                          # <-- FIXED
  ) %>%
  mutate(
    datetime = ymd_hms(datetime),
    datetime = round_date(datetime, "10 minutes")
  )

# Pendants (Temp) -------------------------------------------------------------
# Cleans the pendant data 
clean_pendant <- function(df, depth) {
  
  df <- df %>% 
    rename_with(~ c("datetime", "temp"), .cols = 2:3) %>% # renames datetime and temp to more reasonable strings
    select(datetime, temp) %>%  # selects just those two columns
    mutate(datetime = mdy_hms(datetime)) %>% # converts to datetime format
    complete(datetime = seq(min(datetime), max(datetime), by = "5 min")) %>%  # converts from 15-minute to 5-minute intervals
    arrange(datetime) %>% # arranges in sequential time order
    mutate(temp = approx(datetime, temp, datetime)$y) %>% # Interpolates the missing temps in the new 5-minute format
    mutate(datetime = floor_date(datetime, "10 minutes")) %>%  # round down to nearest 10 min
    group_by(datetime) %>%                                     # group by new 10-min bins
    summarise(temp = mean(temp, na.rm = TRUE)) %>%       # averages to 10-minutes 
    ungroup() %>%
    mutate(temp = (temp - 32) * 5/9) %>% # Converts to Celsius
    rename(!!paste0("temp", "_", depth) := temp) # Renames "temp" column to "temp_depth"
  
  return(df)
}
temp_2m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Upper Brooks/RAW/Temp/UB_2m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("2m") %>% 
  mutate(datetime = datetime + hours(1)) # had to fix this just for this one since it was GMT-7 instead of GMT-6

temp_3m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Upper Brooks/RAW/Temp/UB_3m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("3m")


df_list <- list(do_depth, do_surface, par, temp_2m, temp_3m, battery_depth, battery_surface, weather)

upper_brooks <- reduce(df_list, full_join, by = "datetime") %>% 
  mutate(datetime = force_tz(datetime, tzone = "Etc/GMT+6")) %>%
  filter(
    datetime >= as.POSIXct("2025-06-29 11:30:00", tz = "Etc/GMT+6")) %>% 
  filter(
    datetime <= as.POSIXct("2025-10-14 07:00:00", tz = "Etc/GMT+6")) %>% # Trims the end of the dataset so there's no "out of water time"
  mutate(across(-1, ~ as.numeric(.))) %>% 
  mutate(temp_4m = rowMeans(across(c(temp_4m, temp_4m_2)))) %>% 
  mutate(atm_pressure_mbar = atm_pressure_mmHg*33.8639) %>% 
  select(datetime, 
         do_mgl_1m, 
         do_sat_1m, 
         do_mgl_4m, 
         do_sat_4m, 
         par_4m,
         battery_1m,
         battery_4m,
         temp_1m, 
         temp_2m, 
         temp_3m, 
         temp_4m,
         atm_pressure_mbar,
         par_surface,
         temp_air,
         RH,
         wind_speed,
         gust_speed,
         wind_dir) %>% 
  mutate(atm_pressure_mbar = adjust_pressure_elevation(atm_pressure_mbar, elev_from = 2758, elev_to = 2775, T_C = temp_air)) %>% 
  mutate(do_sat_mgl_1m = o2.at.sat.base(temp_1m, atm_pressure_mbar, model = "garcia-benson"),
         do_sat_mgl_4m = o2.at.sat.base(temp_4m, atm_pressure_mbar, model = "garcia-benson")) %>% 
  mutate(do_sat_1m = do_mgl_1m/do_sat_mgl_1m * 100,
         do_sat_4m = do_mgl_4m/do_sat_mgl_4m * 100) %>% 
  select(-do_sat_mgl_1m, -do_sat_mgl_4m) %>% 
  mutate(datetime = as.character(datetime))

upper_brooks <- convert_to_character(upper_brooks)

write.csv(
  upper_brooks,
  "/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Upper Brooks/UpperBrooks_2025_cleaned.csv",
  row.names = FALSE
)
