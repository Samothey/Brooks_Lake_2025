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
weather <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/Weather/BrooksWeather_2025.csv", skip = 1, header = T) %>% 
  rename_with(~ c("datetime", "atm_pressure_mmHg", "par_surface", "temp_air", "RH", "wind_speed", "gust_speed", "wind_dir"), .cols = 2:9) %>%  # Gives better names 
  select(-X.) %>% # Removes uneccessary index column
  mutate(datetime = mdy_hms(datetime)) %>%  # Changes datetime to date time format
  mutate(datetime = round_date(datetime, "10 minutes")) %>% # Makes sure time stamps are at EVEN 10-minute intervals
  mutate(datetime = datetime + hours(1)) # Corrects timezone offset for weather data 

# Dissolved Oxygen -------------------------------------------------------------
do_surface <- read.table("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/DO/LJ_Surface.TXT", sep=",", header=TRUE, skip = 7) %>% 
  slice(-1) %>% 
  select(Mountain.Standard.Time, Temperature, Dissolved.Oxygen, Dissolved.Oxygen.Saturation) %>% # Selects only necessary columns
  rename(datetime = Mountain.Standard.Time,
         temp_1m = Temperature,
         do_mgl_1m = Dissolved.Oxygen,
         do_sat_1m = Dissolved.Oxygen.Saturation) %>% 
  mutate(datetime = ymd_hms(datetime)) %>%  # Converts to date time
  mutate(datetime = round_date(datetime, "10 minutes")) # Makes sure time stamps are at EVEN 10-minute intervals

do_depth <- read.table("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/DO/LJ_Depth.TXT", sep=",", header=TRUE, skip = 7) %>% 
  slice(-1) %>% 
  select(Mountain.Standard.Time, Temperature, Dissolved.Oxygen, Dissolved.Oxygen.Saturation) %>% 
  rename(datetime = Mountain.Standard.Time,
         temp_14m = Temperature,
         do_mgl_14m = Dissolved.Oxygen,
         do_sat_14m = Dissolved.Oxygen.Saturation) %>% 
  mutate(datetime = ymd_hms(datetime)) %>%  # Converts to date time
  mutate(datetime = round_date(datetime, "10 minutes")) # Makes sure time stamps are at EVEN 10-minute intervals

# PAR --------------------------------------------------------------------------
par <- read.table("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/PAR/LJ_PAR.TXT", sep = ",", header = T, skip = 4) %>% 
  slice(-1) %>% 
  select(Mountain.Standard.Time, Temperature, PAR) %>% 
  rename(datetime = Mountain.Standard.Time,
         temp_9m = Temperature,
         par_9m = PAR) %>% 
  mutate(datetime = ymd_hms(datetime)) %>%  # Converts to date time
  mutate(datetime = round_date(datetime, "10 minutes")) # Makes sure time stamps are at EVEN 10-minute intervals

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

temp_2m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_2m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("2m") # runs the function while reading in the data, everything happens all at once
temp_4m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_4m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("4m")
temp_5m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_5m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("5m")
temp_6m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp//LJ_6m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("6m")
temp_7m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_7m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("7m")
temp_8m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_8m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("8m")
temp_10m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_10m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("10m")
temp_11m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_11m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("11m")
temp_12m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_12m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("12m")
temp_13m <- read.csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/RAW/Temp/LJ_13m.csv", skip = 1, header = TRUE) %>% 
  clean_pendant("13m")

# Create a list of all dataframes we need to merge
df_list <- list(do_depth, do_surface, par, temp_2m, temp_4m, temp_5m, temp_6m, temp_7m, temp_8m, temp_10m, temp_11m, temp_12m, temp_13m, weather)

lower_jade <- reduce(df_list, full_join, by = "datetime") %>%  # Joins all dataframes contained in list by "datetime" column
  mutate(datetime = force_tz(datetime, tzone = "Etc/GMT+6")) %>% # Sets timezone for datetime if it doesn't already have it
  filter(
    datetime >= as.POSIXct("2025-06-29 11:30:00", tz = "Etc/GMT+6")) %>% # Filters out anything before the weather station started recording
  filter(
    datetime <= as.POSIXct("2025-10-14 07:00:00", tz = "Etc/GMT+6")) %>% # Trims the end of the dataset so there's no "out of water time"
  mutate(atm_pressure_mbar = atm_pressure_mmHg*33.8639) %>% # Converts atmospheric pressure from Hg to mbar for 
  # Selects all necessary (and no unnecessary) columns and orders them in an understandable way
  select(datetime, 
         do_mgl_1m, 
         do_sat_1m, 
         do_mgl_14m, 
         do_sat_14m, 
         par_9m,
         temp_1m, 
         temp_2m, 
         temp_4m, 
         temp_5m, 
         temp_6m, 
         temp_7m, 
         temp_8m, 
         temp_9m, 
         temp_10m, 
         temp_11m, 
         temp_12m, 
         temp_13m, 
         temp_14m,
         atm_pressure_mbar,
         par_surface,
         temp_air,
         RH,
         wind_speed,
         gust_speed,
         wind_dir) %>%
  mutate(across(-1, ~ as.numeric(.))) %>% # Makes sure all columns are actually numeric
  mutate(atm_pressure_mbar = adjust_pressure_elevation(atm_pressure_mbar, elev_from = 2758, elev_to = 2876, T_C = temp_air)) %>% # Estimates atmospheric pressure at each sample for the different elevation
  mutate(do_sat_mgl_1m = o2.at.sat.base(temp_1m, atm_pressure_mbar, model = "garcia-benson"), # Determines O2 saturation in mg/L at 1m deep
         do_sat_mgl_14m = o2.at.sat.base(temp_14m, atm_pressure_mbar, model = "garcia-benson")) %>% # Determines O2 saturation in mg/L at 14m deep
  mutate(do_sat_1m = do_mgl_1m/do_sat_mgl_1m * 100,
         do_sat_14m = do_mgl_14m/do_sat_mgl_14m * 100) %>% # Calculates a new saturation percentage based on the saturation amount cacluated in previous lines. 
  select(-do_sat_mgl_1m, -do_sat_mgl_14m) %>% # Removes unneeded columns  
  mutate(datetime = as.character(datetime)) # Converts datetime to character premptively (so we can run it through the conver to character function)


lower_jade <- convert_to_character(lower_jade) # Fills in the missing time slots

out_dir <- "/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/cleaned/"
out_file <- file.path(out_dir, "LowerJade_buoy_2025.csv")

write.csv(lower_jade, out_file, row.names = FALSE)