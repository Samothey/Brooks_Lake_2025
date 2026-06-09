
# =========================================================
# Script: 01_build_lower_jade_buoy.R
# Purpose:
#   Import raw Lower Jade buoy files and build:
#     1) lower_jade_raw_merged  = merged raw-resolution dataset
#     2) lower_jade_hourly      = hourly aggregated dataset
#
# Design choices:
#   - UTC is the master time zone for merging and storage.
#   - Local Mountain time is added only for interpretation/plotting.
#   - HOBO temperature files are recorded in local Mountain time
#     and exported with a GMT-offset style header, so they are
#     parsed as America/Denver and then converted to UTC.
#   - miniDOT files contain both UTC and local time columns;
#     UTC is used directly.
#
# Notes:
#   - HOBO files are 15-minute temperature records.
#   - miniDOT files are 10-minute records for temp + DO.
#   - Some sensors start and end at different times.
#   - This script does NOT apply QA/QC corrections yet.
# =========================================================

# -------------------------------
# 1) Load packages
# -------------------------------
library(dplyr)
library(lubridate)
library(stringr)
library(purrr)
library(readr)
library(xts)
library(dygraphs)
library(RColorBrewer)

# -------------------------------
# 2) Define file paths
# -------------------------------
base_dir <- "~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade"

hobo_files <- c(
  "2"  = file.path(base_dir, "LJ_2m.csv"),
  "4"  = file.path(base_dir, "LJ_4m.csv"),
  "5"  = file.path(base_dir, "LJ_5m.csv"),
  "6"  = file.path(base_dir, "LJ_6m.csv"),
  "7"  = file.path(base_dir, "LJ_7m.csv"),
  "8"  = file.path(base_dir, "LJ_8m.csv"),
  "10" = file.path(base_dir, "LJ_10m.csv"),
  "11" = file.path(base_dir, "LJ_11m.csv"),
  "12" = file.path(base_dir, "LJ_12m.csv"),
  "13" = file.path(base_dir, "LJ_13m.csv")
)

surface_file <- file.path(base_dir, "LJ_Surface.TXT")
depth_file   <- file.path(base_dir, "LJ_Depth.TXT")

# -------------------------------
# 3) Helper: clean HOBO temperature file
# -------------------------------
# HOBO exports:
#   - include 1 metadata/header line to skip
#   - store local Mountain time in 12-hour format with AM/PM
#   - temperature is exported in Fahrenheit
clean_hobo_temp <- function(file, depth_m) {
  
  raw <- read.csv(file, skip = 1, header = TRUE)
  
  datetime_col <- names(raw)[str_detect(names(raw), "Date.Time")]
  temp_col     <- names(raw)[str_detect(names(raw), "Temp")]
  
  if (length(datetime_col) != 1) {
    stop("Could not uniquely identify HOBO datetime column in: ", file)
  }
  
  if (length(temp_col) != 1) {
    stop("Could not uniquely identify HOBO temperature column in: ", file)
  }
  
  cleaned <- raw %>%
    transmute(
      # Parse HOBO clock as local Mountain time, then convert to UTC
      datetime_local = mdy_hms(.data[[datetime_col]], tz = "America/Denver"),
      datetime_utc   = with_tz(datetime_local, "UTC"),
      
      # Convert Fahrenheit to Celsius
      temp_c = (.data[[temp_col]] - 32) * 5 / 9
    ) %>%
    select(datetime_utc, temp_c) %>%
    arrange(datetime_utc)
  
  names(cleaned)[2] <- paste0("temp_", depth_m, "m")
  
  cleaned
}

# -------------------------------
# 4) Helper: clean miniDOT file
# -------------------------------
# miniDOT exports:
#   - contain 7 metadata lines before the table
#   - include a units row as the first imported data row
#   - include both UTC and local time columns
#   - use UTC as the merge time
clean_minidot <- function(file, depth_m) {
  
  raw <- read.table(
    file,
    sep = ",",
    header = TRUE,
    skip = 7
  )
  
  cleaned <- raw %>%
    slice(-1) %>%  # remove units row
    transmute(
      datetime_utc = ymd_hms(trimws(UTC_Date_._Time), tz = "UTC"),
      temp         = as.numeric(trimws(Temperature)),
      do_mgl       = as.numeric(trimws(Dissolved.Oxygen)),
      do_sat       = as.numeric(trimws(Dissolved.Oxygen.Saturation))
    ) %>%
    arrange(datetime_utc)
  
  names(cleaned) <- c(
    "datetime_utc",
    paste0("temp_", depth_m, "m"),
    paste0("do_mgl_", depth_m, "m"),
    paste0("do_sat_", depth_m, "m")
  )
  
  cleaned
}

# -------------------------------
# 5) Import and clean all HOBO temperature files
# -------------------------------
hobo_list <- imap(hobo_files, ~ clean_hobo_temp(file = .x, depth_m = .y))

# Merge all HOBO depths into one wide table
lower_jade_hobo <- reduce(hobo_list, full_join, by = "datetime_utc") %>%
  arrange(datetime_utc)

# -------------------------------
# 6) Import and clean miniDOT surface + deep files
# -------------------------------
do_surface_clean <- clean_minidot(surface_file, depth_m = 1)
do_depth_clean   <- clean_minidot(depth_file,   depth_m = 14)

# -------------------------------
# 7) Merge all buoy components into one raw dataset
# -------------------------------
lower_jade_raw_merged <- lower_jade_hobo %>%
  full_join(do_surface_clean, by = "datetime_utc") %>%
  full_join(do_depth_clean,   by = "datetime_utc") %>%
  arrange(datetime_utc) %>%
  select(
    datetime_utc,
    temp_1m, do_mgl_1m, do_sat_1m,
    temp_2m,
    temp_4m,
    temp_5m,
    temp_6m,
    temp_7m,
    temp_8m,
    temp_10m,
    temp_11m,
    temp_12m,
    temp_13m,
    temp_14m, do_mgl_14m, do_sat_14m
  ) %>%
  mutate(
    datetime_local = with_tz(datetime_utc, "America/Denver")
  ) %>%
  relocate(datetime_local, .after = datetime_utc)

# -------------------------------
# 8) Structural checks on raw merged dataset
# -------------------------------

# Check for duplicate timestamps
duplicate_times <- lower_jade_raw_merged %>%
  count(datetime_utc) %>%
  filter(n > 1)

print(duplicate_times)

# Check time gaps between consecutive rows
time_gaps <- lower_jade_raw_merged %>%
  arrange(datetime_utc) %>%
  mutate(
    dt_mins = as.numeric(difftime(datetime_utc, lag(datetime_utc), units = "mins"))
  )

print(table(time_gaps$dt_mins, useNA = "ifany"))

# Summarize missingness by column
missing_counts <- colSums(is.na(lower_jade_raw_merged))
print(missing_counts)

# -------------------------------
# 9) Aggregate to hourly dataset
# -------------------------------
# This makes the mixed 10-min and 15-min data much easier to use
# for heatmaps, time-series plots, and later analysis.
lower_jade_hourly <- lower_jade_raw_merged %>%
  mutate(datetime_hour_utc = floor_date(datetime_utc, unit = "hour")) %>%
  group_by(datetime_hour_utc) %>%
  summarise(
    across(-c(datetime_utc, datetime_local), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    across(-datetime_hour_utc, ~ ifelse(is.nan(.x), NA, .x)),
    datetime_hour_local = with_tz(datetime_hour_utc, "America/Denver")
  ) %>%
  relocate(datetime_hour_local, .after = datetime_hour_utc)

print(head(lower_jade_hourly))
print(colSums(is.na(lower_jade_hourly)))

# -------------------------------
# 10) Interactive temperature QC plot with dygraphs
# -------------------------------
# Create a temperature-only table for interactive plotting
temp_hourly <- lower_jade_hourly %>%
  select(
    datetime_hour_local,
    temp_1m,
    temp_2m,
    temp_4m,
    temp_5m,
    temp_6m,
    temp_7m,
    temp_8m,
    temp_10m,
    temp_11m,
    temp_12m,
    temp_13m,
    temp_14m
  )

# Convert to xts for dygraphs
temp_xts <- xts(
  x = temp_hourly %>% select(-datetime_hour_local),
  order.by = temp_hourly$datetime_hour_local
)

# Distinct color palette for depths
temp_cols <- c(
  "temp_1m"  = "#D73027",
  "temp_2m"  = "#FC8D59",
  "temp_4m"  = "#FEE08B",
  "temp_5m"  = "#D9EF8B",
  "temp_6m"  = "#91CF60",
  "temp_7m"  = "#66C2A5",
  "temp_8m"  = "#3288BD",
  "temp_10m" = "#5E4FA2",
  "temp_11m" = "#7B3294",
  "temp_12m" = "#C2A5CF",
  "temp_13m" = "#A6CEE3",
  "temp_14m" = "#1F78B4"
)

dygraph(temp_xts, main = "Lower Jade hourly temperatures") %>%
  dyOptions(
    drawPoints = FALSE,
    strokeWidth = 1.2,
    colors = unname(temp_cols)
  ) %>%
  dyAxis("y", label = "Temperature (°C)") %>%
  dyAxis("x", label = "Date") %>%
  dyRangeSelector() %>%
  dyHighlight(
    highlightCircleSize = 3,
    highlightSeriesBackgroundAlpha = 0.2,
    hideOnMouseOut = TRUE
  ) %>%
  dyLegend(show = "follow") %>%
  dySeries("temp_1m",  label = "1 m") %>%
  dySeries("temp_2m",  label = "2 m") %>%
  dySeries("temp_4m",  label = "4 m") %>%
  dySeries("temp_5m",  label = "5 m") %>%
  dySeries("temp_6m",  label = "6 m") %>%
  dySeries("temp_7m",  label = "7 m") %>%
  dySeries("temp_8m",  label = "8 m") %>%
  dySeries("temp_10m", label = "10 m") %>%
  dySeries("temp_11m", label = "11 m") %>%
  dySeries("temp_12m", label = "12 m") %>%
  dySeries("temp_13m", label = "13 m") %>%
  dySeries("temp_14m", label = "14 m")

# -------------------------------
# 11) Optional save outputs
# -------------------------------
saveRDS(lower_jade_raw_merged,"~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/lower_jade_raw_merged.rds"
 )
#
saveRDS(lower_jade_hourly,"~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/lower_jade_raw_hourly_merged.rds")





### code to begin cleaning 
# =========================================================
# Script: 02_clean_lower_jade_buoy.R
# Purpose:
#   Apply basic QA/QC cleaning to the hourly Lower Jade buoy
#   dataset created in 01_build_lower_jade_buoy.R.
#
# Main tasks:
#   1) Load hourly data
#   2) Flag obvious deployment/startup artifacts
#   3) Flag impossible or suspicious values
#   4) Create a cleaned version for downstream plots/analysis
#   5) Save both flagged and cleaned outputs
#
# Notes:
#   - This script should preserve transparency.
#   - Do not overwrite the original hourly dataset.
#   - Keep a column-based record of what was flagged.
# =========================================================

# -------------------------------
# 1) Load packages
# -------------------------------
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(stringr)

# -------------------------------
# 2) Load hourly dataset
# -------------------------------
# Option A: if object is already in memory, skip this section
# Option B: load from saved RDS from the build script

# lower_jade_hourly <- readRDS(
#   "~/Desktop/Project/Brooks_lake_2025/data_processed/buoy/lower_jade_hourly.rds"
# )

# -------------------------------
# 3) Inspect structure
# -------------------------------
glimpse(lower_jade_hourly)
summary(lower_jade_hourly)

# -------------------------------
# 4) Create working copy
# -------------------------------
lower_jade_flagged <- lower_jade_hourly %>%
  mutate(
    # master row-level notes field
    qc_notes = NA_character_,
    
    # deployment/startup flag for deep miniDOT
    flag_startup_14m = FALSE,
    
    # optional retrieval/end-of-record flag if needed later
    flag_end_14m = FALSE
  )

# -------------------------------
# 5) Flag known startup artifact at deep miniDOT
# -------------------------------
# Based on initial inspection, the 14 m miniDOT shows a clear
# equilibration/deployment artifact at the start of the record.
#
# Replace the datetime below once you decide on the final cutoff.
# For now, this is just a starter example.

startup_cutoff_utc <- ymd_hms("2025-06-25 00:00:00", tz = "UTC")

lower_jade_flagged <- lower_jade_flagged %>%
  mutate(
    flag_startup_14m = datetime_hour_utc < startup_cutoff_utc &
      (!is.na(temp_14m) | !is.na(do_mgl_14m) | !is.na(do_sat_14m)),
    
    qc_notes = case_when(
      flag_startup_14m ~ "14 m miniDOT startup/deployment equilibration artifact",
      TRUE ~ qc_notes
    )
  )

# -------------------------------
# 6) Flag simple range problems
# -------------------------------
# These are broad screening thresholds, not final ecological rules.
# Adjust later if needed.

lower_jade_flagged <- lower_jade_flagged %>%
  mutate(
    flag_temp_1m_range  = !is.na(temp_1m)  & (temp_1m  < 0 | temp_1m  > 30),
    flag_temp_2m_range  = !is.na(temp_2m)  & (temp_2m  < 0 | temp_2m  > 30),
    flag_temp_4m_range  = !is.na(temp_4m)  & (temp_4m  < 0 | temp_4m  > 30),
    flag_temp_5m_range  = !is.na(temp_5m)  & (temp_5m  < 0 | temp_5m  > 30),
    flag_temp_6m_range  = !is.na(temp_6m)  & (temp_6m  < 0 | temp_6m  > 30),
    flag_temp_7m_range  = !is.na(temp_7m)  & (temp_7m  < 0 | temp_7m  > 30),
    flag_temp_8m_range  = !is.na(temp_8m)  & (temp_8m  < 0 | temp_8m  > 30),
    flag_temp_10m_range = !is.na(temp_10m) & (temp_10m < 0 | temp_10m > 30),
    flag_temp_11m_range = !is.na(temp_11m) & (temp_11m < 0 | temp_11m > 30),
    flag_temp_12m_range = !is.na(temp_12m) & (temp_12m < 0 | temp_12m > 30),
    flag_temp_13m_range = !is.na(temp_13m) & (temp_13m < 0 | temp_13m > 30),
    flag_temp_14m_range = !is.na(temp_14m) & (temp_14m < 0 | temp_14m > 30),
    
    flag_do_1m_range    = !is.na(do_mgl_1m)  & (do_mgl_1m  < 0 | do_mgl_1m  > 20),
    flag_do_14m_range   = !is.na(do_mgl_14m) & (do_mgl_14m < 0 | do_mgl_14m > 20),
    
    flag_do_sat_1m_range  = !is.na(do_sat_1m)  & (do_sat_1m  < 0 | do_sat_1m  > 200),
    flag_do_sat_14m_range = !is.na(do_sat_14m) & (do_sat_14m < 0 | do_sat_14m > 200)
  )

# -------------------------------
# 7) Flag large hour-to-hour jumps
# -------------------------------
# These are initial screening thresholds. Tune later after plotting.

lower_jade_flagged <- lower_jade_flagged %>%
  arrange(datetime_hour_utc) %>%
  mutate(
    flag_temp_1m_jump  = !is.na(temp_1m)  & !is.na(lag(temp_1m))  & abs(temp_1m  - lag(temp_1m))  > 3,
    flag_temp_14m_jump = !is.na(temp_14m) & !is.na(lag(temp_14m)) & abs(temp_14m - lag(temp_14m)) > 3,
    
    flag_do_1m_jump    = !is.na(do_mgl_1m)  & !is.na(lag(do_mgl_1m))  & abs(do_mgl_1m  - lag(do_mgl_1m))  > 3,
    flag_do_14m_jump   = !is.na(do_mgl_14m) & !is.na(lag(do_mgl_14m)) & abs(do_mgl_14m - lag(do_mgl_14m)) > 3
  )

# -------------------------------
# 8) Make a row-level "any flag" field
# -------------------------------
flag_cols <- names(lower_jade_flagged)[str_detect(names(lower_jade_flagged), "^flag_")]

lower_jade_flagged <- lower_jade_flagged %>%
  mutate(
    any_flag = if_any(all_of(flag_cols), identity)
  )

# -------------------------------
# 9) Create cleaned dataset
# -------------------------------
# Strategy:
#   - Only blank out values clearly affected by the known startup artifact
#   - Keep the rest of the flags for review first
#
# This is intentionally conservative.

lower_jade_clean <- lower_jade_flagged %>%
  mutate(
    temp_14m   = if_else(flag_startup_14m, NA_real_, temp_14m),
    do_mgl_14m = if_else(flag_startup_14m, NA_real_, do_mgl_14m),
    do_sat_14m = if_else(flag_startup_14m, NA_real_, do_sat_14m)
  )

# -------------------------------
# 10) Quick summaries of flagged records
# -------------------------------
flag_summary <- lower_jade_flagged %>%
  summarise(across(all_of(flag_cols), ~ sum(.x, na.rm = TRUE)))

print(flag_summary)

lower_jade_flagged %>%
  filter(any_flag) %>%
  select(datetime_hour_utc, datetime_hour_local, qc_notes, any_flag, all_of(flag_cols)) %>%
  print(n = 50)

# -------------------------------
# 11) Quick visual checks
# -------------------------------

# A. Deep temperature and DO at start of record
lower_jade_clean %>%
  filter(datetime_hour_utc < ymd_hms("2025-06-27 00:00:00", tz = "UTC")) %>%
  ggplot(aes(x = datetime_hour_local)) +
  geom_line(aes(y = temp_14m), linewidth = 0.5) +
  labs(
    title = "Lower Jade 14 m temperature after startup cleaning",
    x = "Date",
    y = "Temperature (°C)"
  ) +
  theme_minimal()

lower_jade_clean %>%
  filter(datetime_hour_utc < ymd_hms("2025-06-27 00:00:00", tz = "UTC")) %>%
  ggplot(aes(x = datetime_hour_local)) +
  geom_line(aes(y = do_mgl_14m), linewidth = 0.5) +
  labs(
    title = "Lower Jade 14 m dissolved oxygen after startup cleaning",
    x = "Date",
    y = "DO (mg/L)"
  ) +
  theme_minimal()

# B. Surface vs deep temperature check
lower_jade_clean %>%
  select(datetime_hour_local, temp_1m, temp_14m) %>%
  pivot_longer(cols = c(temp_1m, temp_14m), names_to = "depth", values_to = "temp_c") %>%
  ggplot(aes(x = datetime_hour_local, y = temp_c, color = depth)) +
  geom_line(linewidth = 0.5) +
  labs(
    title = "Lower Jade surface vs deep temperature",
    x = "Date",
    y = "Temperature (°C)",
    color = "Depth"
  ) +
  theme_minimal()

# -------------------------------
# 12) Optional save outputs
# -------------------------------
# saveRDS(
#   lower_jade_flagged,
#   "~/Desktop/Project/Brooks_lake_2025/data_processed/buoy/lower_jade_hourly_flagged.rds"
# )
#
# saveRDS(
#   lower_jade_clean,
#   "~/Desktop/Project/Brooks_lake_2025/data_processed/buoy/lower_jade_hourly_clean.rds"
# )
#
# write.csv(
#   lower_jade_clean,
#   "~/Desktop/Project/Brooks_lake_2025/data_processed/buoy/lower_jade_hourly_clean.csv",
#   row.names = FALSE
# )