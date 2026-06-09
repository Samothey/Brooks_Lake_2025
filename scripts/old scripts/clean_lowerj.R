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