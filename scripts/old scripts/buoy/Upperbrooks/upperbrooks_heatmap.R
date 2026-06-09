
# ------------------------------------------------------------
# Upper Brooks buoy: prepare hourly temperature matrix + plot
# a temperature heatmap using rLakeAnalyzer
# ------------------------------------------------------------

library(tidyverse)     # dplyr + ggplot + tidyr, etc.
library(lubridate)     # parse_date_time(), floor_date()
library(rLakeAnalyzer) # load.ts(), wtr.heat.map()

# ---- 1) Load the cleaned Upper Brooks buoy dataset ----
upper_brooks <- readRDS(
  "/Users/samanthapena/Desktop/GISProject/upper_brooks_2025_clean.rds"
)

# ---- 2) Convert 10-minute data to hourly means in "var_depth" format ----
# rLakeAnalyzer likes a "wide" format where:
#   - first column is datetime
#   - remaining columns are named var_1, var_2, ... where the number = depth (m)
ub_heat <- upper_brooks %>%
  # Parse the datetime strings in datetime_MST into a true POSIXct datetime
  mutate(
    datetime = parse_date_time(
      datetime_MST,
      orders = c("Ymd HMS", "Y-m-d H:M:S", "Y-m-d H:M"),
      tz = "America/Denver"
    )
  ) %>%
  # Drop rows that failed to parse (avoids downstream errors)
  filter(!is.na(datetime)) %>%
  # Bin timestamps into hourly "buckets" (e.g., 10:00–10:59 -> 10:00)
  mutate(datetime_hr = floor_date(datetime, "hour")) %>%
  # For each hour, compute mean temperature at each depth
  group_by(datetime_hr) %>%
  summarize(
    var_1 = mean(temp_1m, na.rm = TRUE),  # hourly mean at 1 m
    var_2 = mean(temp_2m, na.rm = TRUE),  # hourly mean at 2 m
    var_3 = mean(temp_3m, na.rm = TRUE),  # hourly mean at 3 m
    var_4 = mean(temp_4m, na.rm = TRUE),  # hourly mean at 4 m
    .groups = "drop"
  ) %>%
  # Rename datetime_hr back to datetime so the first column is "datetime"
  rename(datetime = datetime_hr) %>%
  # Remove hours with any missing depth values (strict; may drop lots of rows)
  drop_na()

# ---- 3) Try plotting directly (may fail depending on rLakeAnalyzer version) ----
# Some installs of wtr.heat.map() are picky about input object type.
# If it works for you, great. If it errors, use the load.ts() route below.
wtr.heat.map(ub_heat)

# ---- 4) The most reliable rLakeAnalyzer workflow: write -> load.ts() -> plot ----
# Make a copy and ensure datetime is explicitly POSIXct
ub_heat2 <- ub_heat
ub_heat2$datetime <- as.POSIXct(ub_heat2$datetime, tz = "America/Denver")

# load.ts() expects the datetime column to be a nicely formatted string in the file
ub_out <- ub_heat2
ub_out$datetime <- format(ub_out$datetime, "%Y-%m-%d %H:%M:%S")

# Create a temporary file path (won't clutter your folders)
tmp <- tempfile(fileext = ".txt")

# Write a tab-delimited text file in the exact format load.ts() expects
write.table(
  ub_out, tmp,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Read the file back in as an rLakeAnalyzer time series object
ts_ub <- load.ts(tmp)

# Plot the rLakeAnalyzer heatmap (this should be the dependable one)
wtr.heat.map(ts_ub)






####### mixed and strat 

# ub_heat has columns: datetime, var_1, var_2, var_3, var_4 (hourly means)

ub_mix <- ub_heat %>%
  mutate(
    dT_1_4 = abs(var_1 - var_4),      # surface-bottom temp difference
    dT_1_2 = abs(var_1 - var_2),
    dT_2_3 = abs(var_2 - var_3),
    dT_3_4 = abs(var_3 - var_4),
    max_grad = pmax(dT_1_2, dT_2_3, dT_3_4, na.rm = TRUE),  # steepest layer
    state = case_when(
      dT_1_4 < 0.5 ~ "Mixed",
      dT_1_4 > 1.0 ~ "Stratified",
      TRUE ~ "Transition"
    )
  )

table(ub_mix$state)


ggplot(ub_mix, aes(datetime, dT_1_4)) +
  geom_line() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  labs(
    title = "Upper Brooks mixing indicator (ΔT = |T1m - T4m|)",
    x = NULL,
    y = "ΔT (°C)"
  ) +
  theme_bw()

ggplot(ub_mix, aes(datetime, 1, fill = state)) +
  geom_tile() +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "Upper Brooks mixing state (hourly)", x = NULL, fill = NULL) +
  theme_bw()