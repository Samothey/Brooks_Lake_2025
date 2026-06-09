library(tidyverse)
library(lubridate)
library(rLakeAnalyzer)


# Load data
brooks_wq <- read_csv("~/Desktop/Brooks_Modeling/data_clean/brooks_wq_2025_cleaned.csv",
                      show_col_types = FALSE) %>%
  select(-...1) %>%  # drop index column
  mutate(datetime = as.POSIXct(DATETIME_MST, tz = "America/Denver")) %>%
  filter(!is.na(datetime))

# Hourly temps in var_depth format
brooks_heat <- brooks_wq %>%
  mutate(datetime_hr = floor_date(datetime, "hour")) %>%
  group_by(datetime_hr) %>%
  summarize(
    var_1  = mean(WTemp1m,  na.rm = TRUE),
    var_4  = mean(WTemp4m,  na.rm = TRUE),
    var_7  = mean(WTemp7m,  na.rm = TRUE),
    var_9  = mean(WTemp9m,  na.rm = TRUE),
    var_10 = mean(WTemp10m, na.rm = TRUE),
    var_13 = mean(WTemp13m, na.rm = TRUE),
    var_15 = mean(WTemp15m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(datetime = datetime_hr) %>%
  filter(if_any(starts_with("var_"), ~ !is.na(.x)))

# Convert to load.ts() object via temp file
brooks_out <- brooks_heat %>%
  mutate(datetime = format(as.POSIXct(datetime, tz = "America/Denver"),
                           "%Y-%m-%d %H:%M:%S"))

tmp <- tempfile(fileext = ".txt")
write.table(brooks_out, tmp, sep = "\t", row.names = FALSE, quote = FALSE)

ts_brooks <- load.ts(tmp)

# Plot heatmap
wtr.heat.map(ts_brooks)


# How many NAs per depth column?
sapply(brooks_heat %>% select(starts_with("var_")), \(x) sum(is.na(x)))

# Zoom in on October and see missingness by depth
brooks_heat %>%
  filter(datetime >= as.POSIXct("2025-10-01", tz="America/Denver")) %>%
  summarize(across(starts_with("var_"), ~ sum(is.na(.x))))











##########################################
# ============================================================
# Brooks buoy temperature heatmap (Option 2: fill blanks)
# - Hourly means by depth (var_1, var_4, ...)
# - Complete missing hours on the time axis
# - Fill missing depth values via interpolation (short gaps only)
# - Export to rLakeAnalyzer load.ts() format + plot wtr.heat.map()
# ============================================================

library(tidyverse)
library(lubridate)
library(zoo)
library(rLakeAnalyzer)

# ---- 1) Load cleaned data ----
brooks_wq <- read_csv(
  "~/Desktop/Brooks_Modeling/data_clean/brooks_wq_2025_cleaned.csv",
  show_col_types = FALSE
) %>%
  select(-...1) %>%  # drop index col if present
  mutate(
    datetime = as.POSIXct(DATETIME_MST, tz = "America/Denver")
  ) %>%
  filter(!is.na(datetime))

# ---- 2) Hourly means in var_depth format (ensure NA not NaN) ----
brooks_heat <- brooks_wq %>%
  mutate(datetime_hr = floor_date(datetime, "hour")) %>%
  group_by(datetime_hr) %>%
  summarize(
    var_1  = if (all(is.na(WTemp1m )))  NA_real_ else mean(WTemp1m,  na.rm = TRUE),
    var_4  = if (all(is.na(WTemp4m )))  NA_real_ else mean(WTemp4m,  na.rm = TRUE),
    var_7  = if (all(is.na(WTemp7m )))  NA_real_ else mean(WTemp7m,  na.rm = TRUE),
    var_9  = if (all(is.na(WTemp9m )))  NA_real_ else mean(WTemp9m,  na.rm = TRUE),
    var_10 = if (all(is.na(WTemp10m)))  NA_real_ else mean(WTemp10m, na.rm = TRUE),
    var_13 = if (all(is.na(WTemp13m)))  NA_real_ else mean(WTemp13m, na.rm = TRUE),
    var_15 = if (all(is.na(WTemp15m)))  NA_real_ else mean(WTemp15m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(datetime = datetime_hr) %>%
  arrange(datetime)

# ---- 3) Ensure a complete hourly timeline (prevents missing-hour stripes) ----
brooks_heat <- brooks_heat %>%
  complete(datetime = seq(min(datetime), max(datetime), by = "1 hour")) %>%
  arrange(datetime)

# ---- 4) Fill blanks (interpolate short gaps only) ----
# maxgap = 6 means: only fill runs of <= 6 missing hours.
# If you want to fill longer gaps, bump this up (e.g., 12 or 24).
brooks_heat_filled <- brooks_heat %>%
  mutate(across(starts_with("var_"),
                ~ na.approx(.x, x=datetime, na.rm=FALSE, rule=2)))

# (Optional) If you want *no blanks at all*, uncomment this:
# brooks_heat_filled <- brooks_heat %>%
#   mutate(across(starts_with("var_"),
#                 ~ na.approx(.x, x = datetime, na.rm = FALSE, rule = 2)))

# ---- 5) Convert datetime to character + write temp file for load.ts() ----
brooks_out <- brooks_heat_filled %>%
  mutate(datetime = format(as.POSIXct(datetime, tz = "America/Denver"),
                           "%Y-%m-%d %H:%M:%S"))

tmp <- tempfile(fileext = ".txt")
write.table(brooks_out, tmp, sep = "\t", row.names = FALSE, quote = FALSE)

# ---- 6) Load + plot heatmap ----
ts_brooks <- load.ts(tmp)
wtr.heat.map(ts_brooks)

# ---- 7) Quick QA checks (optional) ----
cat("\nNAs per depth AFTER fill (maxgap=6):\n")
print(sapply(brooks_heat_filled %>% select(starts_with("var_")), \(x) sum(is.na(x))))

cat("\nHours where ALL depths are missing (should be 0 unless full outage):\n")
print(brooks_heat_filled %>%
        summarize(all_depths_missing = sum(if_all(starts_with("var_"), is.na))))