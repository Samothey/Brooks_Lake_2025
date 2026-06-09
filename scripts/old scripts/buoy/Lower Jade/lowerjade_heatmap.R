library(tidyverse)     # dplyr + ggplot + tidyr, etc.
library(lubridate)     # parse_date_time(), floor_date()
library(rLakeAnalyzer) # load.ts(), wtr.heat.map()

# ---- 1) Load the cleaned Lower Jade buoy dataset ----
lower_jade <- read_csv(
  "/Users/samanthapena/Desktop/Project/Brooks_lake_2025/DATA/BUOY/Lower Jade/cleaned/LowerJade_buoy_2025.csv",
  show_col_types = FALSE
)

# ---- 2) Convert 10-minute data to hourly means in "var_depth" format ----
# rLakeAnalyzer likes a "wide" format where:
#   - first column is datetime
#   - remaining columns are named var_1, var_2, ... where the number = depth (m)
lj_heat <- lower_jade %>%
  # Robust datetime parsing (handles "YYYY-mm-dd HH:MM:SS" and "YYYY-mm-dd HH:MM")
  mutate(
    datetime = parse_date_time(
      str_trim(datetime),
      orders = c("Y-m-d H:M:S", "Y-m-d H:M"),
      tz = "America/Denver"
    )
  ) %>%
  filter(!is.na(datetime)) %>%
  mutate(datetime_hr = floor_date(datetime, "hour")) %>%
  group_by(datetime_hr) %>%
  summarize(
    var_1  = mean(temp_1m,  na.rm = TRUE),
    var_2  = mean(temp_2m,  na.rm = TRUE),
    var_4  = mean(temp_4m,  na.rm = TRUE),
    var_5  = mean(temp_5m,  na.rm = TRUE),
    var_6  = mean(temp_6m,  na.rm = TRUE),
    var_7  = mean(temp_7m,  na.rm = TRUE),
    var_8  = mean(temp_8m,  na.rm = TRUE),
    var_9  = mean(temp_9m,  na.rm = TRUE),
    var_10 = mean(temp_10m, na.rm = TRUE),
    var_11 = mean(temp_11m, na.rm = TRUE),
    var_12 = mean(temp_12m, na.rm = TRUE),
    var_13 = mean(temp_13m, na.rm = TRUE),
    var_14 = mean(temp_14m, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(datetime = datetime_hr) %>%
  # Strict option like your template:
  drop_na()

# ---- 3) Try plotting directly (may work depending on rLakeAnalyzer version) ----
wtr.heat.map(lj_heat)

# ---- 4) Reliable workflow: write -> load.ts() -> plot ----
lj_heat2 <- lj_heat
lj_heat2$datetime <- as.POSIXct(lj_heat2$datetime, tz = "America/Denver")

lj_out <- lj_heat2
lj_out$datetime <- format(lj_out$datetime, "%Y-%m-%d %H:%M:%S")

tmp <- tempfile(fileext = ".txt")

write.table(
  lj_out, tmp,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ts_lj <- load.ts(tmp)
wtr.heat.map(ts_lj)





####### mixed and strat (Lower Jade)

# lj_heat has columns: datetime, var_1, var_2, var_4, ..., var_14 (hourly means)

lj_mix <- lj_heat %>%
  mutate(
    dT_1_14 = abs(var_1 - var_14),  # surface-bottom temp difference (1m vs 14m)
    
    # optional: within-column gradients (coarse because we skip 3m, 15m, etc.)
    dT_1_2  = abs(var_1 - var_2),
    dT_2_4  = abs(var_2 - var_4),
    dT_4_5  = abs(var_4 - var_5),
    dT_5_6  = abs(var_5 - var_6),
    dT_6_7  = abs(var_6 - var_7),
    dT_7_8  = abs(var_7 - var_8),
    dT_8_9  = abs(var_8 - var_9),
    dT_9_10 = abs(var_9 - var_10),
    dT_10_11 = abs(var_10 - var_11),
    dT_11_12 = abs(var_11 - var_12),
    dT_12_13 = abs(var_12 - var_13),
    dT_13_14 = abs(var_13 - var_14),
    
    max_grad = pmax(
      dT_1_2, dT_2_4, dT_4_5, dT_5_6, dT_6_7, dT_7_8,
      dT_8_9, dT_9_10, dT_10_11, dT_11_12, dT_12_13, dT_13_14,
      na.rm = TRUE
    ),
    
    state = case_when(
      dT_1_14 < 0.5 ~ "Mixed",
      dT_1_14 > 1.0 ~ "Stratified",
      TRUE ~ "Transition"
    )
  )

table(lj_mix$state)

ggplot(lj_mix, aes(datetime, dT_1_14)) +
  geom_line() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  labs(
    title = "Lower Jade mixing indicator (ΔT = |T1m - T14m|)",
    x = NULL,
    y = "ΔT (°C)"
  ) +
  theme_bw()

ggplot(lj_mix, aes(datetime, 1, fill = state)) +
  geom_tile() +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "Lower Jade mixing state (hourly)", x = NULL, fill = NULL) +
  theme_bw()