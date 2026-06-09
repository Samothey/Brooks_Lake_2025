# ======================================================================
# Rainbow Lake temperature and DO structure plots
#
# Purpose:
#   1. Plot Rainbow buoy temperature data while the full thermistor
#      string was active.
#   2. Remove incomplete hours where only the 7 m miniDOT recorded.
#   3. Plot DEQ vertical profile temperature and DO snapshots for the
#      full season.
#
# Notes:
#   - Rainbow buoy thermistor string failed later in July.
#   - The 7 m miniDOT continued recording, but the full temperature string
#     did not.
#   - For the buoy heatmap, we keep only hours where all 7 temperature
#     depths are present.
#   - Removed hours are compressed out of the x-axis using a time index.
# ======================================================================

library(tidyverse)
library(lubridate)
library(janitor)
library(readxl)
library(scico)
library(scales)

# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis/rainbow"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load Rainbow buoy data
# ======================================================================

Rainbow_WQ_2025_cleaned <- read_csv(
  "data_clean/buoy/Rainbow_WQ_2025_cleaned.csv",
  show_col_types = FALSE
)

rainbow_buoy <- Rainbow_WQ_2025_cleaned %>%
  clean_names() %>%
  mutate(
    datetime = datetime_mst
  )

# ======================================================================
# 3. Buoy temperature data: temperature string only
# ======================================================================

rainbow_temp_hourly <- rainbow_buoy %>%
  select(
    datetime,
    w_temp_c,
    w_temp1m,
    w_temp2m,
    w_temp3m,
    w_temp4m,
    w_temp5_5m,
    w_temp7m
  ) %>%
  pivot_longer(
    cols = starts_with("w_temp"),
    names_to = "depth_raw",
    values_to = "temp_c"
  ) %>%
  mutate(
    depth = case_when(
      depth_raw == "w_temp_c" ~ 0.75,
      depth_raw == "w_temp1m" ~ 1,
      depth_raw == "w_temp2m" ~ 2,
      depth_raw == "w_temp3m" ~ 3,
      depth_raw == "w_temp4m" ~ 4,
      depth_raw == "w_temp5_5m" ~ 5.5,
      depth_raw == "w_temp7m" ~ 7,
      TRUE ~ NA_real_
    ),
    datetime_hour = floor_date(datetime, unit = "hour")
  ) %>%
  filter(
    !is.na(datetime_hour),
    !is.na(depth),
    !is.na(temp_c),
    as.Date(datetime_hour) <= as.Date("2025-07-17")
  ) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 4. Keep only complete buoy hours
# ======================================================================

complete_hours <- rainbow_temp_hourly %>%
  count(datetime_hour, name = "n_depths") %>%
  filter(n_depths == 7) %>%
  arrange(datetime_hour) %>%
  mutate(
    time_index = row_number()
  )

rainbow_temp_plot <- rainbow_temp_hourly %>%
  semi_join(
    complete_hours,
    by = "datetime_hour"
  ) %>%
  left_join(
    complete_hours %>%
      select(datetime_hour, time_index),
    by = "datetime_hour"
  ) %>%
  mutate(
    depth = factor(
      depth,
      levels = rev(c(0.75, 1, 2, 3, 4, 5.5, 7))
    )
  )

# ======================================================================
# 5. QA checks for buoy temperature data
# ======================================================================

rainbow_temp_plot %>%
  summarise(
    first_time = min(datetime_hour),
    last_time = max(datetime_hour),
    n_complete_hours = n_distinct(datetime_hour),
    n_depths = n_distinct(depth),
    min_temp = min(temp_c, na.rm = TRUE),
    max_temp = max(temp_c, na.rm = TRUE)
  )

# Hours removed because not all 7 depths were present
rainbow_temp_hourly %>%
  count(datetime_hour, name = "n_depths") %>%
  filter(n_depths != 7)

# Raw timestamp check
rainbow_buoy %>%
  arrange(datetime) %>%
  mutate(
    dt_diff = as.numeric(
      difftime(
        datetime,
        lag(datetime),
        units = "mins"
      )
    )
  ) %>%
  count(dt_diff)

# Depth-specific observation counts
rainbow_temp_hourly %>%
  group_by(depth) %>%
  summarise(
    n_hours = n(),
    first_time = min(datetime_hour),
    last_time = max(datetime_hour),
    .groups = "drop"
  )

# ======================================================================
# 6. Plot Rainbow buoy temperature heatmap
# ======================================================================

date_breaks <- complete_hours %>%
  mutate(date = as.Date(datetime_hour)) %>%
  filter(
    date %in% seq(
      min(date),
      max(date),
      by = "3 days"
    )
  ) %>%
  group_by(date) %>%
  slice(1) %>%
  ungroup()

p_rainbow_buoy_temp <- ggplot(
  rainbow_temp_plot,
  aes(
    x = time_index,
    y = depth,
    fill = temp_c
  )
) +
  geom_tile(
    width = 1,
    height = 0.95
  ) +
  scale_fill_gradientn(
    colours = c(
      "#313695",
      "#4575B4",
      "#74ADD1",
      "#ABD9E9",
      "#FEE090",
      "#FDAE61",
      "#F46D43",
      "#D73027",
      "#A50026"
    ),
    limits = c(4, 19),
    oob = squish
  ) +
  scale_x_continuous(
    breaks = date_breaks$time_index,
    labels = format(date_breaks$datetime_hour, "%b %d"),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Depth (m)",
    fill = "Temp (°C)",
    title = "Rainbow buoy temperature record"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_rainbow_buoy_temp

ggsave(
  file.path(fig_dir, "rainbow_buoy_temperature_complete_hours.png"),
  p_rainbow_buoy_temp,
  width = 10,
  height = 4,
  dpi = 300
)

# ======================================================================
# 7. Load Rainbow DEQ vertical profile data
# ======================================================================

lake_profile_compiled <- read_excel(
  "data_clean/deq/lake_profile_compiled.xlsx"
)

rainbow_profiles <- lake_profile_compiled %>%
  mutate(
    date = as.Date(date),
    lake = str_to_lower(lake),
    depth = as.numeric(depth),
    temp = as.numeric(temp),
    do_mg_l = as.numeric(do_mg_l)
  ) %>%
  filter(
    lake %in% c("rainbow", "rainbow lake"),
    !is.na(depth)
  ) %>%
  arrange(date, depth)

# ======================================================================
# 8. Rainbow DEQ profile temperature stripes
# ======================================================================

p_rainbow_profile_temp <- ggplot(
  rainbow_profiles %>%
    filter(!is.na(temp)),
  aes(
    x = date,
    y = depth,
    fill = temp
  )
) +
  geom_tile(
    width = 2,
    height = 0.8
  ) +
  scale_y_reverse(
    breaks = seq(0, 8, by = 1)
  ) +
  scale_fill_scico(
    palette = "vik",
    direction = 1,
    limits = c(4, 19),
    oob = scales::squish
  ) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %d"
  ) +
  labs(
    x = NULL,
    y = "Depth (m)",
    fill = "Temp (°C)",
    title = "Rainbow profile temperature snapshots"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_rainbow_profile_temp

ggsave(
  file.path(fig_dir, "rainbow_profile_temperature_stripes.png"),
  p_rainbow_profile_temp,
  width = 10,
  height = 5,
  dpi = 300
)

# ======================================================================
# 9. Rainbow DEQ profile DO stripes
# ======================================================================

p_rainbow_profile_do <- ggplot(
  rainbow_profiles %>%
    filter(!is.na(do_mg_l)),
  aes(
    x = date,
    y = depth,
    fill = do_mg_l
  )
) +
  geom_tile(
    width = 2,
    height = 0.8
  ) +
  scale_y_reverse(
    breaks = seq(0, 8, by = 1)
  ) +
  scale_fill_scico(
    palette = "lajolla",
    limits = c(0, 16),
    oob = scales::squish
  ) +
  scale_x_date(
    date_breaks = "2 weeks",
    date_labels = "%b %d"
  ) +
  labs(
    x = NULL,
    y = "Depth (m)",
    fill = "DO (mg/L)",
    title = "Rainbow profile dissolved oxygen snapshots"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p_rainbow_profile_do

ggsave(
  file.path(fig_dir, "rainbow_profile_do_stripes.png"),
  p_rainbow_profile_do,
  width = 10,
  height = 5,
  dpi = 300
)



