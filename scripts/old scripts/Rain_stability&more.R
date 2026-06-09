# ============================================================
# RAINBOW LAKE: Temperature heatmap, Schmidt stability,
# DO, and nutrient comparison
# ============================================================

library(tidyverse)
library(lubridate)
library(rLakeAnalyzer)
library(zoo)
library(janitor)
library(scico)
library(here)

# ------------------------------------------------------------
# 1. Clean Rainbow buoy data
# ------------------------------------------------------------

rainbow_clean <- rainbow |>
  select(-...1) |>
  rename(WTemp0.75m = WTempC) |>
  mutate(
    datetime = DATETIME_MST
  )

# ------------------------------------------------------------
# 2. Convert temperature data to long format
# ------------------------------------------------------------

wtemp_rn <- rainbow_clean |>
  select(datetime, matches("^WTemp")) |>
  pivot_longer(
    cols = matches("^WTemp"),
    names_to = "depth",
    values_to = "temp_c"
  ) |>
  mutate(
    depth = str_remove(depth, "WTemp"),
    depth = str_remove(depth, "m"),
    depth = as.numeric(depth)
  ) |>
  filter(!is.na(temp_c))

# ------------------------------------------------------------
# 3. Summarize to hourly means
# ------------------------------------------------------------

wtemp_hourly_rn <- wtemp_rn |>
  mutate(datetime_hour = round_date(datetime, unit = "hour")) |>
  group_by(datetime_hour, depth) |>
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4. Interpolate between sensor depths for heatmap
# ------------------------------------------------------------

interp_hourly_rn <- wtemp_hourly_rn |>
  group_by(datetime_hour) |>
  nest() |>
  mutate(n_depths = map_dbl(data, nrow)) |>
  filter(n_depths > 1) |>
  mutate(data = map(data, ~ arrange(.x, depth))) |>
  mutate(
    interp_fun = map(data, ~ approxfun(.x$depth, .x$temp_c)),
    raster = map2(data, interp_fun, function(df, func) {
      tibble(
        depth = seq(min(df$depth), max(df$depth), by = 0.2),
        temp_c = func(depth)
      )
    })
  ) |>
  select(datetime_hour, raster) |>
  unnest(raster) |>
  filter(!is.na(temp_c))

# ------------------------------------------------------------
# 5. Plot Rainbow temperature heatmap
# ------------------------------------------------------------

rainbow_heatmap <- ggplot(
  interp_hourly_rn,
  aes(x = datetime_hour, y = depth, fill = temp_c)
) +
  geom_tile(width = 3600, height = 0.2) +
  scale_y_reverse(breaks = seq(0, 8, by = 1)) +
  scale_fill_scico(
    palette = "romaO",
    limits = c(4, 19),
    oob = scales::squish
  ) +
  scale_x_datetime(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  coord_cartesian(
    ylim = c(7.2, 0.6),
    expand = FALSE
  ) +
  labs(
    x = "Date",
    y = "Depth (m)",
    fill = "Temp (°C)",
    title = ""
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

rainbow_heatmap

# Optional save
out_path <- here("figures/heatmaps")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

ggsave(
  filename = file.path(out_path, "rainbow_temperature_heatmap_2025.pdf"),
  plot = rainbow_heatmap,
  width = 10,
  height = 5.5,
  units = "in"
)

# ------------------------------------------------------------
# 6. Create temperature wide table for Schmidt stability
# ------------------------------------------------------------

rainbow_temp_wide <- wtemp_hourly_rn |>
  mutate(depth_name = paste0("wtr_", depth)) |>
  select(datetime = datetime_hour, depth_name, temp_c) |>
  pivot_wider(
    names_from = depth_name,
    values_from = temp_c
  )

temp_cols_rn <- c(
  "wtr_0.75",
  "wtr_1",
  "wtr_2",
  "wtr_3",
  "wtr_4",
  "wtr_5.5",
  "wtr_7"
)

temp_depths_rn <- c(0.75, 1, 2, 3, 4, 5.5, 7)

rainbow_temp_wide_full <- rainbow_temp_wide |>
  select(datetime, all_of(temp_cols_rn)) |>
  filter(if_all(all_of(temp_cols_rn), ~ !is.na(.x)))

# ------------------------------------------------------------
# 7. Rainbow bathymetry
# ------------------------------------------------------------
# Surface area = 31 acres
# Max depth = 8 m
# Areas are estimated and converted from acres to m2.

acre_to_m2 <- 4046.856

bthD_rn <- c(0.75, 1, 2, 3, 4, 5.5, 7, 8)

bthA_rn <- c(
  31,
  30,
  27,
  23,
  17,
  10,
  3,
  0.5
) * acre_to_m2

# ------------------------------------------------------------
# 8. Calculate Schmidt stability
# ------------------------------------------------------------

schmidt_hourly_rn <- rainbow_temp_wide_full |>
  rowwise() |>
  mutate(
    schmidt_stability = schmidt.stability(
      wtr = c_across(all_of(temp_cols_rn)),
      depths = temp_depths_rn,
      bthD = bthD_rn,
      bthA = bthA_rn
    )
  ) |>
  ungroup() |>
  select(datetime_hour = datetime, schmidt_stability) |>
  mutate(
    stability_24hr = rollmean(
      schmidt_stability,
      k = 24,
      fill = NA,
      align = "center"
    )
  )

# ------------------------------------------------------------
# 9. Define lake-specific relative stability states
# ------------------------------------------------------------

schmidt_hourly_rn <- schmidt_hourly_rn |>
  mutate(
    strat_state = case_when(
      stability_24hr <= quantile(stability_24hr, 0.33, na.rm = TRUE) ~ "Low stability",
      stability_24hr <= quantile(stability_24hr, 0.66, na.rm = TRUE) ~ "Moderate stability",
      stability_24hr > quantile(stability_24hr, 0.66, na.rm = TRUE) ~ "High stability",
      TRUE ~ NA_character_
    ),
    strat_state = factor(
      strat_state,
      levels = c("Low stability", "Moderate stability", "High stability")
    )
  )

# ------------------------------------------------------------
# 10. Plot hourly raw + smoothed stability
# ------------------------------------------------------------

sampling_dates_rn <- as.POSIXct(
  paste0(c(
    "2025-07-08"
  ), " 12:00:00"),
  tz = "America/Denver"
)

rainbow_stability_plot <- ggplot(schmidt_hourly_rn, aes(x = datetime_hour)) +
  geom_line(
    aes(y = schmidt_stability),
    linewidth = 0.35,
    alpha = 0.2,
    color = "gray60"
  ) +
  geom_line(
    aes(y = stability_24hr),
    linewidth = 1.1,
    color = "black"
  ) +
  geom_point(
    aes(y = stability_24hr, color = strat_state),
    size = 1.2,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = sampling_dates_rn,
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Low stability" = "#2c7bb6",
      "Moderate stability" = "#fdae61",
      "High stability" = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  scale_x_datetime(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  labs(
    x = "Date",
    y = expression("Schmidt stability (J m"^{-2}*")"),
    color = "Stability state",
    title = ""
  ) +
  theme_classic(base_family = "Helvetica")

rainbow_stability_plot

# ------------------------------------------------------------
# 11. Dissolved oxygen
# ------------------------------------------------------------

do_rn <- rainbow_clean |>
  select(datetime, odomgL, odomgL_7m) |>
  pivot_longer(
    cols = c(odomgL, odomgL_7m),
    names_to = "depth",
    values_to = "do_mgl"
  ) |>
  mutate(
    depth = case_when(
      depth == "odomgL" ~ "Surface",
      depth == "odomgL_7m" ~ "7 m",
      TRUE ~ depth
    )
  ) |>
  filter(!is.na(do_mgl))

do_hourly_rn <- do_rn |>
  mutate(datetime_hour = round_date(datetime, unit = "hour")) |>
  group_by(datetime_hour, depth) |>
  summarise(
    do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

rainbow_do_plot <- ggplot(
  do_hourly_rn,
  aes(x = datetime_hour, y = do_mgl, color = depth)
) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(
    values = c(
      "Surface" = "#1f78b4",
      "7 m" = "#d73027"
    )
  ) +
  scale_y_continuous(limits = c(0, 18)) +
  scale_x_datetime(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  labs(
    x = "Date",
    y = "Dissolved oxygen (mg/L)",
    color = "Depth"
  ) +
  theme_classic(base_family = "Helvetica")

rainbow_do_plot

# ------------------------------------------------------------
# 12. Join bottom DO and stability - hourly
# ------------------------------------------------------------

bottom_do_hourly_rn <- do_hourly_rn |>
  filter(depth == "7 m") |>
  select(datetime_hour, bottom_do_mgl = do_mgl)

bottom_do_strat_hourly_rn <- bottom_do_hourly_rn |>
  left_join(schmidt_hourly_rn, by = "datetime_hour") |>
  filter(!is.na(strat_state), !is.na(bottom_do_mgl))

bottom_do_boxplot_hourly_rn <- ggplot(
  bottom_do_strat_hourly_rn,
  aes(x = strat_state, y = bottom_do_mgl)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(
    aes(color = strat_state),
    width = 0.15,
    size = 1.5,
    alpha = 0.55
  ) +
  scale_color_manual(
    values = c(
      "Low stability" = "#2c7bb6",
      "Moderate stability" = "#fdae61",
      "High stability" = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Relative stability state",
    y = "Bottom DO at 7 m (mg/L)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position = "none")

bottom_do_boxplot_hourly_rn

# ------------------------------------------------------------
# 13. Daily bottom DO and daily stability
# ------------------------------------------------------------

bottom_do_daily_rn <- do_hourly_rn |>
  filter(depth == "7 m") |>
  mutate(date = as.Date(datetime_hour)) |>
  group_by(date) |>
  summarise(
    bottom_do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

schmidt_daily_rn <- schmidt_hourly_rn |>
  mutate(date = as.Date(datetime_hour)) |>
  filter(!is.na(stability_24hr)) |>
  group_by(date) |>
  summarise(
    stability_daily = mean(stability_24hr, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    strat_state = case_when(
      stability_daily <= quantile(stability_daily, 0.33, na.rm = TRUE) ~ "Low stability",
      stability_daily <= quantile(stability_daily, 0.66, na.rm = TRUE) ~ "Moderate stability",
      stability_daily > quantile(stability_daily, 0.66, na.rm = TRUE) ~ "High stability",
      TRUE ~ NA_character_
    ),
    strat_state = factor(
      strat_state,
      levels = c("Low stability", "Moderate stability", "High stability")
    )
  )

bottom_do_strat_daily_rn <- bottom_do_daily_rn |>
  left_join(schmidt_daily_rn, by = "date") |>
  filter(!is.na(strat_state), !is.na(bottom_do_mgl))

table(bottom_do_strat_daily_rn$strat_state)

bottom_do_boxplot_daily_rn <- ggplot(
  bottom_do_strat_daily_rn,
  aes(x = strat_state, y = bottom_do_mgl)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(
    aes(color = strat_state),
    width = 0.15,
    size = 2,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = c(
      "Low stability" = "#2c7bb6",
      "Moderate stability" = "#fdae61",
      "High stability" = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Relative stability state",
    y = "Daily mean bottom DO at 7 m (mg/L)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

bottom_do_boxplot_daily_rn

# ------------------------------------------------------------
# 14. Bring in nutrients
# ------------------------------------------------------------

deq_nutrients_clean_2025 <- readRDS(
  "~/Desktop/Project/Brooks_lake_2025/data_clean/deq/deq_nutrients_clean_2025.rds"
)

nutrients_rn <- deq_nutrients_clean_2025 |>
  clean_names() |>
  filter(str_to_lower(lake) %in% c("rainbow", "rainbow lake")) |>
  mutate(
    date = as.Date(date),
    depth_zone = case_when(
      str_detect(str_to_lower(depth), "surface|surf") ~ "Surface",
      str_detect(str_to_lower(depth), "bottom|bot|deep") ~ "Deep",
      TRUE ~ NA_character_
    ),
    depth_zone = factor(
      depth_zone,
      levels = c("Surface", "Deep")
    )
  ) |>
  filter(!is.na(depth_zone)) |>
  filter(date >= as.Date("2025-04-01"))

nutrients_long_rn <- nutrients_rn |>
  select(date, depth_zone, ammonia, tp, tn) |>
  pivot_longer(
    cols = c(ammonia, tp, tn),
    names_to = "nutrient",
    values_to = "concentration"
  ) |>
  filter(!is.na(concentration)) |>
  group_by(date, depth_zone, nutrient) |>
  summarise(
    concentration = mean(concentration, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    nutrient = recode(
      nutrient,
      ammonia = "Ammonia",
      tp = "Total Phosphorus",
      tn = "Total Nitrogen"
    )
  )

# ------------------------------------------------------------
# 15. Bottom nutrients by stability state
# ------------------------------------------------------------

bottom_nutrients_rn <- nutrients_long_rn |>
  filter(depth_zone == "Deep") |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily_rn, by = "date") |>
  filter(!is.na(strat_state), !is.na(concentration))

bottom_nutrient_boxplot_rn <- ggplot(
  bottom_nutrients_rn,
  aes(x = strat_state, y = concentration)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(
    aes(color = strat_state),
    width = 0.15,
    size = 2,
    alpha = 0.75
  ) +
  facet_wrap(~ nutrient, scales = "free_y") +
  scale_color_manual(
    values = c(
      "Low stability" = "#2c7bb6",
      "Moderate stability" = "#fdae61",
      "High stability" = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Relative stability state",
    y = "Bottom nutrient concentration (µg/L)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

bottom_nutrient_boxplot_rn