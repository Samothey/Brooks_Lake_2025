# ============================================================
# UPPER BROOKS: Temperature heatmap, Schmidt stability,
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
# 1. Clean Upper Brooks data
# ------------------------------------------------------------

upper_brooks <- upper_brooks_2025_clean |>
  clean_names() |>
  mutate(
    datetime = as.POSIXct(datetime_mst, tz = "America/Denver")
  )

# ------------------------------------------------------------
# 2. Convert temperature data to long format
# ------------------------------------------------------------

wtemp_ub <- upper_brooks |>
  select(datetime, temp_1m, temp_2m, temp_3m, temp_4m) |>
  pivot_longer(
    cols = starts_with("temp_"),
    names_to = "depth",
    values_to = "temp_c"
  ) |>
  mutate(
    depth = str_remove(depth, "temp_"),
    depth = str_remove(depth, "m"),
    depth = as.numeric(depth)
  ) |>
  filter(!is.na(temp_c))

# ------------------------------------------------------------
# 3. Summarize to hourly means
# ------------------------------------------------------------

wtemp_hourly_ub <- wtemp_ub |>
  mutate(datetime_hour = round_date(datetime, unit = "hour")) |>
  group_by(datetime_hour, depth) |>
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 4. Interpolate between sensor depths for heatmap
# ------------------------------------------------------------

interp_hourly_ub <- wtemp_hourly_ub |>
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
# 5. Plot Upper Brooks temperature heatmap
# ------------------------------------------------------------

upper_brooks_heatmap <- ggplot(
  interp_hourly_ub,
  aes(x = datetime_hour, y = depth, fill = temp_c)
) +
  geom_tile(
    width = 3600,
    height = 0.2
  ) +
  scale_y_reverse(
    breaks = seq(0, 5, by = 1)
  ) +
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
    ylim = c(4.2, 0.8),
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

upper_brooks_heatmap

# Optional save
out_path <- here("figures/heatmaps")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

ggsave(
  filename = file.path(out_path, "upper_brooks_temperature_heatmap_2025.pdf"),
  plot = upper_brooks_heatmap,
  width = 10,
  height = 5.5,
  units = "in"
)

# ------------------------------------------------------------
# 6. Create temperature wide table for Schmidt stability
# ------------------------------------------------------------

upper_brooks_temp_wide <- wtemp_hourly_ub |>
  mutate(
    depth_name = paste0("wtr_", depth)
  ) |>
  select(datetime = datetime_hour, depth_name, temp_c) |>
  pivot_wider(
    names_from = depth_name,
    values_from = temp_c
  )

temp_cols_ub <- c(
  "wtr_1",
  "wtr_2",
  "wtr_3",
  "wtr_4"
)

temp_depths_ub <- c(1, 2, 3, 4)

upper_brooks_temp_wide_full <- upper_brooks_temp_wide |>
  select(datetime, all_of(temp_cols_ub)) |>
  filter(if_all(all_of(temp_cols_ub), ~ !is.na(.x)))

# ------------------------------------------------------------
# 7. Upper Brooks bathymetry
# ------------------------------------------------------------
# Surface area = 24 acres
# Max depth = 5 m
# schmidt.stability() needs area in m2

acre_to_m2 <- 4046.856

bthD_ub <- c(0, 1, 2, 3, 4, 5)

bthA_ub <- c(
  24,
  22,
  18,
  12,
  6,
  0.5
) * acre_to_m2


# ------------------------------------------------------------
# 8. Calculate Schmidt stability
# ------------------------------------------------------------

schmidt_hourly_ub <- upper_brooks_temp_wide_full |>
  rowwise() |>
  mutate(
    schmidt_stability = schmidt.stability(
      wtr = c_across(all_of(temp_cols_ub)),
      depths = temp_depths_ub,
      bthD = bthD_ub,
      bthA = bthA_ub
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

schmidt_hourly_ub <- schmidt_hourly_ub |>
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
# 9. Plot hourly raw + smoothed stability
# ------------------------------------------------------------

sampling_dates_ub <- as.POSIXct(
  paste0(c(
    "2025-06-23",
    "2025-07-08",
    "2025-07-22",
    "2025-08-05",
    "2025-08-19",
    "2025-09-02",
    "2025-09-16",
    "2025-09-30",
    "2025-10-14"
  ), " 12:00:00"),
  tz = "America/Denver"
)

upper_brooks_stability_plot <- ggplot(schmidt_hourly_ub, aes(x = datetime_hour)) +
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
    xintercept = sampling_dates_ub,
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

upper_brooks_stability_plot

# ------------------------------------------------------------
# 10. Dissolved oxygen
# ------------------------------------------------------------

do_ub <- upper_brooks |>
  select(datetime, do_mgl_1m, do_mgl_4m) |>
  pivot_longer(
    cols = c(do_mgl_1m, do_mgl_4m),
    names_to = "depth",
    values_to = "do_mgl"
  ) |>
  mutate(
    depth = case_when(
      depth == "do_mgl_1m" ~ "Surface",
      depth == "do_mgl_4m" ~ "4 m",
      TRUE ~ depth
    )
  ) |>
  filter(!is.na(do_mgl))

do_hourly_ub <- do_ub |>
  mutate(datetime_hour = round_date(datetime, unit = "hour")) |>
  group_by(datetime_hour, depth) |>
  summarise(
    do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

upper_brooks_do_plot <- ggplot(
  do_hourly_ub,
  aes(x = datetime_hour, y = do_mgl, color = depth)
) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(
    values = c(
      "Surface" = "#1f78b4",
      "4 m" = "#d73027"
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

upper_brooks_do_plot

# ------------------------------------------------------------
# 11. Join deep DO and stability - hourly
# ------------------------------------------------------------

deep_do_hourly_ub <- do_hourly_ub |>
  filter(depth == "4 m") |>
  select(datetime_hour, deep_do_mgl = do_mgl)

deep_do_strat_hourly_ub <- deep_do_hourly_ub |>
  left_join(schmidt_hourly_ub, by = "datetime_hour") |>
  filter(!is.na(strat_state), !is.na(deep_do_mgl))

deep_do_boxplot_hourly_ub <- ggplot(
  deep_do_strat_hourly_ub,
  aes(x = strat_state, y = deep_do_mgl)
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
    x = "Stratification state",
    y = "Deep dissolved oxygen (mg/L)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(legend.position = "none")

deep_do_boxplot_hourly_ub


# ------------------------------------------------------------
# 12. Daily DO at 4 m and daily stability
# ------------------------------------------------------------

deep_do_daily_ub <- do_hourly_ub |>
  filter(depth == "4 m") |>
  mutate(date = as.Date(datetime_hour)) |>
  group_by(date) |>
  summarise(
    deep_do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

schmidt_daily_ub <- schmidt_hourly_ub |>
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

deep_do_strat_daily_ub <- deep_do_daily_ub |>
  left_join(schmidt_daily_ub, by = "date") |>
  filter(!is.na(strat_state), !is.na(deep_do_mgl))

table(deep_do_strat_daily_ub$strat_state)

deep_do_boxplot_daily_ub <- ggplot(
  deep_do_strat_daily_ub,
  aes(x = strat_state, y = deep_do_mgl)
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
    y = "Daily mean DO at 4 m (mg/L)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

deep_do_boxplot_daily_ub
# ------------------------------------------------------------
# 13. Bring in nutrients
# ------------------------------------------------------------

nutrients_ub <- deq_nutrients_clean_2025 |>
  clean_names() |>
  filter(str_to_lower(lake) %in% c("upper brooks lake")) |>
  mutate(
    date = as.Date(date),
    depth_zone = case_when(
      str_detect(str_to_lower(depth), "surf") ~ "Surface",
      str_detect(str_to_lower(depth), "mid") ~ "Middle / thermocline",
      str_detect(str_to_lower(depth), "therm") ~ "Middle / thermocline",
      str_detect(str_to_lower(depth), "bot") ~ "Deep",
      str_detect(str_to_lower(depth), "deep") ~ "Deep",
      TRUE ~ NA_character_
    ),
    depth_zone = factor(
      depth_zone,
      levels = c("Surface", "Middle / thermocline", "Deep")
    )
  ) |>
  filter(!is.na(depth_zone)) |>
  filter(date >= as.Date("2025-04-01"))

nutrients_long_ub <- nutrients_ub |>
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
# 14. Bottom nutrients by stratification state
# ------------------------------------------------------------

bottom_nutrients_ub <- nutrients_long_ub |>
  filter(depth_zone == "Deep") |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily_ub, by = "date") |>
  filter(!is.na(strat_state), !is.na(concentration))

bottom_nutrient_boxplot_ub <- ggplot(
  bottom_nutrients_ub,
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
    x = "Stratification state",
    y = "Deep nutrient concentration (µg/L)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

bottom_nutrient_boxplot_ub

# ------------------------------------------------------------
# 15. Deep nutrients joined to daily stability
# ------------------------------------------------------------

deep_nutrients_overlay_ub <- nutrients_long_ub |>
  filter(depth_zone == "Deep") |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily_ub, by = "date") |>
  filter(!is.na(strat_state), !is.na(concentration))

# ------------------------------------------------------------
# 16. Ammonia overlay
# ------------------------------------------------------------

ammonia_overlay_ub <- deep_nutrients_overlay_ub |>
  filter(nutrient == "Ammonia") |>
  mutate(
    concentration_scaled =
      concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily_ub$stability_daily, na.rm = TRUE)
  )

ggplot() +
  geom_line(
    data = schmidt_daily_ub,
    aes(x = date, y = stability_daily),
    linewidth = 1,
    color = "black"
  ) +
  geom_point(
    data = ammonia_overlay_ub,
    aes(
      x = date,
      y = concentration_scaled,
      color = strat_state,
      size = concentration
    ),
    alpha = 0.85
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
    x = "Date",
    y = expression("Schmidt stability (J m"^{-2}*")"),
    color = "Stability state",
    size = "Ammonia (µg/L)"
  ) +
  theme_classic(base_family = "Helvetica")

# ------------------------------------------------------------
# 17. Total nitrogen overlay
# ------------------------------------------------------------

tn_overlay_ub <- deep_nutrients_overlay_ub |>
  filter(nutrient == "Total Nitrogen") |>
  mutate(
    concentration_scaled =
      concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily_ub$stability_daily, na.rm = TRUE)
  )

ggplot() +
  geom_line(
    data = schmidt_daily_ub,
    aes(x = date, y = stability_daily),
    linewidth = 1,
    color = "black"
  ) +
  geom_point(
    data = tn_overlay_ub,
    aes(
      x = date,
      y = concentration_scaled,
      color = strat_state,
      size = concentration
    ),
    alpha = 0.85
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
    x = "Date",
    y = expression("Schmidt stability (J m"^{-2}*")"),
    color = "Stability state",
    size = "Total Nitrogen (µg/L)"
  ) +
  theme_classic(base_family = "Helvetica")



## Tp Overlay
tp_overlay_ub <- deep_nutrients_overlay_ub |>
  filter(nutrient == "Total Phosphorus") |>
  mutate(
    concentration_scaled =
      concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily_ub$stability_daily, na.rm = TRUE)
  )

ggplot() +
  geom_line(
    data = schmidt_daily_ub,
    aes(x = date, y = stability_daily),
    linewidth = 1,
    color = "black"
  ) +
  geom_point(
    data = tp_overlay_ub,
    aes(
      x = date,
      y = concentration_scaled,
      color = strat_state,
      size = concentration
    ),
    alpha = 0.85
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
    x = "Date",
    y = expression("Schmidt stability (J m"^{-2}*")"),
    color = "Stability state",
    size = "Total Phosphorus (µg/L)"
  ) +
  theme_classic(base_family = "Helvetica")
