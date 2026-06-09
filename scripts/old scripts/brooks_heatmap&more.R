library(tidyverse)
library(lubridate)
library(scico)
library(rLakeAnalyzer)
library(zoo)

theme_set(theme_bw(base_family = "Helvetica"))

# Load data ---------------------------------------------------------------

brooks <- read_csv("/Users/samanthapena/Desktop/Project/Brooks_lake_2025/data_clean/buoy/brooks_wq_2025_cleaned.csv")
view(brooks)
colnames(brooks)


brooks <- brooks |>
  select(-...1) |>
  rename(WTemp0.75m = WTempC) |>
  mutate(
    datetime = DATETIME_MST
  )

# Convert temperature data to long format --------------------------------

wtemp <- brooks |>
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

# Summarize to hourly means -----------------------------------------------

wtemp_hourly <- wtemp |>
  mutate(datetime_hour = round_date(datetime, unit = "hour")) |>
  group_by(datetime_hour, depth) |>
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# Interpolate between sensor depths --------------------------------------

interp_hourly <- wtemp_hourly |>
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

# Plot heatmap ------------------------------------------------------------

brooks_heatmap <- ggplot(interp_hourly, aes(x = datetime_hour, y = depth, fill = temp_c)) +
  geom_tile(
    width = 3600,
    height = 0.2
  ) +
  scale_y_reverse(
    breaks = seq(0, 15, by = 2)
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
    ylim = c(14.9, 0.6),
    expand = FALSE
  ) +
  labs(
    x = "Date",
    y = "Depth (m)",
    title = ""
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

brooks_heatmap


# Create output folder if it doesn't exist
brooks_heatmap <- here("figures/heatmaps")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}
ggsave(
  filename = file.path(out_path, "brooks_temperature_heatmap_2025.pdf"),
  plot = brooks_heatmap,
  width = 10,
  height = 5.5,
  units = "in"
)


# Create temperature wide table for Schmidt stability ----------------------

brooks_temp_wide <- wtemp_hourly |>
  mutate(
    depth_name = paste0("wtr_", depth)
  ) |>
  select(datetime = datetime_hour, depth_name, temp_c) |>
  pivot_wider(
    names_from = depth_name,
    values_from = temp_c
  )

# Temperature columns in correct depth order ------------------------------

temp_cols <- c(
  "wtr_0.75",
  "wtr_1",
  "wtr_4",
  "wtr_7",
  "wtr_9",
  "wtr_10",
  "wtr_13",
  "wtr_15"
)

temp_depths <- c(0.75, 1, 4, 7, 9, 10, 13, 15)

brooks_temp_wide_full <- brooks_temp_wide |>
  select(datetime, all_of(temp_cols)) |>
  filter(if_all(all_of(temp_cols), ~ !is.na(.x)))
# Estimated Brooks bathymetry ---------------------------------------------

bthD <- c(0.75, 1, 4, 7, 9, 10, 13, 15)

bthA <- c(
  866030,
  820000,
  650000,
  430000,
  300000,
  230000,
  90000,
  10000
)

# Calculate Schmidt stability --------------------------------------------

schmidt_hourly <- brooks_temp_wide_full |>
  rowwise() |>
  mutate(
    schmidt_stability = schmidt.stability(
      wtr = c_across(all_of(temp_cols)),
      depths = temp_depths,
      bthD = bthD,
      bthA = bthA
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
    ),
    strat_state = case_when(
      stability_24hr <= 25 ~ "Mixed/Weak",
      stability_24hr <= 50 ~ "Moderate",
      stability_24hr > 50 ~ "Strong",
      TRUE ~ NA_character_
    ),
    strat_state = factor(
      strat_state,
      levels = c("Mixed/Weak", "Moderate", "Strong")
    )
  )


# Plot hourly raw + smoothed stability ------------------------------------------
sampling_dates <- as.POSIXct(paste0(c(
  "2025-06-23",
  "2025-07-08",
  "2025-07-22",
  "2025-08-05",
  "2025-08-19",
  "2025-09-02",
  "2025-09-16",
  "2025-09-30",
  "2025-10-14"
) ," 12:00:00"),  tz = "America/Denver")


## plotting stability with colored lines 
brooks_stability_plot <- ggplot(schmidt_hourly, aes(x = datetime_hour)) +
  geom_line(
    aes(y = schmidt_stability),
    linewidth = 0.35,
    alpha = 0.2,
    color = "gray60"
  ) +
  geom_line(
    aes(y = stability_24hr, color = strat_state),
    linewidth = 1.1,
    color = "black"
  ) +
  geom_point(
    aes(y = stability_24hr, color = strat_state),
    size = 1.2,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = sampling_dates,
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate" = "#fdae61",
      "Strong" = "#d7191c"
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

brooks_stability_plot



## temporal trend in DO over time 


# Dissolved oxygen --------------------------------------------------------

do_brooks <- brooks |>
  select(datetime, odomgL, odomgL_15m) |>
  pivot_longer(
    cols = c(odomgL, odomgL_15m),
    names_to = "depth",
    values_to = "do_mgl"
  ) |>
  mutate(
    depth = case_when(
      depth == "odomgL" ~ "Surface",
      depth == "odomgL_15m" ~ "15 m",
      TRUE ~ depth
    )
  ) |>
  filter(!is.na(do_mgl))

do_hourly <- do_brooks |>
  mutate(datetime_hour = round_date(datetime, unit = "hour")) |>
  group_by(datetime_hour, depth) |>
  summarise(
    do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

## time series
brooks_do_plot <- ggplot(do_hourly, aes(x = datetime_hour, y = do_mgl, color = depth)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(
    values = c(
      "Surface" = "#1f78b4",
      "15 m" = "#d73027"
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

brooks_do_plot


## join DO and stability 
deep_do_hourly <- do_hourly |>
  filter(depth == "15 m") |>
  select(datetime_hour, deep_do_mgl = do_mgl)

deep_do_strat_hourly <- deep_do_hourly |>
  left_join(schmidt_hourly, by = "datetime_hour") |>
  filter(!is.na(strat_state), !is.na(deep_do_mgl))


# Deep DO hourly, boxplot
deep_do_boxplot_hourly <- ggplot(
  deep_do_strat_hourly,
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
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"   = "#fdae61",
      "Strong"     = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Stratification state",
    y = "Deep dissolved oxygen (mg/L)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none"
  )

deep_do_boxplot_hourly


## just to compare lets look at daily strat state and daily DO
# Daily deep DO ------------------------------------------------------------

deep_do_daily <- do_hourly |>
  filter(depth == "15 m") |>
  mutate(date = as.Date(datetime_hour)) |>
  group_by(date) |>
  summarise(
    deep_do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

# Daily stability state ---------------------------------------------------

schmidt_daily <- schmidt_hourly |>
  mutate(date = as.Date(datetime_hour)) |>
  filter(!is.na(stability_24hr)) |>
  group_by(date) |>
  summarise(
    stability_daily = mean(stability_24hr, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    strat_state = case_when(
      stability_daily <= 25 ~ "Mixed/Weak",
      stability_daily <= 50 ~ "Moderate",
      stability_daily > 50  ~ "Strong",
      TRUE ~ NA_character_
    ),
    strat_state = factor(
      strat_state,
      levels = c("Mixed/Weak", "Moderate", "Strong")
    )
  )

# Join daily DO + stratification state ------------------------------------

deep_do_strat_daily <- deep_do_daily |>
  left_join(schmidt_daily, by = "date") |>
  filter(!is.na(strat_state), !is.na(deep_do_mgl))

# Check sample sizes -------------------------------------------------------

table(deep_do_strat_daily$strat_state)

deep_do_boxplot_daily <- ggplot(
  deep_do_strat_daily,
  aes(x = strat_state, y = deep_do_mgl)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = strat_state),
    width = 0.15,
    size = 2,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"   = "#fdae61",
      "Strong"     = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Stratification state",
    y = "Daily mean deep DO (mg/L)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

deep_do_boxplot_daily



#### now bring ins nutrients
## load nutrients

deq_nutrients_clean_2025 <- readRDS("~/Desktop/Project/Brooks_lake_2025/data_clean/deq/deq_nutrients_clean_2025.rds"
                                    
nutrients_brooks <- deq_nutrients_clean_2025 |>
  filter(str_to_lower(lake) == "brooks lake") |>  # adjust if exact name differs
  mutate(
    date = as.Date(date),
    depth_zone = case_when(
      str_detect(str_to_lower(depth), "surf") ~ "Surface",
      str_detect(str_to_lower(depth), "mid") ~ "Middle / thermocline",
      str_detect(str_to_lower(depth), "therm") ~ "Middle / thermocline",
      str_detect(str_to_lower(depth), "bot") ~ "Deep",
      TRUE ~ NA_character_
    ),
    depth_zone = factor(
      depth_zone,
      levels = c("Surface", "Middle / thermocline", "Deep")
    )
  ) |>
  filter(!is.na(depth_zone)) |>
  filter(date >= as.Date("2025-04-01"))

nutrients_long <- nutrients_brooks |>
  select(date, depth_zone, ammonia, tp, tn) |>
  pivot_longer(
    cols = c(ammonia,tp, tn),
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


bottom_nutrients <- nutrients_long |>
  filter(depth_zone == "Deep") |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily, by = "date") |>
  filter(!is.na(strat_state), !is.na(concentration))

bottom_nutrient_boxplot <- ggplot(
  bottom_nutrients,
  aes(x = strat_state, y = concentration)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = strat_state),
    width = 0.15,
    size = 2,
    alpha = 0.75
  ) +
  facet_wrap(~ nutrient, scales = "free_y") +
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"   = "#fdae61",
      "Strong"     = "#d7191c"
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

bottom_nutrient_boxplot

##  ## plot nutrients with stability 
# Deep nutrients joined to daily stability --------------------------------

deep_nutrients_overlay <- nutrients_long |>
  filter(depth_zone == "Deep") |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily, by = "date") |>
  filter(!is.na(strat_state), !is.na(concentration))

# Scale nutrient concentration to plot on stability axis -------------------

deep_nutrients_overlay <- deep_nutrients_overlay |>
  group_by(nutrient) |>
  mutate(
    concentration_scaled = concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  ) |>
  ungroup()

ammonia_overlay <- deep_nutrients_overlay |>
  filter(nutrient == "Ammonia")

ggplot() +
  geom_line(
    data = schmidt_daily,
    aes(x = date, y = stability_daily),
    linewidth = 1,
    color = "black"
  ) +
  geom_point(
    data = ammonia_overlay,
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
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"   = "#fdae61",
      "Strong"     = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Date",
    y = "Schmidt stability (J m⁻²)",
    color = "Stability state",
    size = "Ammonia (µg/L)"
  ) +
  theme_classic(base_family = "Helvetica")

 
 
 
 