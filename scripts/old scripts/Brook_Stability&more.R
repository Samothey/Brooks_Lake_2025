library(tidyverse)
library(lubridate)
library(scico)
library(rLakeAnalyzer)
library(zoo)
library(janitor)

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
    stability_24hr = as.numeric(
      rollmean(
      schmidt_stability,
      k = 24,
      fill = NA,
      align = "center"
      )
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

## plot without colroed lines 
brooks_stability_plot <- ggplot(
  schmidt_hourly,
  aes(x = datetime_hour)
) +
  
  # raw hourly
  geom_line(
    aes(y = schmidt_stability),
    linewidth = 0.35,
    alpha = 0.2,
    color = "gray60"
  ) +
  
  # smoothed 24 hr
  geom_line(
    aes(y = stability_24hr),
    linewidth = 1.1,
    color = "black"
  ) +
  
  geom_vline(
    xintercept = sampling_dates,
    linetype = "dashed",
    color = "gray50",
    alpha = 0.6
  ) +
  
  scale_x_datetime(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  
  labs(
    x = "Date",
    y = expression("Schmidt stability (J m"^{-2}*")")
  ) +
  
  theme_classic(base_family = "Helvetica")

brooks_stability_plot


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

deq_nutrients_clean_2025 <- readRDS("~/Desktop/Project/Brooks_lake_2025/data_clean/deq/deq_nutrients_clean_2025.rds")
                                    
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
##  nutriets with stability 

# Deep nutrients joined to daily stability --------------------------------

deep_nutrients_overlay <- nutrients_long |>
  filter(depth_zone == "Deep") |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily, by = "date") |>
  filter(!is.na(strat_state), !is.na(concentration))

## overlay for ammonia
ammonia_overlay <- deep_nutrients_overlay |>
  filter(nutrient == "Ammonia") |>
  mutate(
    concentration_scaled =
      concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  )
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
      "Mixed/Weak (0-25)" = "#2c7bb6",
      "Moderate (25-50)"   = "#fdae61",
      "Strong (>50) "     = "#d7191c"
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


## tn overlay 
tn_overlay <- deep_nutrients_overlay |>
  filter(nutrient == "Total Nitrogen") |>
  mutate(
    concentration_scaled =
      concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  )
ggplot() +
  geom_line(
    data = schmidt_daily,
    aes(x = date, y = stability_daily),
    linewidth = 1,
    color = "black"
  ) +
  geom_point(
    data = tn_overlay,
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
    size = "Total Nitrogen (µg/L)"
  ) +
  theme_classic(base_family = "Helvetica")

## Tp Overlay 
tp_overlay <- deep_nutrients_overlay |>
  filter(nutrient == "Total Phosphorus") |>
  mutate(
    concentration_scaled =
      concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  )
ggplot() +
  geom_line(
    data = schmidt_daily,
    aes(x = date, y = stability_daily),
    linewidth = 1,
    color = "black"
  ) +
  geom_point(
    data = tp_overlay,
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
    size = "Total Phosphorus (µg/L)"
  ) +
  theme_classic(base_family = "Helvetica")


## toxin stuff

grab_class_totals <- readRDS("~/Desktop/Project/Brooks_lake_2025/data_clean/toxins/grab_class_totals.rds")

grab_toxins_brooks <- grab_class_totals |>
  mutate(
    date = as.Date(date),
    lake = str_to_lower(lake),
    site_type = str_to_lower(site_type),
    depth_category = str_to_lower(depth_category),
    toxin_class = str_to_lower(toxin_class)
  ) |>
  filter(
    lake == "brooks",
    method == "grab",
    site_type == "buoy",
    depth_category %in% c("surface", "bottom"),
    date >= as.Date("2025-04-01")
  ) |>
  left_join(schmidt_daily, by = "date") |>
  filter(!is.na(strat_state), !is.na(total))


## boxplot
bottom_toxins <- grab_toxins_brooks |>
  filter(depth_category == "bottom")

bottom_toxin_boxplot <- ggplot(
  bottom_toxins,
  aes(x = strat_state, y = total)
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
  facet_wrap(~ toxin_class, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
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
    y = "Bottom toxin concentration, total class sum",
    title = "Bottom grab toxins by stratification state"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
## overlay with strat
bottom_toxin_boxplot
bottom_toxin_overlay <- grab_toxins_brooks |>
  filter(
    depth_category == "bottom",
    toxin_class == "microcystin"
  ) |>
  mutate(
    total_scaled =
      total / max(total, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  )

bottom_mc_overlay_plot <- ggplot() +
  geom_line(
    data = schmidt_daily,
    aes(x = date, y = stability_daily),
    linewidth = 1,
    color = "black"
  ) +
  geom_point(
    data = bottom_toxin_overlay,
    aes(
      x = date,
      y = total_scaled,
      color = strat_state,
      size = total
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
    size = "Bottom microcystin"
  ) +
  theme_classic(base_family = "Helvetica")

bottom_mc_overlay_plot

# surface and depth buoy
toxin_depth_difference <- grab_toxins_brooks |>
  filter(toxin_class == "microcystin") |>
  group_by(date, depth_category, strat_state, stability_daily) |>
  summarise(
    total = mean(total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = depth_category,
    values_from = total
  ) |>
  filter(!is.na(surface), !is.na(bottom)) |>
  mutate(
    bottom_minus_surface = bottom - surface
  )

depth_difference_plot <- ggplot(
  toxin_depth_difference,
  aes(x = stability_daily, y = bottom_minus_surface)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(color = strat_state),
    size = 3,
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
    x = "Schmidt stability (J m⁻²)",
    y = "Bottom - surface microcystin",
    color = "Stability state",
    title = "Vertical toxin difference across stability gradient"
  ) +
  theme_classic(base_family = "Helvetica")

depth_difference_plot

## surface and depth adn stability. 
surface_bottom_toxins_scaled <- grab_toxins_brooks |>
  filter(toxin_class == "microcystin") |>
  group_by(date, depth_category, strat_state, stability_daily) |>
  summarise(
    total = mean(total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    toxin_scaled =
      total / max(total, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  )

surface_bottom_stability_plot <- ggplot() +
  geom_line(
    data = schmidt_daily,
    aes(x = date, y = stability_daily),
    color = "black",
    linewidth = 1
  ) +
  geom_line(
    data = surface_bottom_toxins_scaled,
    aes(x = date, y = toxin_scaled, color = depth_category, group = depth_category),
    linewidth = 1
  ) +
  geom_point(
    data = surface_bottom_toxins_scaled,
    aes(
      x = date,
      y = toxin_scaled,
      color = depth_category,
      shape = strat_state,
      size = total
    ),
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "surface" = "#1b9e77",
      "bottom"  = "#7570b3"
    )
  ) +
  labs(
    x = "Date",
    y = "Schmidt stability, with toxin scaled to stability axis",
    color = "Depth",
    shape = "Stratification state",
    size = "Microcystin",
    title = "Surface and bottom microcystin overlaid on Schmidt stability"
  ) +
  theme_classic(base_family = "Helvetica")

surface_bottom_stability_plot
colnames(spatt_presence)
## spatt
spatt_richness <- spatt_presence |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily, by = "date") |>
  filter(
    lake == "brooks",
    site_type == "buoy",
    depth_category %in% c("surface", "bottom")
  ) |>
  group_by(
    date,
    depth_category,
    strat_state
  ) |>
  summarise(
    toxin_classes_detected = sum(detected, na.rm = TRUE),
    .groups = "drop"
  )
ggplot(
  spatt_richness,
  aes(
    x = date,
    y = toxin_classes_detected,
    color = depth_category
  )
) +
  geom_line(linewidth = 1) +
  geom_point(
    aes(shape = strat_state),
    size = 3
  ) +
  scale_y_continuous(
    breaks = 0:6
  ) +
  labs(
    x = "Date",
    y = "Number of toxin classes detected",
    color = "Depth",
    shape = "Stratification state",
    title = "SPATT toxin detections through time"
  ) +
  theme_classic(base_family = "Helvetica")
## plot 
spatt_richness_scaled <- spatt_richness |>
  mutate(
    richness_scaled =
      toxin_classes_detected /
      max(toxin_classes_detected, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  )

ggplot() +
  geom_line(
    data = schmidt_daily,
    aes(x = date, y = stability_daily),
    color = "black",
    linewidth = 1
  ) +
  geom_line(
    data = spatt_richness_scaled,
    aes(
      x = date,
      y = richness_scaled,
      color = depth_category,
      group = depth_category
    ),
    linewidth = 1
  ) +
  geom_point(
    data = spatt_richness_scaled,
    aes(
      x = date,
      y = richness_scaled,
      color = depth_category,
      shape = strat_state,
      size = toxin_classes_detected
    ),
    alpha = 0.9
  ) +
  labs(
    x = "Date",
    y = "Schmidt stability (scaled overlay)",
    color = "Depth",
    shape = "Stratification state",
    size = "Toxin classes detected",
    title = "SPATT toxin richness and lake stability"
  ) +
  theme_classic(base_family = "Helvetica")

## plot

grab_boxplot_data <- grab_class_totals |>
  mutate(
    date = as.Date(date),
    lake = str_to_lower(lake),
    site_type = str_to_lower(site_type),
    depth_category = str_to_lower(depth_category),
    toxin_class = str_to_lower(toxin_class)
  ) |>
  filter(
    lake == "brooks",
    method == "grab",
    site_type == "buoy",
    depth_category %in% c("surface", "bottom"),
    toxin_class == "microcystin"
  ) |>
  left_join(schmidt_daily, by = "date") |>
  filter(!is.na(strat_state))
##plot 
surface_bottom_stability_boxplot <- ggplot(
  grab_boxplot_data,
  aes(
    x = depth_category,
    y = total,
    fill = strat_state
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.6,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    aes(color = strat_state),
    alpha = 0.8,
    size = 2,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    )
  ) +
  scale_y_continuous(trans = "log1p") +
  scale_fill_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"   = "#fdae61",
      "Strong"     = "#d7191c"
    )
  ) +
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"   = "#fdae61",
      "Strong"     = "#d7191c"
    )
  ) +
  labs(
    x = "Depth",
    y = "Microcystin concentration",
    fill = "Stability state",
    color = "Stability state",
    title = "Surface vs bottom microcystins across stratification states"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14)
  )

surface_bottom_stability_boxplot


# nutrietns overlapy.
# Deep nutrients joined to daily stability --------------------------------

deep_nutrients_overlay <- nutrients_long |>
  filter(depth_zone == "Deep") |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily, by = "date") |>
  filter(
    !is.na(strat_state),
    !is.na(concentration),
    nutrient %in% c("Ammonia", "Total Nitrogen", "Total Phosphorus")
  ) |>
  group_by(nutrient) |>
  mutate(
    concentration_scaled =
      concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  ) |>
  ungroup()

ggplot() +
  geom_line(
    data = schmidt_daily,
    aes(
      x = date,
      y = stability_daily,
      color = strat_state
    ),
    linewidth = 1.2
  ) +
  geom_point(
    data = deep_nutrients_overlay,
    aes(
      x = date,
      y = concentration_scaled,
      shape = nutrient,
      size = concentration
    ),
    color = "black",
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Mixed/Weak (0-25)" = "#2c7bb6",
      "Moderate (25-50)"  = "#fdae61",
      "Strong (>50)"      = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Ammonia" = 16,
      "Total Nitrogen" = 17,
      "Total Phosphorus" = 15
    )
  ) +
  labs(
    x = "Date",
    y = "Schmidt stability (J m⁻²)",
    color = "Stability state",
    shape = "Deep nutrient",
    size = "Nutrient concentration"
  ) +
  theme_classic(base_family = "Helvetica")




## try again, 
# Deep nutrients joined to daily stability -------------------------------

deep_nutrients_overlay <- nutrients_long |>
  filter(depth_zone == "Deep") |>
  mutate(date = as.Date(date)) |>
  left_join(schmidt_daily, by = "date") |>
  filter(
    !is.na(strat_state),
    !is.na(concentration),
    nutrient %in% c(
      "Ammonia",
      "Total Nitrogen",
      "Total Phosphorus"
    )
  ) |>
  group_by(nutrient) |>
  mutate(
    concentration_scaled =
      concentration / max(concentration, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  ) |>
  ungroup() |>
  arrange(date)

# make sure stability data is ordered
schmidt_daily_clean <- schmidt_daily |>
  arrange(date)

ggplot() +
  
  # stability line
  geom_line(
    data = schmidt_daily_clean,
    aes(
      x = date,
      y = stability_daily,
      color = strat_state,
      group = 1
    ),
    linewidth = 1.3,
    alpha = 0.9
  ) +
  
  # nutrient connecting lines
  geom_line(
    data = deep_nutrients_overlay,
    aes(
      x = date,
      y = concentration_scaled,
      group = nutrient,
      linetype = nutrient
    ),
    color = "grey40",
    linewidth = 0.7,
    alpha = 0.7
  ) +
  
  # nutrient points
  geom_point(
    data = deep_nutrients_overlay,
    aes(
      x = date,
      y = concentration_scaled,
      shape = nutrient,
      size = concentration
    ),
    color = "black",
    alpha = 0.85
  ) +
  
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"  = "#fdae61",
      "Strong"      = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  
  scale_shape_manual(
    values = c(
      "Ammonia" = 16,
      "Total Nitrogen" = 17,
      "Total Phosphorus" = 17
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "Ammonia" = "solid",
      "Total Nitrogen" = "dashed",
      "Total Phosphorus" = "dotted"
    )
  ) +
  
  labs(
    x = "Date",
    y = "Schmidt stability (J m⁻²)",
    color = "Stability state",
    shape = "Deep nutrient",
    linetype = "Deep nutrient",
    size = "Nutrient concentration"
  ) +
  
  theme_classic(base_family = "Helvetica")



surface_bottom_toxins_scaled <- grab_toxins_brooks |>
  filter(toxin_class == "microcystin") |>
  group_by(date, depth_category) |>
  summarise(
    total = mean(total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    toxin_scaled =
      total / max(total, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  ) |>
  arrange(date)

schmidt_daily_clean <- schmidt_daily |>
  mutate(date = as.Date(date)) |>
  arrange(date)

surface_bottom_stability_plot <- ggplot() +
  
  # colored stability line
  geom_line(
    data = schmidt_daily_clean,
    aes(
      x = date,
      y = stability_daily,
      color = strat_state,
      group = 1
    ),
    linewidth = 1.2
  ) +
  
  # toxin lines
  geom_line(
    data = surface_bottom_toxins_scaled,
    aes(
      x = date,
      y = toxin_scaled,
      group = depth_category,
      linetype = depth_category
    ),
    color = "grey35",
    linewidth = 0.8,
    alpha = 0.8
  ) +
  
  # toxin points colored by depth
  geom_point(
    data = surface_bottom_toxins_scaled,
    aes(
      x = date,
      y = toxin_scaled,
      color = depth_category,
      size = total
    ),
    alpha = 0.85
  ) +
  
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"  = "#fdae61",
      "Strong"      = "#d7191c",
      "surface" = "#1b9e77",
      "bottom"  = "#7570b3"
    )
  ) +
  
  labs(
    x = "Date",
    y = "Schmidt stability, with microcystin scaled to stability axis",
    color = "Variable",
    linetype = "Depth",
    size = "Microcystin",
    title = "Surface and bottom microcystin overlaid on Schmidt stability"
  ) +
  
  theme_classic(base_family = "Helvetica")

surface_bottom_stability_plot


surface_bottom_toxins_scaled <- grab_toxins_brooks |>
  filter(toxin_class == "microcystin") |>
  group_by(date, depth_category) |>
  summarise(
    total = mean(total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    toxin_scaled =
      total / max(total, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  ) |>
  arrange(date)

schmidt_daily_clean <- schmidt_daily |>
  mutate(date = as.Date(date)) |>
  arrange(date)

surface_bottom_stability_plot <- ggplot() +
  
  # colored stability line
  geom_line(
    data = schmidt_daily_clean,
    aes(
      x = date,
      y = stability_daily,
      color = strat_state,
      group = 1
    ),
    linewidth = 1.2,
    alpha = 0.9
  ) +
  
  # toxin lines
  geom_line(
    data = surface_bottom_toxins_scaled,
    aes(
      x = date,
      y = toxin_scaled,
      color = depth_category,
      group = depth_category
    ),
    linewidth = 1
  ) +
  
  # toxin points
  geom_point(
    data = surface_bottom_toxins_scaled,
    aes(
      x = date,
      y = toxin_scaled,
      color = depth_category,
      size = total
    ),
    alpha = 0.9
  ) +
  
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"  = "#fdae61",
      "Strong"      = "#d7191c",
      "surface"           = "#1b9e77",
      "bottom"            = "#7570b3"
    )
  ) +
  
  labs(
    x = "Date",
    y = "Schmidt stability, with microcystin scaled to stability axis",
    color = "Variable",
    size = "Microcystin",
    title = "Surface and bottom microcystin overlaid on Schmidt stability"
  ) +
  
  theme_classic(base_family = "Helvetica")

surface_bottom_stability_plot

spatt_richness_scaled <- spatt_richness |>
  mutate(
    richness_scaled =
      toxin_classes_detected /
      max(toxin_classes_detected, na.rm = TRUE) *
      max(schmidt_daily$stability_daily, na.rm = TRUE)
  ) |>
  arrange(date)

schmidt_daily_clean <- schmidt_daily |>
  mutate(date = as.Date(date)) |>
  arrange(date)

ggplot() +
  
  # colored stability line
  geom_line(
    data = schmidt_daily_clean,
    aes(
      x = date,
      y = stability_daily,
      color = strat_state,
      group = 1
    ),
    linewidth = 1.2,
    alpha = 0.9
  ) +
  
  # SPATT richness lines
  geom_line(
    data = spatt_richness_scaled,
    aes(
      x = date,
      y = richness_scaled,
      color = depth_category,
      group = depth_category
    ),
    linewidth = 1
  ) +
  
  # SPATT richness points
  geom_point(
    data = spatt_richness_scaled,
    aes(
      x = date,
      y = richness_scaled,
      color = depth_category,
      size = toxin_classes_detected
    ),
    alpha = 0.9
  ) +
  
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate"  = "#fdae61",
      "Strong"      = "#d7191c",
      "surface"           = "#1b9e77",
      "bottom"            = "#7570b3"
    )
  ) +
  
  labs(
    x = "Date",
    y = "Schmidt stability (scaled overlay)",
    color = "Variable",
    size = "Toxin classes detected",
    title = "SPATT toxin richness and lake stability"
  ) +
  
  theme_classic(base_family = "Helvetica")

# FIGURE 4: Time series — MC over time -------------------------------

grab_mc_ts <- grab_mc %>%
  mutate(
    sample_date = as.Date(date),
    
    depth_plot = case_when(
      site_type == "shore" ~ "Shore",
      site_type == "buoy" & depth_category == "surface" ~ "Buoy surface",
      site_type == "buoy" & depth_category == "bottom"  ~ "Buoy bottom",
      TRUE ~ NA_character_
    ),
    
    line_id = if_else(
      site_type == "buoy",
      paste(site_id, depth_plot, sep = "_"),
      site_id
    )
  ) %>%
  filter(!is.na(depth_plot))

p_ts <- ggplot(
  grab_mc_ts,
  aes(
    x = sample_date,
    y = total,
    group = line_id,
    color = depth_plot
  )
) +
  
  geom_line(
    linewidth = 0.9,
    alpha = 0.8
  ) +
  
  geom_point(
    size = 2.2,
    alpha = 0.9
  ) +
  
  geom_hline(
    yintercept = mc_threshold,
    color = "red",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  facet_wrap(~ lake, scales = "free_y") +
  
  scale_y_continuous(trans = "log1p") +
  
  scale_color_manual(
    values = c(
      "Shore" = "#1b9e77",
      "Buoy surface" = "#1f78b4",
      "Buoy bottom" = "#6a3d9a"
    )
  ) +
  
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    color = "",
    title = "Microcystins over time by lake (GRAB)",
    subtitle = "Lines represent individual sites; buoy split into surface vs bottom"
  ) +
  
  theme_classic(base_family = "Helvetica") +
  
  theme(
    strip.background = element_blank(),
    panel.spacing = unit(1, "lines")
  )

p_ts
# FIGURE: Brooks Lake MC time series by site ----------------------------

grab_mc_ts <- grab_class_totals %>%
  mutate(
    sample_date = as.Date(date),
    
    depth_plot = case_when(
      site_type == "shore" ~ "Shore",
      site_type == "buoy" & depth_category == "surface" ~ "Buoy surface",
      site_type == "buoy" & depth_category == "bottom"  ~ "Buoy bottom",
      TRUE ~ NA_character_
    ),
    
    # unique plotting ID
    line_id = case_when(
      site_type == "shore" ~ site_id,
      site_type == "buoy" ~ paste(site_id, depth_plot, sep = "_"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(depth_plot),
    str_to_lower(lake) == "brooks"
  )

p_ts <- ggplot(
  grab_mc_ts,
  aes(
    x = sample_date,
    y = total,
    group = line_id,
    color = line_id
  )
) +
  
  geom_line(
    linewidth = 1,
    alpha = 0.85
  ) +
  
  geom_point(
    size = 2.4,
    alpha = 0.9
  ) +
  
  geom_hline(
    yintercept = mc_threshold,
    color = "red",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  scale_y_continuous(trans = "log1p") +
  
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    color = "Site",
    title = "Brooks Lake microcystins over time (GRAB)",
    subtitle = "Individual site trajectories"
  ) +
  
  theme_classic(base_family = "Helvetica") +
  
  theme(
    legend.position = "right"
  )

p_ts


grab_mc_ts <- grab_class_totals %>%
  filter(
    str_to_lower(toxin_class) == "microcystin"
  ) %>%
  mutate(
    sample_date = as.Date(date),
    
    depth_plot = case_when(
      site_type == "shore" ~ "Shore",
      site_type == "buoy" & depth_category == "surface" ~ "Buoy surface",
      site_type == "buoy" & depth_category == "bottom"  ~ "Buoy bottom",
      TRUE ~ NA_character_
    ),
    
    # unique plotting ID
    line_id = case_when(
      site_type == "shore" ~ site_id,
      site_type == "buoy" ~ paste(site_id, depth_plot, sep = "_"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(depth_plot),
    str_to_lower(lake) == "brooks"
  )

p_ts <- ggplot(
  grab_mc_ts,
  aes(
    x = sample_date,
    y = total,
    group = line_id,
    color = line_id
  )
) +
  
  geom_line(
    linewidth = 1,
    alpha = 0.85
  ) +
  
  geom_point(
    size = 2.4,
    alpha = 0.9
  ) +
  
  geom_hline(
    yintercept = mc_threshold,
    color = "red",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  scale_y_continuous(trans = "log1p") +
  
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    color = "Site",
    title = "Brooks Lake microcystins over time (GRAB)",
    subtitle = "Individual site trajectories"
  ) +
  
  theme_classic(base_family = "Helvetica") +
  
  theme(
    legend.position = "right"
  )

p_ts



## phytoplanlton data 

library(readxl)
library(tidyverse)
library(lubridate)
library(janitor)
library(scales)

# ------------------------------------------------------------
# Clean phytoplankton data
# ------------------------------------------------------------
WYDEQ_2025_Lander_Phytoplankton_Heterocyst_Flat_Data <- read_excel("data_raw/phytoplankton/WYDEQ 2025 Lander Phytoplankton Heterocyst Flat Data.xlsx")
View(WYDEQ_2025_Lander_Phytoplankton_Heterocyst_Flat_Data)  

phyto_clean <- WYDEQ_2025_Lander_Phytoplankton_Heterocyst_Flat_Data |>
  clean_names() |>
  select(
    sample_station_name,
    sample_date_collected,
    division,
    taxon,
    per_taxa_id,
    heterocyst_count,
    ncu_counted,
    cells_counted,
    live_cells,
    dead_cells,
    ncu_counted_l,
    cells_counted_l,
    live_cells_l,
    dead_cells_l,
    total_cells_l,
    fields_of_view_analyzed,
    volume_analyzed_l,
    analyst,
    heterocyst_density_l
  ) |>
  rename(
    lake_site = sample_station_name,
    date = sample_date_collected,
    taxa_id = per_taxa_id,
    heterocysts = heterocyst_count,
    total_cells = total_cells_l,
    live_cells_density = live_cells_l,
    dead_cells_density = dead_cells_l,
    volume_analyzed = volume_analyzed_l,
    heterocyst_density = heterocyst_density_l
  ) |>
  mutate(
    date = mdy(date),
    
    lake_site = str_to_lower(lake_site),
    lake_site = str_replace_all(lake_site, "\\s*-\\s*", " - "),
    
    division = str_to_lower(division),
    taxon = str_squish(taxon),
    analyst = str_squish(analyst),
    
    total_cells = as.numeric(total_cells),
    live_cells_density = as.numeric(live_cells_density),
    dead_cells_density = as.numeric(dead_cells_density),
    heterocyst_density = as.numeric(heterocyst_density),
    
    lake = case_when(
      str_detect(lake_site, "^upper brooks lake") ~ "Upper Brooks Lake",
      str_detect(lake_site, "^brooks lake") ~ "Brooks Lake",
      str_detect(lake_site, "^rainbow lake") ~ "Rainbow Lake",
      str_detect(lake_site, "^lower jade lake") ~ "Lower Jade Lake",
      str_detect(lake_site, "^upper jade lake") ~ "Upper Jade Lake",
      TRUE ~ "Other"
    ),
    
    sample_type = case_when(
      str_detect(lake_site, "dup") ~ "duplicate",
      TRUE ~ "regular"
    )
  )

# ------------------------------------------------------------
# Brooks Lake only
# ------------------------------------------------------------

brooks_phyto <- phyto_clean |>
  filter(
    lake == "Brooks Lake",
    sample_type == "regular"
  )

unique(brooks_phyto$lake_site)
unique(brooks_phyto$date)

# ------------------------------------------------------------
# Community composition by division
# ------------------------------------------------------------

brooks_division_summary <- brooks_phyto |>
  group_by(date, division) |>
  summarize(
    total_cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(date) |>
  mutate(
    relative_abundance = total_cells / sum(total_cells, na.rm = TRUE)
  ) |>
  ungroup()

ggplot(
  brooks_division_summary,
  aes(x = date, y = relative_abundance, fill = division)
) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Date",
    y = "Relative abundance",
    fill = "Division",
    title = "Brooks Lake phytoplankton community composition"
  ) +
  theme_classic()

# ------------------------------------------------------------
# Cyanophyta only
# ------------------------------------------------------------

brooks_cyano <- brooks_phyto |>
  filter(division == "cyanophyta")

brooks_cyano_summary <- brooks_cyano |>
  group_by(date) |>
  summarize(
    cyano_cells_l = sum(total_cells, na.rm = TRUE),
    heterocyst_density_l = sum(heterocyst_density, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  brooks_cyano_summary,
  aes(x = date, y = cyano_cells_l)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_y_log10(labels = comma) +
  labs(
    x = "Date",
    y = "Cyanobacteria cells/L",
    title = "Brooks Lake cyanobacteria abundance through time"
  ) +
  theme_classic()

# ------------------------------------------------------------
# Dominant cyanobacteria taxa
# ------------------------------------------------------------

brooks_cyano_taxa <- brooks_cyano |>
  group_by(date, taxon) |>
  summarize(
    total_cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  )

top_cyano_taxa <- brooks_cyano_taxa |>
  group_by(taxon) |>
  summarize(
    total_cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  ) |>
  slice_max(total_cells, n = 8)

brooks_cyano_top <- brooks_cyano_taxa |>
  filter(taxon %in% top_cyano_taxa$taxon)

ggplot(
  brooks_cyano_top,
  aes(x = date, y = total_cells, color = taxon)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_y_log10(labels = comma) +
  labs(
    x = "Date",
    y = "Cells/L",
    color = "Taxon",
    title = "Dominant cyanobacteria taxa in Brooks Lake"
  ) +
  theme_classic()

## cyano by abundace 
brooks_cyano_relative_taxa <- brooks_cyano |>
  
  group_by(date, taxon) |>
  summarize(
    total_cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  ) |>
  
  group_by(date) |>
  mutate(
    relative_abundance =
      total_cells /
      sum(total_cells, na.rm = TRUE)
  ) |>
  
  ungroup()

top_relative_taxa <- brooks_cyano_relative_taxa |>
  
  group_by(taxon) |>
  
  summarize(
    mean_relative_abundance =
      mean(relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  
  slice_max(mean_relative_abundance, n = 8)
brooks_cyano_relative_top <- brooks_cyano_relative_taxa |>
  filter(taxon %in% top_relative_taxa$taxon)
ggplot(
  brooks_cyano_relative_top,
  aes(
    x = date,
    y = relative_abundance,
    fill = taxon
  )
) +
  
  geom_col() +
  
  scale_y_continuous(labels = percent_format()) +
  
  labs(
    x = "Date",
    y = "Relative abundance within cyanobacteria",
    fill = "Taxon",
    title = "Dominant cyanobacteria taxa in Brooks Lake"
  ) +
  
  theme_classic()


# ------------------------------------------------------------
# Cyanobacteria abundance by stratification state
# ------------------------------------------------------------

cyano_stability <- brooks_cyano_summary |>
  mutate(date = as.Date(date)) |>
  left_join(
    schmidt_daily |>
      select(date, stability_daily, strat_state),
    by = "date"
  ) |>
  filter(
    !is.na(strat_state),
    !is.na(cyano_cells_l)
  )

glimpse(cyano_stability)
table(cyano_stability$strat_state)

# ------------------------------------------------------------
# Boxplot: cyano abundance by stratification state
# ------------------------------------------------------------

cyano_stability_boxplot <- ggplot(
  cyano_stability,
  aes(x = strat_state, y = cyano_cells_l)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = strat_state),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  scale_y_continuous(
    trans = "log1p",
    labels = scales::comma
  ) +
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate" = "#fdae61",
      "Strong" = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Stratification state",
    y = "Cyanobacteria abundance (cells/L)",
    title = "Brooks Lake cyanobacteria abundance by stratification state"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

cyano_stability_boxplot
# ------------------------------------------------------------
# Model: cyano abundance as a function of stratification state
# ------------------------------------------------------------

cyano_stability_lm <- lm(
  log1p(cyano_cells_l) ~ strat_state,
  data = cyano_stability
)

summary(cyano_stability_lm)

cyano_stability_results <- broom::tidy(
  cyano_stability_lm,
  conf.int = TRUE
)

cyano_stability_results
cyano_stability_results_clean <- cyano_stability_results |>
  mutate(
    term = recode(
      term,
      `(Intercept)` = "Mixed/Weak",
      `strat_stateModerate` = "Moderate vs Mixed/Weak",
      `strat_stateStrong` = "Strong vs Mixed/Weak"
    )
  ) |>
  select(
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

cyano_stability_results_clean

# ------------------------------------------------------------
# Time series: cyano abundance colored by stratification state
# ------------------------------------------------------------

cyano_stability_timeseries <- ggplot(
  cyano_stability,
  aes(x = date, y = cyano_cells_l)
) +
  geom_line(
    linewidth = 1,
    color = "gray40"
  ) +
  geom_point(
    aes(color = strat_state),
    size = 3
  ) +
  scale_y_continuous(
    trans = "log1p",
    labels = scales::comma
  ) +
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate" = "#fdae61",
      "Strong" = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Date",
    y = "Cyanobacteria abundance (cells/L)",
    color = "Stratification state",
    title = "Brooks Lake cyanobacteria abundance through time"
  ) +
  theme_classic(base_family = "Helvetica")

cyano_stability_timeseries

###
# ------------------------------------------------------------
## ------------------------------------------------------------
# Overlay: stability + cyano abundance + surface/bottom MC
# ------------------------------------------------------------

schmidt_daily_clean <- schmidt_daily |>
  mutate(date = as.Date(date)) |>
  arrange(date)

# scaling factor so cyano can be shown on secondary axis
cyano_scale_factor <- max(schmidt_daily_clean$stability_daily, na.rm = TRUE) /
  max(brooks_cyano_summary$cyano_cells_l, na.rm = TRUE)

cyano_overlay <- brooks_cyano_summary |>
  mutate(
    date = as.Date(date),
    cyano_scaled = cyano_cells_l * cyano_scale_factor
  ) |>
  left_join(
    schmidt_daily_clean |>
      select(date, strat_state, stability_daily),
    by = "date"
  )

# surface and bottom microcystin, kept separate
mc_depth_overlay <- grab_toxins_brooks |>
  filter(toxin_class == "microcystin") |>
  group_by(date, depth_category) |>
  summarise(
    microcystin = mean(total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    schmidt_daily_clean |>
      select(date, strat_state, stability_daily),
    by = "date"
  ) |>
  filter(
    depth_category %in% c("surface", "bottom"),
    !is.na(microcystin)
  ) |>
  mutate(
    mc_scaled = microcystin /
      max(microcystin, na.rm = TRUE) *
      max(schmidt_daily_clean$stability_daily, na.rm = TRUE)
  )

stability_cyano_mc_overlay <- ggplot() +
  
  # stability line colored by state
  geom_line(
    data = schmidt_daily_clean,
    aes(
      x = date,
      y = stability_daily,
      color = strat_state,
      group = 1
    ),
    linewidth = 1.3
  ) +
  
  # cyano abundance scaled to stability axis
  geom_line(
    data = cyano_overlay,
    aes(x = date, y = cyano_scaled),
    color = "darkgreen",
    linewidth = 1
  ) +
  geom_point(
    data = cyano_overlay,
    aes(x = date, y = cyano_scaled),
    color = "darkgreen",
    size = 3,
    alpha = 0.85
  ) +
  
  # microcystin surface/bottom as shapes, size = concentration
  geom_point(
    data = mc_depth_overlay,
    aes(
      x = date,
      y = mc_scaled,
      shape = depth_category,
      size = microcystin
    ),
    color = "purple4",
    alpha = 0.85
  ) +
  
  scale_color_manual(
    values = c(
      "Mixed/Weak" = "#2c7bb6",
      "Moderate" = "#fdae61",
      "Strong" = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  
  scale_shape_manual(
    values = c(
      "surface" = 16,
      "bottom" = 17
    )
  ) +
  
  scale_size_continuous(
    range = c(2, 8)
  ) +
  
  scale_y_continuous(
    name = expression("Schmidt stability (J m"^{-2}*")"),
    sec.axis = sec_axis(
      ~ . / cyano_scale_factor,
      name = "Cyanobacteria abundance (cells/L)",
      labels = scales::comma
    )
  ) +
  
  labs(
    x = "Date",
    color = "Stability state",
    shape = "Microcystin depth",
    size = "Microcystin concentration",
    title = "Brooks Lake stability, cyanobacteria abundance, and microcystin",
    subtitle = "Microcystin point size reflects concentration; cyano abundance is shown on the secondary axis"
  ) +
  
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

stability_cyano_mc_overlay



##### models 
library(tidyverse)
library(broom)
library(gtsummary)

# ------------------------------------------------------------
# MODEL 1:
# Does cyanobacteria abundance vary by stability state?
# ThermalStability → CyanobacteriaAbundance
# ------------------------------------------------------------

model_cyano_stability <- lm(
  log1p(cyano_cells_l) ~ strat_state,
  data = cyano_stability
)

summary(model_cyano_stability)

tbl_cyano_stability <- tidy(
  model_cyano_stability,
  conf.int = TRUE
)

tbl_cyano_stability
# Clean table

tbl_cyano_stability_clean <- tbl_cyano_stability |>
  mutate(
    model = "Cyanobacteria abundance ~ stability state",
    term = recode(
      term,
      `(Intercept)` = "Mixed/Weak",
      `strat_stateModerate` = "Moderate vs Mixed/Weak",
      `strat_stateStrong` = "Strong vs Mixed/Weak"
    )
  ) |>
  select(model, term, estimate, conf.low, conf.high, p.value)

tbl_cyano_stability_clean

# ------------------------------------------------------------
# MODEL 2:
# Do deep nutrients vary by stability state?
# ThermalStability → Nutrients
# ------------------------------------------------------------

nutrient_stability_models <- bottom_nutrients |>
  group_by(nutrient) |>
  nest() |>
  mutate(
    model = map(
      data,
      ~ lm(log1p(concentration) ~ strat_state, data = .x)
    ),
    results = map(model, ~ tidy(.x, conf.int = TRUE))
  ) |>
  select(nutrient, results) |>
  unnest(results)

nutrient_stability_models

# Clean nutrient model table

nutrient_stability_results_clean <- nutrient_stability_models |>
  mutate(
    model = paste("Deep", nutrient, "~ stability state"),
    term = recode(
      term,
      `(Intercept)` = "Mixed/Weak",
      `strat_stateModerate` = "Moderate vs Mixed/Weak",
      `strat_stateStrong` = "Strong vs Mixed/Weak"
    )
  ) |>
  select(model, nutrient, term, estimate, conf.low, conf.high, p.value)

nutrient_stability_results_clean
# ------------------------------------------------------------
# MODEL 3:
# Does microcystin vary with cyanobacteria abundance and stability?
# CyanobacteriaAbundance + ThermalStability → Microcystin
# ------------------------------------------------------------

# First make a modeling dataframe
mc_cyano_model_data <- grab_toxins_brooks |>
  filter(
    toxin_class == "microcystin",
    depth_category %in% c("surface", "bottom")
  ) |>
  group_by(date, strat_state, stability_daily) |>
  summarize(
    microcystin = mean(total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    brooks_cyano_summary,
    by = "date"
  ) |>
  filter(
    !is.na(microcystin),
    !is.na(cyano_cells_l),
    !is.na(strat_state)
  )

model_mc_cyano_stability <- lm(
  log1p(microcystin) ~ log1p(cyano_cells_l) + strat_state,
  data = mc_cyano_model_data
)

summary(model_mc_cyano_stability)

tbl_mc_cyano_stability <- tidy(
  model_mc_cyano_stability,
  conf.int = TRUE
)

tbl_mc_cyano_stability
# Clean microcystin model table

tbl_mc_cyano_stability_clean <- tbl_mc_cyano_stability |>
  mutate(
    model = "Microcystin ~ cyanobacteria abundance + stability state",
    term = recode(
      term,
      `(Intercept)` = "Mixed/Weak",
      `log1p(cyano_cells_l)` = "Cyanobacteria abundance",
      `strat_stateModerate` = "Moderate vs Mixed/Weak",
      `strat_stateStrong` = "Strong vs Mixed/Weak"
    )
  ) |>
  select(model, term, estimate, conf.low, conf.high, p.value)

tbl_mc_cyano_stability_clean

final_model_table <- bind_rows(
  tbl_cyano_stability_clean,
  nutrient_stability_results_clean,
  tbl_mc_cyano_stability_clean
)

final_model_table


ggplot(
  mc_cyano_model_data,
  aes(
    x = cyano_cells_l,
    y = microcystin,
    color = strat_state
  )
) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  theme_classic()


### overlay it all 
# Surface overlay data ----------------------------------------------------

surface_nutrients <- nutrients_long |>
  filter(depth_zone == "Surface") |>
  transmute(
    date,
    parameter = nutrient,
    value = concentration
  )

surface_mc <- grab_toxins_brooks |>
  filter(toxin_class == "microcystin",
         depth_category == "surface") |>
  group_by(date) |>
  summarise(value = mean(total, na.rm = TRUE), .groups = "drop") |>
  mutate(parameter = "Surface microcystin")

surface_cyano <- brooks_cyano_summary |>
  transmute(
    date = as.Date(date),
    parameter = "Cyanobacteria abundance",
    value = cyano_cells_l
  )

stability_overlay <- schmidt_daily |>
  transmute(
    date,
    parameter = "Schmidt stability",
    value = stability_daily
  )

surface_overlay <- bind_rows(
  surface_nutrients,
  surface_mc,
  surface_cyano,
  stability_overlay
) |>
  filter(!is.na(value)) |>
  group_by(parameter) |>
  mutate(value_z = as.numeric(scale(value))) |>
  ungroup()

ggplot(surface_overlay,
       aes(x = date, y = value_z, color = parameter)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Date",
    y = "Standardized value (z-score)",
    color = "Parameter",
    title = "Brooks Lake surface samples: standardized overlay"
  ) +
  theme_classic(base_family = "Helvetica")
## bottom 
bottom_nutrients_overlay <- nutrients_long |>
  filter(depth_zone == "Deep") |>
  transmute(
    date,
    parameter = nutrient,
    value = concentration
  )

bottom_mc <- grab_toxins_brooks |>
  filter(toxin_class == "microcystin",
         depth_category == "bottom") |>
  group_by(date) |>
  summarise(value = mean(total, na.rm = TRUE), .groups = "drop") |>
  mutate(parameter = "Bottom microcystin")

bottom_do <- deep_do_daily |>
  transmute(
    date,
    parameter = "Bottom DO",
    value = deep_do_mgl
  )

bottom_overlay <- bind_rows(
  bottom_nutrients_overlay,
  bottom_mc,
  bottom_do,
  stability_overlay
) |>
  filter(!is.na(value)) |>
  group_by(parameter) |>
  mutate(value_z = as.numeric(scale(value))) |>
  ungroup()

ggplot(bottom_overlay,
       aes(x = date, y = value_z, color = parameter)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray70") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Date",
    y = "Standardized value (z-score)",
    color = "Parameter",
    title = "Brooks Lake bottom samples: standardized overlay"
  ) +
  theme_classic(base_family = "Helvetica")

colnames(deq_nutrients_clean_2025)
colnames(LakeData_Brooks2009_2024_Rstudio)
