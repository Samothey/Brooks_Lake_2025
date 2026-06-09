# ======================================================================
# 02_lowerjade_explore_patterns.R
#
# Purpose:
#   Explore Lower Jade seasonal patterns:
#   - physical drivers
#   - nutrients
#   - Chl-a / Secchi
#   - cyanobacteria
#   - toxins
# ======================================================================

library(tidyverse)
library(lubridate)
library(scales)

# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis/lowerjade"
data_out_dir <- "data_clean/analysis/lowerjade"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load data
# ======================================================================

lowerjade_daily_drivers <- readRDS(
  "data_clean/analysis/lowerjade/lowerjade_daily_drivers_2025.rds"
)

deq_nutrients_clean_2025 <- readRDS(
  "data_clean/deq/deq_nutrients_clean_2025.rds"
)

phyto_clean <- readRDS(
  "data_clean/phytoplankton/phyto_clean.rds"
)

grab_tox <- readRDS(
  "data_clean/toxins/grab_tox.rds"
)



schmidt_hourly_lj <- readRDS(
  "data_clean/analysis/lowerjade/schmidt_hourly_lowerjade_2025.rds"
)

schmidt_daily_lj <- readRDS(
  "data_clean/analysis/lowerjade/schmidt_daily_lowerjade_2025.rds"
)

deep_do_daily_lj <- readRDS(
  "data_clean/analysis/lowerjade/deep_do_daily_lowerjade_2025.rds"
)

lowerjade_daily_drivers <- readRDS(
  "data_clean/analysis/lowerjade/lowerjade_daily_drivers_2025.rds"
)

# ----------------------------------------------------------------------
# 3. Sampling dates
# ----------------------------------------------------------------------

sampling_dates_lj <- as.POSIXct(
  paste0(
    c(
      "2025-06-23",
      "2025-07-08",
      "2025-07-22",
      "2025-08-05",
      "2025-08-19",
      "2025-09-02",
      "2025-09-16",
      "2025-09-30",
      "2025-10-14"
    ),
    " 12:00:00"
  ),
  tz = "America/Denver"
)

# ----------------------------------------------------------------------
# 4. Hourly stability plot
# ----------------------------------------------------------------------

lowerjade_stability_plot <- ggplot(
  schmidt_hourly_lj,
  aes(x = datetime_hour)
) +
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
    aes(y = stability_24hr, color = strat_state_hourly),
    size = 1.2,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = sampling_dates_lj,
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
    title = "Lower Jade Schmidt stability"
  ) +
  theme_classic(base_family = "Helvetica")

lowerjade_stability_plot

ggsave(
  filename = file.path(fig_dir, "lowerjade_schmidt_stability_hourly.png"),
  plot = lowerjade_stability_plot,
  width = 10,
  height = 4.5,
  dpi = 300
)

# ----------------------------------------------------------------------
# 5. Daily stability and deep DO
# ----------------------------------------------------------------------

lowerjade_stability_do_daily <- lowerjade_daily_drivers %>%
  mutate(date = as.Date(date)) %>%
  select(
    date,
    stability_daily,
    stability_min,
    stability_max,
    strat_state,
    deep_do_mgl
  )

p_lj_daily_stability_do <- ggplot(
  lowerjade_stability_do_daily,
  aes(x = date)
) +
  geom_line(
    aes(y = stability_daily),
    linewidth = 1,
    color = "black",
    na.rm = TRUE
  ) +
  geom_point(
    aes(y = stability_daily, color = strat_state),
    size = 2,
    na.rm = TRUE
  ) +
  geom_line(
    aes(
      y = deep_do_mgl /
        max(deep_do_mgl, na.rm = TRUE) *
        max(stability_daily, na.rm = TRUE)
    ),
    linewidth = 1,
    linetype = "dashed",
    color = "firebrick",
    na.rm = TRUE
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
    x = NULL,
    y = "Schmidt stability; deep DO scaled",
    color = "Stability state",
    title = "Lower Jade daily stability and deep DO"
  ) +
  theme_bw()

p_lj_daily_stability_do

ggsave(
  filename = file.path(fig_dir, "lowerjade_daily_stability_deep_do.png"),
  plot = p_lj_daily_stability_do,
  width = 10,
  height = 4.5,
  dpi = 300
)

# ----------------------------------------------------------------------
# 6. Daily deep DO by stratification state
# ----------------------------------------------------------------------

deep_do_strat_daily_lj <- lowerjade_daily_drivers %>%
  mutate(date = as.Date(date)) %>%
  filter(
    !is.na(strat_state),
    !is.na(deep_do_mgl)
  )

table(deep_do_strat_daily_lj$strat_state)

deep_do_boxplot_daily_lj <- ggplot(
  deep_do_strat_daily_lj,
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
      "Moderate" = "#fdae61",
      "Strong" = "#d7191c"
    ),
    na.translate = FALSE
  ) +
  labs(
    x = "Stratification state",
    y = "Daily mean deep DO (mg/L)",
    title = "Lower Jade deep DO by stratification state"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

deep_do_boxplot_daily_lj

ggsave(
  filename = file.path(fig_dir, "lowerjade_deep_do_by_strat_state.png"),
  plot = deep_do_boxplot_daily_lj,
  width = 6,
  height = 4.5,
  dpi = 300
)



# ======================================================================
# 3. Physical drivers timeline
# ======================================================================

drivers_long_lj <- lowerjade_daily_drivers %>%
  select(
    date,
    stability_daily,
    deep_do_mgl,
    air_temp_mean_c,
    wind_speed_mean_ms,
    gust_speed_max_ms,
    snotel_prcp_mm
  ) %>%
  pivot_longer(
    cols = -date,
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  group_by(variable) %>%
  mutate(
    scaled_value = value / max(value, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    variable = recode(
      variable,
      stability_daily = "Schmidt stability",
      deep_do_mgl = "Deep DO",
      air_temp_mean_c = "Air temperature",
      wind_speed_mean_ms = "Mean wind speed",
      gust_speed_max_ms = "Max gust speed",
      snotel_prcp_mm = "Precipitation"
    )
  )

p_lj_physical_drivers <- ggplot(
  drivers_long_lj,
  aes(date, scaled_value, color = variable)
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  labs(
    title = "Lower Jade physical drivers",
    x = NULL,
    y = "Scaled seasonal value",
    color = NULL
  ) +
  theme_bw()

p_lj_physical_drivers

ggsave(
  file.path(fig_dir, "lowerjade_physical_drivers_scaled.png"),
  p_lj_physical_drivers,
  width = 10,
  height = 5,
  dpi = 300
)

# ======================================================================
# 4. Lower Jade nutrients
# ======================================================================

nutrients_lj <- deq_nutrients_clean_2025 %>%
  mutate(
    date = as.Date(date),
    type = str_to_lower(type)
  ) %>%
  filter(
    lake == "Lower Jade Lake",
    type %in% c("surface", "bottom")
  ) %>%
  select(
    date,
    type,
    ammonia,
    tn,
    tp,
    chla,
    secchi
  ) %>%
  arrange(date, type)

surface_biology_lj <- nutrients_lj %>%
  filter(type == "surface") %>%
  select(
    date,
    chla,
    secchi
  )

# ======================================================================
# 5. Lower Jade cyanobacteria
# ======================================================================

cyano_lj <- phyto_clean %>%
  mutate(
    date = as.Date(date),
    division = str_to_lower(division)
  ) %>%
  filter(
    lake == "Lower Jade Lake",
    division %in% c("cyanophyta", "cyanobacteria"),
    sample_type != "duplicate"
  ) %>%
  group_by(date) %>%
  summarise(
    cyano_cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 6. Lower Jade toxins: surface and depth separate
# ======================================================================

tox_lj_depth <- grab_tox %>%
  mutate(
    date = as.Date(sample_date),
    lake = str_to_lower(lake),
    site_type = str_trim(str_to_lower(site_type))
  ) %>%
  filter(
    lake == "lower jade",
    site_type %in% c("buoy_surface", "buoy_depth"),
    sample_type != "duplicate"
  ) %>%
  mutate(
    type = case_when(
      site_type == "buoy_surface" ~ "surface",
      site_type == "buoy_depth" ~ "bottom",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(type)) %>%
  group_by(date, type) %>%
  summarise(
    total_mc = mean(total_mc, na.rm = TRUE),
    max_total_mc = max(total_mc, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 7. Nearest toxin join within +/- 2 days
# ======================================================================

join_nearest_toxin <- function(story_data, toxin_data, max_days = 2) {
  
  story_data %>%
    mutate(row_id = row_number()) %>%
    left_join(
      toxin_data,
      by = "type",
      relationship = "many-to-many"
    ) %>%
    mutate(
      day_diff = abs(as.numeric(date.x - date.y))
    ) %>%
    filter(
      is.na(date.y) | day_diff <= max_days
    ) %>%
    group_by(row_id) %>%
    slice_min(day_diff, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      date = date.x,
      type,
      ammonia,
      tn,
      tp,
      chla,
      secchi,
      cyano_cells,
      stability_daily,
      stability_min,
      stability_max,
      strat_state,
      deep_do_mgl,
      total_mc,
      max_total_mc,
      toxin_sample_date = date.y,
      toxin_day_diff = day_diff
    )
}

# ======================================================================
# 8. Build Lower Jade story tables
# ======================================================================

story_surface_lj_base <- nutrients_lj %>%
  filter(type == "surface") %>%
  left_join(
    cyano_lj,
    by = "date"
  ) %>%
  left_join(
    lowerjade_daily_drivers,
    by = "date"
  )

story_bottom_lj_base <- nutrients_lj %>%
  filter(type == "bottom") %>%
  select(-chla, -secchi) %>%
  left_join(
    surface_biology_lj,
    by = "date"
  ) %>%
  left_join(
    cyano_lj,
    by = "date"
  ) %>%
  left_join(
    lowerjade_daily_drivers,
    by = "date"
  )

story_surface_lj <- join_nearest_toxin(
  story_surface_lj_base,
  tox_lj_depth %>% filter(type == "surface"),
  max_days = 2
)

story_bottom_lj <- join_nearest_toxin(
  story_bottom_lj_base,
  tox_lj_depth %>% filter(type == "bottom"),
  max_days = 2
)

story_surface_lj %>% count(toxin_day_diff)
story_bottom_lj %>% count(toxin_day_diff)

saveRDS(
  story_surface_lj,
  file.path(data_out_dir, "story_surface_lowerjade_2025.rds")
)

saveRDS(
  story_bottom_lj,
  file.path(data_out_dir, "story_bottom_lowerjade_2025.rds")
)

# ======================================================================
# 9. Scaled story plot function
# ======================================================================

make_scaled_story_plot <- function(data, title_text) {
  
  plot_long <- data %>%
    select(
      date,
      stability_daily,
      deep_do_mgl,
      ammonia,
      tn,
      tp,
      chla,
      secchi,
      cyano_cells,
      total_mc
    ) %>%
    pivot_longer(
      cols = -date,
      names_to = "variable",
      values_to = "value"
    ) %>%
    filter(!is.na(value)) %>%
    group_by(variable) %>%
    mutate(
      value_plot = case_when(
        variable == "secchi" ~ max(value, na.rm = TRUE) - value,
        variable == "deep_do_mgl" ~ max(value, na.rm = TRUE) - value,
        TRUE ~ value
      ),
      scaled_value = value_plot / max(value_plot, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      variable = recode(
        variable,
        stability_daily = "Schmidt stability",
        deep_do_mgl = "Lower deep DO",
        ammonia = "Ammonia",
        tn = "Total nitrogen",
        tp = "Total phosphorus",
        chla = "Chl-a",
        secchi = "Lower clarity",
        cyano_cells = "Cyanobacteria cells",
        total_mc = "Total microcystins"
      )
    )
  
  ggplot(
    plot_long,
    aes(
      x = date,
      y = scaled_value,
      color = variable,
      group = variable
    )
  ) +
    geom_line(linewidth = 1.1, alpha = 0.9, na.rm = TRUE) +
    geom_point(size = 2, alpha = 0.9, na.rm = TRUE) +
    labs(
      title = title_text,
      x = NULL,
      y = "Scaled seasonal value",
      color = NULL
    ) +
    theme_bw()
}

# ======================================================================
# 10. Make Lower Jade story plots
# ======================================================================

p_lj_surface_story <- make_scaled_story_plot(
  story_surface_lj,
  "Lower Jade surface seasonal trends"
)

p_lj_bottom_story <- make_scaled_story_plot(
  story_bottom_lj,
  "Lower Jade bottom seasonal trends"
)

p_lj_surface_story
p_lj_bottom_story

ggsave(
  file.path(fig_dir, "lowerjade_surface_scaled_story.png"),
  p_lj_surface_story,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "lowerjade_bottom_scaled_story.png"),
  p_lj_bottom_story,
  width = 10,
  height = 5,
  dpi = 300
)

# ======================================================================
# 11. Raw surface vs bottom nutrients
# ======================================================================

nutrients_lj_long <- nutrients_lj %>%
  pivot_longer(
    cols = c(ammonia, tn, tp),
    names_to = "nutrient",
    values_to = "concentration"
  ) %>%
  mutate(
    nutrient = recode(
      nutrient,
      ammonia = "Ammonia",
      tn = "Total nitrogen",
      tp = "Total phosphorus"
    )
  )

p_lj_nutrients_raw <- ggplot(
  nutrients_lj_long,
  aes(date, concentration, color = type, group = type)
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(
    ~ nutrient,
    ncol = 1,
    scales = "free_y"
  ) +
  labs(
    title = "Lower Jade surface and bottom nutrients",
    x = NULL,
    y = "Concentration",
    color = NULL
  ) +
  theme_bw()

p_lj_nutrients_raw

ggsave(
  file.path(fig_dir, "lowerjade_raw_surface_bottom_nutrients.png"),
  p_lj_nutrients_raw,
  width = 9,
  height = 7,
  dpi = 300
)

# ======================================================================
# 12. Raw surface vs bottom toxins
# ======================================================================

p_lj_tox_raw <- ggplot(
  tox_lj_depth,
  aes(date, total_mc, color = type, group = type)
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  labs(
    title = "Lower Jade buoy toxins",
    x = NULL,
    y = "Total microcystins",
    color = NULL
  ) +
  theme_bw()

p_lj_tox_raw

ggsave(
  file.path(fig_dir, "lowerjade_raw_surface_bottom_toxins.png"),
  p_lj_tox_raw,
  width = 9,
  height = 4,
  dpi = 300
)