# ======================================================================
# 02_rainbow_explore_story.R
#
# Purpose:
#   Build Rainbow Lake exploratory story plots using:
#   - profile-based physical drivers
#   - surface/bottom nutrients
#   - Chl-a and Secchi
#   - cyanobacteria cells
#   - surface/depth grab toxins
#
# Notes:
#   Rainbow uses profile-based stability from DEQ vertical profiles,
#   not daily buoy-derived stability.
# ======================================================================

library(tidyverse)
library(lubridate)
library(scales)

# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis/rainbow"
data_out_dir <- "data_clean/analysis"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load cleaned / built data
# ======================================================================

rainbow_profile_drivers <- readRDS(
  file.path(data_out_dir, "rainbow_profile_drivers_2025.rds")
)

grab_tox <- readRDS(
  "data_clean/toxins/grab_tox.rds"
)

phyto_clean <- readRDS(
  "data_clean/phytoplankton/phyto_clean.rds"
)

deq_nutrients_clean_2025 <- readRDS(
  "data_clean/deq/deq_nutrients_clean_2025.rds"
)


# ======================================================================
# Rainbow 7 m DO from buoy file
# ======================================================================

library(tidyverse)
library(lubridate)
library(janitor)

Rainbow_WQ_2025_cleaned <- read_csv(
  "data_clean/buoy/Rainbow_WQ_2025_cleaned.csv",
  show_col_types = FALSE
)

rainbow_deep_do_daily <- Rainbow_WQ_2025_cleaned %>%
  clean_names() %>%
  mutate(
    date = as.Date(datetime_mst)
  ) %>%
  group_by(date) %>%
  summarise(
    deep_do_mgl = mean(odomg_l_7m, na.rm = TRUE),
    deep_do_min_mgl = min(odomg_l_7m, na.rm = TRUE),
    deep_do_max_mgl = max(odomg_l_7m, na.rm = TRUE),
    deep_do_sat_pct = mean(odosat_7m, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(
  rainbow_deep_do_daily,
  "data_clean/analysis/rainbow_deep_do_daily_2025.rds"
)

rainbow_profile_drivers <- rainbow_profile_drivers %>%
  select(-any_of(c("bottom_do_mgl", "deep_do_mgl"))) %>%
  left_join(
    rainbow_deep_do_daily,
    by = "date"
  )
# ======================================================================
# 3. Physical driver timeline
# ======================================================================

physical_long_rainbow <- rainbow_profile_drivers %>%
  select(
    date,
    stability_profile,
    temp_diff,
    do_diff,
    bottom_do_mgl,
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
    value_plot = case_when(
      variable == "bottom_do_mgl" ~ max(value, na.rm = TRUE) - value,
      TRUE ~ value
    ),
    scaled_value = value_plot / max(value_plot, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    variable = recode(
      variable,
      stability_profile = "Schmidt stability",
      temp_diff = "Temp difference",
      do_diff = "DO difference",
      bottom_do_mgl = "Low bottom DO",
      air_temp_mean_c = "Air temperature",
      wind_speed_mean_ms = "Mean wind speed",
      gust_speed_max_ms = "Max gust speed",
      snotel_prcp_mm = "Precipitation"
    )
  )

p_rainbow_physical_scaled <- ggplot(
  physical_long_rainbow,
  aes(date, scaled_value, color = variable)
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  labs(
    x = NULL,
    y = "Scaled seasonal value",
    color = NULL,
    title = "Rainbow Lake physical drivers"
  ) +
  theme_bw()

p_rainbow_physical_scaled

ggsave(
  file.path(fig_dir, "rainbow_physical_drivers_scaled_timeline.png"),
  p_rainbow_physical_scaled,
  width = 10,
  height = 5,
  dpi = 300
)

# ======================================================================
# 4. Build Rainbow nutrients, cyano, and toxin tables
# ======================================================================

nutrients_rainbow <- deq_nutrients_clean_2025 %>%
  mutate(
    date = as.Date(date),
    type = str_to_lower(type)
  ) %>%
  filter(
    lake == "Rainbow Lake",
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

surface_biology_rainbow <- nutrients_rainbow %>%
  filter(type == "surface") %>%
  select(
    date,
    chla,
    secchi
  )

cyano_density_rainbow <- phyto_clean %>%
  mutate(
    date = as.Date(date),
    division = str_to_lower(division)
  ) %>%
  filter(
    lake == "Rainbow Lake",
    division %in% c("cyanophyta", "cyanobacteria"),
    sample_type != "duplicate"
  ) %>%
  group_by(date) %>%
  summarise(
    cyano_cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  )

tox_rainbow_depth <- grab_tox %>%
  mutate(
    date = as.Date(sample_date),
    lake = str_to_lower(lake),
    site_type = str_trim(str_to_lower(site_type))
  ) %>%
  filter(
    lake %in% c("rainbow", "rainbow lake"),
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

stability_rainbow_plot <- rainbow_profile_drivers %>%
  mutate(date = as.Date(date)) %>%
  select(
    date,
    stability_profile,
    temp_diff,
    do_diff,
    bottom_do_mgl,
    strat_state
  )

# ======================================================================
# 5. Nearest toxin join within +/- 2 days
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
      stability_profile,
      temp_diff,
      do_diff,
      bottom_do_mgl,
      strat_state,
      total_mc,
      max_total_mc,
      toxin_sample_date = date.y,
      toxin_day_diff = day_diff
    )
}

# ======================================================================
# 6. Build Rainbow surface and bottom story tables
# ======================================================================

story_surface_rainbow_base <- nutrients_rainbow %>%
  filter(type == "surface") %>%
  left_join(
    cyano_density_rainbow,
    by = "date"
  ) %>%
  left_join(
    stability_rainbow_plot,
    by = "date"
  )

story_bottom_rainbow_base <- nutrients_rainbow %>%
  filter(type == "bottom") %>%
  select(-chla, -secchi) %>%
  left_join(
    surface_biology_rainbow,
    by = "date"
  ) %>%
  left_join(
    cyano_density_rainbow,
    by = "date"
  ) %>%
  left_join(
    stability_rainbow_plot,
    by = "date"
  )

story_surface_rainbow <- join_nearest_toxin(
  story_surface_rainbow_base,
  tox_rainbow_depth %>% filter(type == "surface"),
  max_days = 2
)

story_bottom_rainbow <- join_nearest_toxin(
  story_bottom_rainbow_base,
  tox_rainbow_depth %>% filter(type == "bottom"),
  max_days = 2
)

story_surface_rainbow %>% count(toxin_day_diff)
story_bottom_rainbow %>% count(toxin_day_diff)

saveRDS(
  story_surface_rainbow,
  file.path(data_out_dir, "rainbow_story_surface_2025.rds")
)

saveRDS(
  story_bottom_rainbow,
  file.path(data_out_dir, "rainbow_story_bottom_2025.rds")
)

# ======================================================================
# 7. Scaled surface and bottom story plots
# ======================================================================

make_scaled_story_plot_rainbow <- function(data, title_text) {
  
  plot_long <- data %>%
    select(
      date,
      stability_profile,
      temp_diff,
      do_diff,
      bottom_do_mgl,
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
        variable == "bottom_do_mgl" ~ max(value, na.rm = TRUE) - value,
        TRUE ~ value
      ),
      scaled_value = value_plot / max(value_plot, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      variable = recode(
        variable,
        stability_profile = "Schmidt stability",
        temp_diff = "Temp difference",
        do_diff = "DO difference",
        bottom_do_mgl = "Low bottom DO",
        ammonia = "Ammonia",
        tn = "Total nitrogen",
        tp = "Total phosphorus",
        chla = "Chl-a",
        secchi = "Lower clarity",
        cyano_cells = "Cyanobacteria cells",
        total_mc = "Total microcystins"
      ),
      line_type = case_when(
        variable %in% c(
          "Schmidt stability",
          "Temp difference",
          "DO difference",
          "Low bottom DO"
        ) ~ "solid",
        variable %in% c(
          "Ammonia",
          "Total nitrogen",
          "Total phosphorus"
        ) ~ "dashed",
        TRUE ~ "dotted"
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
    geom_line(
      aes(linetype = line_type),
      linewidth = 1.1,
      alpha = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      size = 2,
      alpha = 0.9,
      na.rm = TRUE
    ) +
    scale_linetype_identity() +
    scale_color_manual(
      values = c(
        "Schmidt stability" = "black",
        "Temp difference" = "#984ea3",
        "DO difference" = "#377eb8",
        "Low bottom DO" = "#a65628",
        "Ammonia" = "#1b9e77",
        "Total nitrogen" = "#66a61e",
        "Total phosphorus" = "#d95f02",
        "Chl-a" = "#7570b3",
        "Lower clarity" = "#e7298a",
        "Cyanobacteria cells" = "#1f78b4",
        "Total microcystins" = "#e31a1c"
      )
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = "Scaled seasonal value",
      color = NULL
    ) +
    theme_bw()
}

p_rainbow_surface_story <- make_scaled_story_plot_rainbow(
  story_surface_rainbow,
  "Rainbow surface seasonal trends"
)

p_rainbow_bottom_story <- make_scaled_story_plot_rainbow(
  story_bottom_rainbow,
  "Rainbow bottom seasonal trends"
)

p_rainbow_surface_story
p_rainbow_bottom_story

ggsave(
  file.path(fig_dir, "rainbow_surface_scaled_story_plot.png"),
  p_rainbow_surface_story,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "rainbow_bottom_scaled_story_plot.png"),
  p_rainbow_bottom_story,
  width = 10,
  height = 5,
  dpi = 300
)

# ======================================================================
# 8. Raw surface vs bottom nutrient plots
# ======================================================================

nutrients_rainbow_long <- nutrients_rainbow %>%
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

p_rainbow_nutrients_raw <- ggplot(
  nutrients_rainbow_long,
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
    x = NULL,
    y = "Concentration",
    color = NULL,
    title = "Rainbow Lake surface and bottom nutrients"
  ) +
  theme_bw()

p_rainbow_nutrients_raw

ggsave(
  file.path(fig_dir, "rainbow_raw_surface_bottom_nutrients.png"),
  p_rainbow_nutrients_raw,
  width = 9,
  height = 7,
  dpi = 300
)

# ======================================================================
# 9. Raw surface vs bottom toxins
# ======================================================================

p_rainbow_tox_raw <- ggplot(
  tox_rainbow_depth,
  aes(date, total_mc, color = type, group = type)
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  labs(
    x = NULL,
    y = "Total microcystins",
    color = NULL,
    title = "Rainbow Lake buoy toxins"
  ) +
  theme_bw()

p_rainbow_tox_raw

ggsave(
  file.path(fig_dir, "rainbow_raw_surface_bottom_toxins.png"),
  p_rainbow_tox_raw,
  width = 9,
  height = 4,
  dpi = 300
)

# ======================================================================
# 10. Stability / profile-driver relationship plots
# ======================================================================

p_rainbow_stab_temp <- ggplot(
  rainbow_profile_drivers,
  aes(temp_diff, stability_profile)
) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Surface-bottom temperature difference (°C)",
    y = "Schmidt stability"
  ) +
  theme_bw()

p_rainbow_stab_do <- ggplot(
  rainbow_profile_drivers,
  aes(stability_profile, do_diff)
) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Schmidt stability",
    y = "Surface-bottom DO difference (mg/L)"
  ) +
  theme_bw()

p_rainbow_bottom_do <- ggplot(
  rainbow_profile_drivers,
  aes(stability_profile, bottom_do_mgl)
) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Schmidt stability",
    y = "Bottom DO (mg/L)"
  ) +
  theme_bw()

p_rainbow_stab_temp
p_rainbow_stab_do
p_rainbow_bottom_do