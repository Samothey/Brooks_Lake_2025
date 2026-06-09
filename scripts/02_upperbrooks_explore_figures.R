# ======================================================================
# 02_upperbrooks_explore_patterns.R
# ======================================================================

library(tidyverse)
library(lubridate)
library(scales)

# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis/upperbrooks"
data_out_dir <- "data_clean/analysis/upperbrooks"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load data
# ======================================================================

upperbrooks_daily_drivers <- readRDS(
  "data_clean/analysis/upperbrooks/upperbrooks_daily_drivers_2025.rds"
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

# ======================================================================
# 3. Physical drivers timeline
# ======================================================================

drivers_long_ub <- upperbrooks_daily_drivers %>%
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

p_ub_physical_drivers <- ggplot(
  drivers_long_ub,
  aes(date, scaled_value, color = variable)
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  labs(
    title = "Upper Brooks physical drivers",
    x = NULL,
    y = "Scaled seasonal value",
    color = NULL
  ) +
  theme_bw()

p_ub_physical_drivers

ggsave(
  file.path(fig_dir, "upperbrooks_physical_drivers_scaled.png"),
  p_ub_physical_drivers,
  width = 10,
  height = 5,
  dpi = 300
)

# ======================================================================
# 4. Upper Brooks nutrients
# ======================================================================

nutrients_ub <- deq_nutrients_clean_2025 %>%
  mutate(
    date = as.Date(date),
    type = str_to_lower(type)
  ) %>%
  filter(
    lake == "Upper Brooks Lake",
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

surface_biology_ub <- nutrients_ub %>%
  filter(type == "surface") %>%
  select(
    date,
    chla,
    secchi
  )

# ======================================================================
# 5. Upper Brooks cyanobacteria
# ======================================================================

cyano_ub <- phyto_clean %>%
  mutate(
    date = as.Date(date),
    division = str_to_lower(division)
  ) %>%
  filter(
    lake == "Upper Brooks Lake",
    division %in% c("cyanophyta", "cyanobacteria"),
    sample_type != "duplicate"
  ) %>%
  group_by(date) %>%
  summarise(
    cyano_cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 6. Upper Brooks toxins: surface and depth separate
# ======================================================================

tox_ub_depth <- grab_tox %>%
  mutate(
    date = as.Date(sample_date),
    lake = str_to_lower(lake),
    site_type = str_trim(str_to_lower(site_type))
  ) %>%
  filter(
    lake == "upper brooks",
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
# 8. Build Upper Brooks story tables
# ======================================================================

story_surface_ub_base <- nutrients_ub %>%
  filter(type == "surface") %>%
  left_join(
    cyano_ub,
    by = "date"
  ) %>%
  left_join(
    upperbrooks_daily_drivers,
    by = "date"
  )

story_bottom_ub_base <- nutrients_ub %>%
  filter(type == "bottom") %>%
  select(-chla, -secchi) %>%
  left_join(
    surface_biology_ub,
    by = "date"
  ) %>%
  left_join(
    cyano_ub,
    by = "date"
  ) %>%
  left_join(
    upperbrooks_daily_drivers,
    by = "date"
  )

story_surface_ub <- join_nearest_toxin(
  story_surface_ub_base,
  tox_ub_depth %>% filter(type == "surface"),
  max_days = 2
)

story_bottom_ub <- join_nearest_toxin(
  story_bottom_ub_base,
  tox_ub_depth %>% filter(type == "bottom"),
  max_days = 2
)

story_surface_ub %>% count(toxin_day_diff)
story_bottom_ub %>% count(toxin_day_diff)

saveRDS(
  story_surface_ub,
  file.path(data_out_dir, "story_surface_upperbrooks_2025.rds")
)

saveRDS(
  story_bottom_ub,
  file.path(data_out_dir, "story_bottom_upperbrooks_2025.rds")
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
    scale_color_manual(
      values = c(
        "Schmidt stability" = "black",
        "Lower deep DO" = "#a6cee3",
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

# ======================================================================
# 10. Make Upper Brooks story plots
# ======================================================================

p_ub_surface_story <- make_scaled_story_plot(
  story_surface_ub,
  "Upper Brooks surface seasonal trends"
)

p_ub_bottom_story <- make_scaled_story_plot(
  story_bottom_ub,
  "Upper Brooks bottom seasonal trends"
)

p_ub_surface_story
p_ub_bottom_story

ggsave(
  file.path(fig_dir, "upperbrooks_surface_scaled_story.png"),
  p_ub_surface_story,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "upperbrooks_bottom_scaled_story.png"),
  p_ub_bottom_story,
  width = 10,
  height = 5,
  dpi = 300
)

# ======================================================================
# 11. Raw surface vs bottom nutrients
# ======================================================================

nutrients_ub_long <- nutrients_ub %>%
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

p_ub_nutrients_raw <- ggplot(
  nutrients_ub_long,
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
    title = "Upper Brooks surface and bottom nutrients",
    x = NULL,
    y = "Concentration",
    color = NULL
  ) +
  theme_bw()

p_ub_nutrients_raw

ggsave(
  file.path(fig_dir, "upperbrooks_raw_surface_bottom_nutrients.png"),
  p_ub_nutrients_raw,
  width = 9,
  height = 7,
  dpi = 300
)

# ======================================================================
# 12. Raw surface vs bottom toxins
# ======================================================================

p_ub_tox_raw <- ggplot(
  tox_ub_depth,
  aes(date, total_mc, color = type, group = type)
) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  labs(
    title = "Upper Brooks buoy toxins",
    x = NULL,
    y = "Total microcystins",
    color = NULL
  ) +
  theme_bw()

p_ub_tox_raw

ggsave(
  file.path(fig_dir, "upperbrooks_raw_surface_bottom_toxins.png"),
  p_ub_tox_raw,
  width = 9,
  height = 4,
  dpi = 300
)