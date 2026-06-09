# ======================================================================
# Brooks Lake integrated analysis
#
# Purpose:
#   1. Load cleaned Brooks datasets
#   2. Recalculate Schmidt stability from buoy temperature profiles
#   3. Create daily stability and DO summaries
#   4. Join weather, toxins, phyto, nutrients, and stability
# ======================================================================


# ======================================================================
# 0. Libraries
# ======================================================================

library(tidyverse)
library(lubridate)
library(janitor)
library(zoo)
library(rLakeAnalyzer)
library(readr)
library(readxl)
library(scales)
<<<<<<< HEAD
=======
library(plotly)
>>>>>>> fffb24a (:))


# ======================================================================
# 1. File paths
# ======================================================================

fig_dir <- "figures/integrated_analysis"
data_out_dir <- "data_clean/analysis"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)


# ======================================================================
# 2. Load cleaned datasets
# ======================================================================

grab_tox <- readRDS(
  "data_clean/toxins/grab_tox.rds"
)

<<<<<<< HEAD
Spatt_tox <- readRDS(
  "data_clean/toxins/Spatt_tox.rds"
)

tox_plot <- readRDS(
  "data_clean/toxins/tox_plot.rds"
)

=======

brooks_daily_drivers <- readRDS("~/Desktop/Project/Brooks_lake_2025/data_clean/analysis")
>>>>>>> fffb24a (:))

brooks_weather_noaa_daily_2025 <- readRDS(
  "data_clean/weather/brooks_weather_daily_2025.rds"
)

phyto_clean <- readRDS(
  "data_clean/phytoplankton/phyto_clean.rds"
)

deq_nutrients_clean_2025 <- readRDS(
  "data_clean/deq/deq_nutrients_clean_2025.rds"
)

brooks_wq_2025_cleaned <- read_csv(
  "data_clean/buoy/brooks_wq_2025_cleaned.csv",
  show_col_types = FALSE
)
# ======================================================================
# 3. Prepare Brooks buoy temperature profile data
# ======================================================================

brooks_buoy <- brooks_wq_2025_cleaned %>%
  select(-any_of("...1")) %>%
  rename(WTemp0.75m = WTempC) %>%
  mutate(
    datetime = DATETIME_MST,
    datetime = ymd_hms(datetime, tz = "America/Denver")
  )

wtemp <- brooks_buoy %>%
  select(datetime, matches("^WTemp")) %>%
  pivot_longer(
    cols = matches("^WTemp"),
    names_to = "depth",
    values_to = "temp_c"
  ) %>%
  mutate(
    depth = str_remove(depth, "WTemp"),
    depth = str_remove(depth, "m"),
    depth = as.numeric(depth)
  ) %>%
  filter(!is.na(temp_c))

wtemp_hourly <- wtemp %>%
  mutate(datetime_hour = round_date(datetime, unit = "hour")) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 4. Interpolate temperature profiles for heatmap
# ======================================================================

interp_hourly <- wtemp_hourly %>%
  group_by(datetime_hour) %>%
  nest() %>%
  mutate(n_depths = map_dbl(data, nrow)) %>%
  filter(n_depths > 1) %>%
  mutate(data = map(data, ~ arrange(.x, depth))) %>%
  mutate(
    interp_fun = map(data, ~ approxfun(.x$depth, .x$temp_c)),
    raster = map2(data, interp_fun, function(df, func) {
      tibble(
        depth = seq(min(df$depth), max(df$depth), by = 0.2),
        temp_c = func(depth)
      )
    })
  ) %>%
  select(datetime_hour, raster) %>%
  unnest(raster) %>%
  filter(!is.na(temp_c))

# ======================================================================
# 5. Calculate hourly Schmidt stability
# ======================================================================

brooks_temp_wide <- wtemp_hourly %>%
  mutate(
    depth_name = paste0("wtr_", depth)
  ) %>%
  select(datetime = datetime_hour, depth_name, temp_c) %>%
  pivot_wider(
    names_from = depth_name,
    values_from = temp_c
  )

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

brooks_temp_wide_full <- brooks_temp_wide %>%
  select(datetime, all_of(temp_cols)) %>%
  filter(if_all(all_of(temp_cols), ~ !is.na(.x)))

# Estimated Brooks bathymetry
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

schmidt_hourly <- brooks_temp_wide_full %>%
  rowwise() %>%
  mutate(
    schmidt_stability = schmidt.stability(
      wtr = c_across(all_of(temp_cols)),
      depths = temp_depths,
      bthD = bthD,
      bthA = bthA
    )
  ) %>%
  ungroup() %>%
  select(datetime_hour = datetime, schmidt_stability) %>%
  arrange(datetime_hour) %>%
  mutate(
    stability_24hr = as.numeric(
      rollmean(
        schmidt_stability,
        k = 24,
        fill = NA,
        align = "center"
      )
    ),
    strat_state_hourly = case_when(
      stability_24hr <= 25 ~ "Mixed/Weak",
      stability_24hr <= 50 ~ "Moderate",
      stability_24hr > 50 ~ "Strong",
      TRUE ~ NA_character_
    ),
    strat_state_hourly = factor(
      strat_state_hourly,
      levels = c("Mixed/Weak", "Moderate", "Strong")
    )
  )
# ======================================================================
# 6. Daily Schmidt stability
# ======================================================================

schmidt_daily <- schmidt_hourly %>%
  mutate(date = as.Date(datetime_hour)) %>%
  filter(!is.na(stability_24hr)) %>%
  group_by(date) %>%
  summarise(
    stability_daily = mean(stability_24hr, na.rm = TRUE),
    stability_max = max(stability_24hr, na.rm = TRUE),
    stability_min = min(stability_24hr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    strat_state = case_when(
      stability_daily <= 25 ~ "Mixed/Weak",
      stability_daily <= 50 ~ "Moderate",
      stability_daily > 50 ~ "Strong",
      TRUE ~ NA_character_
    ),
    strat_state = factor(
      strat_state,
      levels = c("Mixed/Weak", "Moderate", "Strong")
    )
  )

# ======================================================================
# 7. Daily dissolved oxygen
# ======================================================================

do_brooks <- brooks_buoy %>%
  select(datetime, odomgL, odomgL_15m) %>%
  pivot_longer(
    cols = c(odomgL, odomgL_15m),
    names_to = "depth",
    values_to = "do_mgl"
  ) %>%
  mutate(
    depth = case_when(
      depth == "odomgL" ~ "Surface",
      depth == "odomgL_15m" ~ "15 m",
      TRUE ~ depth
    )
  ) %>%
  filter(!is.na(do_mgl))

do_hourly <- do_brooks %>%
  mutate(datetime_hour = round_date(datetime, unit = "hour")) %>%
  group_by(datetime_hour, depth) %>%
  summarise(
    do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

deep_do_daily <- do_hourly %>%
  filter(depth == "15 m") %>%
  mutate(date = as.Date(datetime_hour)) %>%
  group_by(date) %>%
  summarise(
    deep_do_mgl = mean(do_mgl, na.rm = TRUE),
    .groups = "drop"
  )

# ======================================================================
# 8. Create Brooks daily driver table
# ======================================================================

<<<<<<< HEAD
brooks_daily_drivers <- brooks_weather_noaa_daily_2025 %>%
  mutate(date = as.Date(date)) %>%
  left_join(schmidt_daily, by = "date") %>%
  left_join(deep_do_daily, by = "date")

saveRDS(
  brooks_daily_drivers,
  file.path(data_out_dir, "brooks_daily_drivers_2025.rds"))

colnames(brooks_daily_drivers)
  
 ## start looking for trends with stability  
=======

# ----------------------------------------------------------------------
# Join daily drivers
# ----------------------------------------------------------------------

##brooks_daily_drivers <- brooks_weather_noaa_daily_2025 %>%
 # mutate(date = as.Date(date)) %>%
 # left_join(schmidt_daily, by = "date") %>%
 # left_join(deep_do_daily, by = "date")

##saveRDS(
 # brooks_daily_drivers,
 # file.path(data_out_dir, "brooks_daily_drivers_2025.rds")
#)

# ----------------------------------------------------------------------
# Timeline of drivers as interavtive 
# ----------------------------------------------------------------------

drivers_long <- brooks_daily_drivers %>%
  select(
    date,
    stability_daily,
    air_temp_mean_c,
    wind_speed_mean_ms,
    gust_speed_max_ms,
    snotel_prcp_mm,
    deep_do_mgl
  ) %>%
  pivot_longer(
    cols = -date,
    names_to = "variable",
    values_to = "value"
  )

p_drivers <- ggplot(
  drivers_long,
  aes(
    x = date,
    y = value,
    group = 1,
    text = paste(
      "Date:", date,
      "<br>Value:", round(value, 3)
    )
  )
) +
  geom_line(
    linewidth = 1,
    color = "black"
  ) +
  facet_wrap(
    ~ variable,
    ncol = 1,
    scales = "free_y"
  ) +
  theme_bw()

ggplotly(
  p_drivers,
  tooltip = "text"
)
p_drivers_interactive


  
 ## start looking for trends with stability 
## Timeline of drivers 
>>>>>>> fffb24a (:))
  drivers_long <- brooks_daily_drivers %>%
    select(
      date,
      stability_daily,
      air_temp_mean_c,
      wind_speed_mean_ms,
      gust_speed_max_ms,
      snotel_prcp_mm
    ) %>%
    pivot_longer(
      cols = -date,
      names_to = "variable",
      values_to = "value"
    )
  
  ggplot(
    drivers_long,
    aes(date, value)
  ) +
    geom_line() +
    facet_wrap(
      ~ variable,
      ncol = 1,
      scales = "free_y"
    ) +
<<<<<<< HEAD
    theme_bw()
=======
    theme_bw()
  
   ## stabiltiy vs date 
  ggplot(
    brooks_daily_drivers,
    aes(
      date,
      stability_daily
    )
  ) +
    geom_point() +
    geom_smooth(
      method = "gam"
    ) +
    theme_bw()
  
  ## stability vs air temp 
  #Does stability increase as air temperature increases?
  ggplot(
  brooks_daily_drivers,
  aes(
    air_temp_mean_c,
    stability_daily
  )
  ) +
  geom_point() +
  geom_smooth(
    method = "gam"
  ) +
  theme_bw()
  

  # Stability vs mean wind
  
  ggplot(
    brooks_daily_drivers,
    aes(
      wind_speed_mean_ms,
      stability_daily
    )
  ) +
    geom_point() +
    geom_smooth(
      method = "gam"
    ) +
    theme_bw()
  
  
  # Stability vs gust speed
  
  ggplot(
    brooks_daily_drivers,
    aes(
      gust_speed_max_ms,
      stability_daily
    )
  ) +
    geom_point() +
    geom_smooth(
      method = "gam"
    ) +
    theme_bw()
  # stability 
  drivers_cor <- brooks_daily_drivers %>%
    select(
      stability_daily,
      air_temp_mean_c,
      wind_speed_mean_ms,
      gust_speed_max_ms,
      snotel_prcp_mm,
      deep_do_mgl
    )
  
  cor(
    drivers_cor,
    use = "pairwise.complete.obs"
  )
  ## stability and wind mean over time 
  ggplot(
    brooks_daily_drivers,
    aes(date)
  ) +
    geom_line(aes(y = stability_daily)) +
    geom_line(
      aes(y = wind_speed_mean_ms * 10),
      color = "red"
    ) +
    theme_bw()


  
  ####
  # ======================================================================
  # Surface and bottom scaled story plots
  # Nutrients + stability + Chl-a/Secchi + cyano cells + toxins
  # Toxins are matched to nutrient dates within +/- 2 days
  # ======================================================================
  

  
  # ----------------------------------------------------------------------
  # 1. Brooks nutrients: surface and bottom only
  # ----------------------------------------------------------------------
  
  nutrients_brooks <- deq_nutrients_clean_2025 %>%
    mutate(
      date = as.Date(date),
      type = str_to_lower(type)
    ) %>%
    filter(
      lake == "Brooks Lake",
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
  
  # ----------------------------------------------------------------------
  # 2. Daily stability
  # ----------------------------------------------------------------------
  
  stability_daily_plot <- schmidt_daily %>%
    mutate(date = as.Date(date)) %>%
    select(
      date,
      stability_daily,
      strat_state
    )
  
  # ----------------------------------------------------------------------
  # 3. Cyanobacteria cell density
  # ----------------------------------------------------------------------
  
  cyano_density <- phyto_clean %>%
    mutate(
      date = as.Date(date),
      division = str_to_lower(division)
    ) %>%
    filter(
      lake == "Brooks Lake",
      division %in% c("cyanophyta", "cyanobacteria"),
      sample_type != "duplicate"
    ) %>%
    group_by(date) %>%
    summarise(
      cyano_cells = sum(total_cells, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ----------------------------------------------------------------------
  # 4. Toxins: keep surface and depth separate
  # ----------------------------------------------------------------------
  
  tox_brooks_depth <- grab_tox %>%
    mutate(
      date = as.Date(sample_date),
      lake = str_to_lower(lake),
      site_type = str_trim(str_to_lower(site_type))
    ) %>%
    filter(
      lake %in% c("brooks", "brooks lake"),
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
  
  # ----------------------------------------------------------------------
  # 5. Function: nearest toxin date join within +/- 2 days
  # ----------------------------------------------------------------------
  
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
        strat_state,
        total_mc,
        max_total_mc,
        toxin_sample_date = date.y,
        toxin_day_diff = day_diff
      )
  }
  
  # ----------------------------------------------------------------------
  # ----------------------------------------------------------------------
  # 6. Surface-only Chl-a and Secchi
  # ----------------------------------------------------------------------
  
  surface_biology <- nutrients_brooks %>%
    filter(type == "surface") %>%
    select(
      date,
      chla,
      secchi
    )
  
  # ----------------------------------------------------------------------
  # 7. Build story tables without toxins
  # ----------------------------------------------------------------------
  
  story_surface_base <- nutrients_brooks %>%
    filter(type == "surface") %>%
    left_join(
      cyano_density,
      by = "date"
    ) %>%
    left_join(
      stability_daily_plot,
      by = "date"
    )
  
  story_bottom_base <- nutrients_brooks %>%
    filter(type == "bottom") %>%
    select(
      -chla,
      -secchi
    ) %>%
    left_join(
      surface_biology,
      by = "date"
    ) %>%
    left_join(
      cyano_density,
      by = "date"
    ) %>%
    left_join(
      stability_daily_plot,
      by = "date"
    )
  
  # ----------------------------------------------------------------------
  # 8. Join nearest toxin samples within +/- 2 days
  # ----------------------------------------------------------------------
  
  story_surface <- join_nearest_toxin(
    story_surface_base,
    tox_brooks_depth %>% filter(type == "surface"),
    max_days = 2
  )
  
  story_bottom <- join_nearest_toxin(
    story_bottom_base,
    tox_brooks_depth %>% filter(type == "bottom"),
    max_days = 2
  )
  
  # Check matching offsets
  story_surface %>% count(toxin_day_diff)
  story_bottom %>% count(toxin_day_diff)
  # ----------------------------------------------------------------------
  # 8. Function: make scaled seasonal plot
  # ----------------------------------------------------------------------
  
  make_scaled_story_plot <- function(data, title_text) {
    
    plot_long <- data %>%
      select(
        date,
        stability_daily,
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
      
      # Reverse Secchi so higher scaled values mean lower clarity
      group_by(variable) %>%
      mutate(
        value_plot = case_when(
          variable == "secchi" ~ max(value, na.rm = TRUE) - value,
          TRUE ~ value
        ),
        scaled_value = value_plot / max(value_plot, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      mutate(
        variable = recode(
          variable,
          stability_daily = "Schmidt stability",
          ammonia = "Ammonia",
          tn = "Total nitrogen",
          tp = "Total phosphorus",
          chla = "Chl-a",
          secchi = "Lower clarity",
          cyano_cells = "Cyanobacteria cells",
          total_mc = "Total microcystins"
        ),
        line_type = case_when(
          variable == "Schmidt stability" ~ "solid",
          variable %in% c("Ammonia", "Total nitrogen", "Total phosphorus") ~ "dashed",
          variable %in% c(
            "Chl-a",
            "Lower clarity",
            "Cyanobacteria cells",
            "Total microcystins"
          ) ~ "dotted",
          TRUE ~ "solid"
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
  
  p_surface_story <- make_scaled_story_plot(
    story_surface,
    "Brooks Surface seasonal trends"
  )
  
  p_bottom_story <- make_scaled_story_plot(
    story_bottom,
    "Brooks Bottom seasonal trends"
  )
  
  p_surface_story
  p_bottom_story
>>>>>>> fffb24a (:))
