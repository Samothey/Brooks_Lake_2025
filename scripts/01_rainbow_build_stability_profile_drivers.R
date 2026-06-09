
# ======================================================================
# 01_rainbow_build_profile_drivers.R
#
# Purpose:
#   Build Rainbow Lake profile-based physical driver dataset:
#   - Schmidt stability from DEQ vertical profiles
#   - surface-bottom temperature difference
#   - surface-bottom DO difference
#   - shared weather drivers
#
# Note:
#   Rainbow buoy died early, so profile-based stability is used instead
#   of daily high-frequency buoy stability.
# ======================================================================

library(tidyverse)
library(lubridate)
library(rLakeAnalyzer)

# ======================================================================
# 1. File paths
# ======================================================================

data_out_dir <- "data_clean/analysis/rainbow"
fig_dir <- "figures/integrated_analysis/rainbow"

dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load data
# ======================================================================

brooks_weather_noaa_daily_2025 <- readRDS(
  "data_clean/weather/brooks_weather_noaa_daily_2025.rds"
)

# assumes lake_profile_compiled is already saved as RDS
  lake_profile_compiled <- read_excel("data_clean/deq/lake_profile_compiled.xlsx")
  

# ======================================================================
# 3. Rainbow profile data
# ======================================================================

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
    !is.na(depth),
    !is.na(temp)
  ) %>%
  arrange(date, depth)

# ======================================================================
# 4. Rainbow bathymetry
# ======================================================================

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

# ======================================================================
# 5. Calculate profile-based Schmidt stability
# ======================================================================

rainbow_profile_stability <- rainbow_profiles %>%
  group_by(date) %>%
  group_modify(~ {
    
    profile <- .x %>%
      arrange(depth) %>%
      filter(
        !is.na(depth),
        !is.na(temp)
      )
    
    tibble(
      n_depths = nrow(profile),
      min_depth = min(profile$depth, na.rm = TRUE),
      max_depth = max(profile$depth, na.rm = TRUE),
      stability_profile = if_else(
        nrow(profile) >= 3,
        schmidt.stability(
          wtr = profile$temp,
          depths = profile$depth,
          bthD = bthD_rn,
          bthA = bthA_rn
        ),
        NA_real_
      )
    )
  }) %>%
  ungroup()

# ======================================================================
# 6. Profile summary: temp and DO gradients
# ======================================================================

rainbow_profile_summary <- rainbow_profiles %>%
  group_by(date) %>%
  summarise(
    surface_depth = min(depth, na.rm = TRUE),
    bottom_depth = max(depth, na.rm = TRUE),
    
    surface_temp = temp[which.min(depth)],
    bottom_temp = temp[which.max(depth)],
    temp_diff = surface_temp - bottom_temp,
    
    surface_do_mgl = do_mg_l[which.min(depth)],
    bottom_do_mgl = do_mg_l[which.max(depth)],
    do_diff = surface_do_mgl - bottom_do_mgl,
    
    .groups = "drop"
  )

# ======================================================================
# 7. Join profile stability + profile gradients + weather
# ======================================================================

rainbow_profile_drivers <- rainbow_profile_summary %>%
  left_join(
    rainbow_profile_stability,
    by = "date"
  ) %>%
  left_join(
    brooks_weather_noaa_daily_2025 %>%
      mutate(date = as.Date(date)),
    by = "date"
  ) %>%
  mutate(
    stability_profile = as.numeric(stability_profile)
  )

stability_q_rn <- quantile(
  rainbow_profile_drivers$stability_profile,
  probs = c(0.33, 0.66),
  na.rm = TRUE
)

rainbow_profile_drivers <- rainbow_profile_drivers %>%
  mutate(
    strat_state = case_when(
      stability_profile <= stability_q_rn[[1]] ~ "Low stability",
      stability_profile <= stability_q_rn[[2]] ~ "Moderate stability",
      stability_profile > stability_q_rn[[2]] ~ "High stability",
      TRUE ~ NA_character_
    ),
    strat_state = factor(
      strat_state,
      levels = c(
        "Low stability",
        "Moderate stability",
        "High stability"
      )
    )
  )

# ======================================================================
# 8. Save outputs
# ======================================================================

saveRDS(
  rainbow_profiles,
  file.path(data_out_dir, "rainbow_profiles_2025.rds")
)

saveRDS(
  rainbow_profile_stability,
  file.path(data_out_dir, "rainbow_profile_stability_2025.rds")
)

saveRDS(
  rainbow_profile_drivers,
  file.path(data_out_dir, "rainbow_profile_drivers_2025.rds")
)

# ======================================================================
# 9. Quick checks
# ======================================================================

glimpse(rainbow_profile_drivers)

rainbow_profile_drivers %>%
  select(
    date,
    n_depths,
    min_depth,
    max_depth,
    stability_profile,
    temp_diff,
    do_diff,
    strat_state
  )
usethis:use_github()
