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
#   of continuous daily buoy stability.
# ======================================================================

library(tidyverse)
library(lubridate)
library(rLakeAnalyzer)

# ======================================================================
# 1. File paths
# ======================================================================

data_out_dir <- "data_clean/analysis"
fig_dir <- "figures/integrated_analysis/rainbow"

dir.create(data_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 2. Load data
# ======================================================================

brooks_weather_noaa_daily_2025 <- readRDS(
  "data_clean/weather/brooks_weather_noaa_daily_2025.rds"
)

lake_profile_compiled <- read_excel("data_clean/deq/lake_profile_compiled.xlsx")

# ======================================================================
# 3. Clean Rainbow profile data
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
# 4. Rainbow bathymetry and target depths
# ======================================================================

acre_to_m2 <- 4046.856

target_depths_rn <- c(
  0.75, 1, 2, 3, 4, 5.5, 7, 8
)

bthD_rn <- c(
  0.75, 1, 2, 3, 4, 5.5, 7, 8
)

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
# 5. Profile summary: temperature and DO gradients
# ======================================================================

rainbow_profile_summary <- rainbow_profiles %>%
  group_by(date) %>%
  summarise(
    n_depths = n(),
    
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
# 6. Profile-based Schmidt stability
# ======================================================================

rainbow_profile_stability <- rainbow_profiles %>%
  group_by(date) %>%
  group_modify(~ {
    
    profile <- .x %>%
      arrange(depth) %>%
      filter(
        !is.na(depth),
        !is.na(temp)
      ) %>%
      distinct(depth, .keep_all = TRUE)
    
    min_profile_depth <- min(profile$depth, na.rm = TRUE)
    max_profile_depth <- max(profile$depth, na.rm = TRUE)
    
    usable_depths <- target_depths_rn[
      target_depths_rn >= min_profile_depth &
        target_depths_rn <= max_profile_depth
    ]
    
    if (nrow(profile) < 3 || length(usable_depths) < 3) {
      return(
        tibble(
          n_profile_depths = nrow(profile),
          min_profile_depth = min_profile_depth,
          max_profile_depth = max_profile_depth,
          stability_profile = NA_real_
        )
      )
    }
    
    temp_interp <- approx(
      x = profile$depth,
      y = profile$temp,
      xout = usable_depths,
      rule = 2
    )$y
    
    stability_value <- tryCatch(
      {
        schmidt.stability(
          wtr = temp_interp,
          depths = usable_depths,
          bthD = bthD_rn,
          bthA = bthA_rn
        )
      },
      error = function(e) NA_real_
    )
    
    tibble(
      n_profile_depths = nrow(profile),
      min_profile_depth = min_profile_depth,
      max_profile_depth = max_profile_depth,
      stability_profile = as.numeric(stability_value)[1]
    )
  }) %>%
  ungroup()
# ======================================================================
# 7. Join profile drivers with weather
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

# ======================================================================
# 8. Define relative profile stability states
# ======================================================================

stability_q_rn <- quantile(
  rainbow_profile_drivers$stability_profile,
  probs = c(0.33, 0.66),
  na.rm = TRUE
)

rainbow_profile_drivers <- rainbow_profile_drivers %>%
  mutate(
    strat_state = case_when(
      is.na(stability_profile) ~ NA_character_,
      stability_profile <= stability_q_rn[1] ~ "Low stability",
      stability_profile <= stability_q_rn[2] ~ "Moderate stability",
      stability_profile > stability_q_rn[2] ~ "High stability",
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
# 9. Save outputs
# ======================================================================

saveRDS(
  rainbow_profiles,
  file.path(data_out_dir, "rainbow_profiles_2025.rds")
)

saveRDS(
  rainbow_profile_summary,
  file.path(data_out_dir, "rainbow_profile_summary_2025.rds")
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
# 10. Quick checks
# ======================================================================
