# - brooks lake buoy data #

library(dplyr)
library(stringr)
library(lubridate)
library(rLakeAnalyzer)

brooks_clean <- read.csv("~/Desktop/Project/Brooks_lake_2025/data_clean/Brooks_WQ_2025_cleaned.csv")

brooks <- brooks_clean %>%
  filter(str_detect(DATETIME_MST, ":")) %>%   # <- filter MST column
  mutate(
    dt_mst = ymd_hms(DATETIME_MST, tz = "Etc/GMT+7")
  )

start_mst <- ymd_hms("2025-06-12 02:05:00", tz = "Etc/GMT+7")

brooks <- brooks %>%
  filter(dt_mst >= start_mst)

min(brooks$dt_mst)

temp_cols <- c("WTemp1m","WTemp4m","WTemp7m","WTemp9m","WTemp10m","WTemp13m","WTemp15m")

brooks_daily <- brooks %>%
  mutate(date = as.Date(dt_mst)) %>%
  group_by(date) %>%
  summarize(
    across(all_of(temp_cols), ~ mean(.x, na.rm = TRUE)),
    DO_15m = mean(odomgL_15m, na.rm = TRUE),
    .groups = "drop"
  )



depths_m <- c(1, 4, 7, 9, 10, 13, 15)

wtr_mat <- as.matrix(brooks_daily %>% select(all_of(temp_cols)))

thermo <- apply(wtr_mat, 1, function(wtr) {
  if (sum(!is.na(wtr)) < 4) return(NA_real_)
  thermo.depth(wtr = wtr, depths = depths_m)
})

max_dTdz <- apply(wtr_mat, 1, function(wtr) {
  if (sum(!is.na(wtr)) < 4) return(NA_real_)
  grads <- diff(wtr) / diff(depths_m)
  max(abs(grads), na.rm = TRUE)
})

brooks_state <- brooks_daily %>%
  mutate(
    thermo_depth_m = thermo,
    deltaT_C = WTemp1m - WTemp15m,
    max_dTdz_C_per_m = max_dTdz,
    strat_state = case_when(
      !is.na(thermo_depth_m) & deltaT_C >= 1 ~ "Stratified",
      TRUE ~ "Mixed/Weak"
    )
  )

summary(brooks_state$thermo_depth_m)
table(brooks_state$strat_state)

library(ggplot2)

ggplot(brooks_state, aes(date, thermo_depth_m)) +
  geom_line() +
  labs(title = "Brooks Lake Thermocline Depth", y = "Depth (m)", x = NULL) +
  theme_minimal()