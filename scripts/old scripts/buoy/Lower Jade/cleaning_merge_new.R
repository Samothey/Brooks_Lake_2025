library(dplyr)
library(lubridate)
library(readr)
library(stringr)

clean_hobo_temp <- function(file, depth_m) {
  raw <- read.csv(file, skip = 1, header = TRUE)
  
  datetime_col <- names(raw)[str_detect(names(raw), "Date.Time")]
  temp_col     <- names(raw)[str_detect(names(raw), "Temp")]
  
  cleaned <- raw %>%
    transmute(
      datetime_local = mdy_hms(.data[[datetime_col]], tz = "America/Denver"),
      datetime_utc   = with_tz(datetime_local, "UTC"),
      temp_c         = (.data[[temp_col]] - 32) * 5/9
    ) %>%
    select(datetime_utc, temp_c)
  
  names(cleaned)[2] <- paste0("temp_", depth_m, "m")
  
  cleaned
}

temp_2m_clean  <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_2m.csv",  2)
temp_4m_clean  <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_4m.csv",  4)
temp_5m_clean  <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_5m.csv",  5)
temp_6m_clean  <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_6m.csv",  6)
temp_7m_clean  <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_7m.csv",  7)
temp_8m_clean  <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_8m.csv",  8)
temp_10m_clean <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_10m.csv", 10)
temp_11m_clean <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_11m.csv", 11)
temp_12m_clean <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_12m.csv", 12)
temp_13m_clean <- clean_hobo_temp("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_13m.csv", 13)

lower_jade_hobo <- temp_2m_clean %>%
  full_join(temp_4m_clean,  by = "datetime_utc") %>%
  full_join(temp_5m_clean,  by = "datetime_utc") %>%
  full_join(temp_6m_clean,  by = "datetime_utc") %>%
  full_join(temp_7m_clean,  by = "datetime_utc") %>%
  full_join(temp_8m_clean,  by = "datetime_utc") %>%
  full_join(temp_10m_clean, by = "datetime_utc") %>%
  full_join(temp_11m_clean, by = "datetime_utc") %>%
  full_join(temp_12m_clean, by = "datetime_utc") %>%
  full_join(temp_13m_clean, by = "datetime_utc") %>%
  arrange(datetime_utc)


do_surface <- read.table("~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_Surface.TXT", sep=",", header=TRUE, skip = 7) 

do_surface_clean <- do_surface %>%
  slice(-1) %>%
  transmute(
    datetime_utc = ymd_hms(trimws(UTC_Date_._Time), tz = "UTC"),
    temp_1m      = as.numeric(trimws(Temperature)),
    do_mgl_1m    = as.numeric(trimws(Dissolved.Oxygen)),
    do_sat_1m    = as.numeric(trimws(Dissolved.Oxygen.Saturation))
  ) %>%
  arrange(datetime_utc)


head(do_surface_clean)
str(do_surface_clean)
range(do_surface_clean$datetime_utc, na.rm = TRUE)
summary(do_surface_clean$temp_1m)
summary(do_surface_clean$do_mgl_1m)
summary(do_surface_clean$do_sat_1m)


do_depth <- read.table(
  "~/Desktop/Project/Brooks_lake_2025/data_raw/buoy/lowerjade/LJ_Depth.TXT",
  sep = ",",
  header = TRUE,
  skip = 7
)

do_depth_clean <- do_depth %>%
  slice(-1) %>%
  transmute(
    datetime_utc = ymd_hms(trimws(UTC_Date_._Time), tz = "UTC"),
    temp_14m     = as.numeric(trimws(Temperature)),
    do_mgl_14m   = as.numeric(trimws(Dissolved.Oxygen)),
    do_sat_14m   = as.numeric(trimws(Dissolved.Oxygen.Saturation))
  ) %>%
  arrange(datetime_utc)


head(do_depth_clean)
str(do_depth_clean)
range(do_depth_clean$datetime_utc, na.rm = TRUE)
summary(do_depth_clean$temp_14m)
summary(do_depth_clean$do_mgl_14m)
summary(do_depth_clean$do_sat_14m)


library(dplyr)
library(purrr)

hobo_list <- list(
  temp_2m_clean,
  temp_4m_clean,
  temp_5m_clean,
  temp_6m_clean,
  temp_7m_clean,
  temp_8m_clean,
  temp_10m_clean,
  temp_11m_clean,
  temp_12m_clean,
  temp_13m_clean
)

lower_jade_hobo <- reduce(hobo_list, full_join, by = "datetime_utc") %>%
  arrange(datetime_utc)

names(lower_jade_hobo)
head(lower_jade_hobo)
range(lower_jade_hobo$datetime_utc, na.rm = TRUE)

lower_jade_raw_merged <- lower_jade_hobo %>%
  full_join(do_surface_clean, by = "datetime_utc") %>%
  full_join(do_depth_clean,   by = "datetime_utc") %>%
  arrange(datetime_utc)

lower_jade_raw_merged <- lower_jade_raw_merged %>%
  select(
    datetime_utc,
    temp_1m, do_mgl_1m, do_sat_1m,
    temp_2m,
    temp_4m,
    temp_5m,
    temp_6m,
    temp_7m,
    temp_8m,
    temp_10m,
    temp_11m,
    temp_12m,
    temp_13m,
    temp_14m, do_mgl_14m, do_sat_14m
  )



lower_jade_raw_merged <- lower_jade_raw_merged %>%
  mutate(datetime_mst = with_tz(datetime_utc, "America/Denver")) %>%
  relocate(datetime_mst, .after = datetime_utc)


lower_jade_raw_merged %>%
  count(datetime_utc) %>%
  filter(n > 1)

time_gaps <- lower_jade_raw_merged %>%
  arrange(datetime_utc) %>%
  mutate(dt_mins = as.numeric(difftime(datetime_utc, lag(datetime_utc), units = "mins")))

table(time_gaps$dt_mins, useNA = "ifany")

colSums(is.na(lower_jade_raw_merged))



library(dplyr)
library(lubridate)

lower_jade_hourly <- lower_jade_raw_merged %>%
  mutate(datetime_hour_utc = floor_date(datetime_utc, unit = "hour")) %>%
  group_by(datetime_hour_utc) %>%
  summarise(
    across(-c(datetime_utc, datetime_mst), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    across(-datetime_hour_utc, ~ ifelse(is.nan(.x), NA, .x)),
    datetime_hour_mst = with_tz(datetime_hour_utc, "America/Denver")
  ) %>%
  relocate(datetime_hour_mst, .after = datetime_hour_utc)



head(lower_jade_hourly)
colSums(is.na(lower_jade_hourly))

ggplot(lower_jade_hourly, aes(x = datetime_hour_mst)) +
  geom_line(aes(y = temp_1m), linewidth = 0.4) +
  geom_line(aes(y = temp_2m), linewidth = 0.4) +
  geom_line(aes(y = temp_14m), linewidth = 0.4) +
  labs(
    title = "Lower Jade Hourly Temperatures",
    x = "Date",
    y = "Temperature (°C)"
  ) +
  theme_minimal()