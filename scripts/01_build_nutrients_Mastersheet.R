```r
library(tidyverse)
library(lubridate)
library(janitor)

# -------------------------------------------------------------------
# Load DEQ dataset
# -------------------------------------------------------------------

deq_raw <- read_csv("data_raw/deq/LakeData_Brooks2009-2024_Rstudio.csv")

# -------------------------------------------------------------------
# Clean + standardize columns
# -------------------------------------------------------------------

deq_nutrients_clean_2025 <- deq_raw |>
  clean_names() |>
  mutate(
    date = mdy(date)
  ) |>
  filter(
    date >= as.Date("2025-06-24")
  ) |>
  select(
    lake = waterbody,
    depth = sample_depth,
    date,
    ammonia,
    ammonia_cen,
    no2no3,
    no2no3_cen,
    tn,
    tn_cen,
    tp,
    tp_cen,
    tn_tp,
    chla,
    chla_cen,
    type,
    secchi
  ) |>
  distinct()
 

# -------------------------------------------------------------------
# Quick checks
# -------------------------------------------------------------------

glimpse(deq_nutrients_clean_2025)

summary(deq_nutrients_clean_2025)

colnames(deq_nutrients_clean_2025)

# -------------------------------------------------------------------
# Save cleaned dataset
# -------------------------------------------------------------------

saveRDS(
  deq_nutrients_clean_2025,
  "~/Desktop/Project/Brooks_lake_2025/data_clean/deq/deq_nutrients_clean_2025.rds"
)

readRDS(
  "~/Desktop/Project/Brooks_lake_2025/data_clean/deq/deq_nutrients_clean_2025.rds"
)