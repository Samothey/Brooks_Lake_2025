# Load required packages -----------------------------------------------

# readxl = reads Excel files
# dplyr = data wrangling tools
# janitor = cleans messy column names
# writexl = writes Excel output files

library(readxl)
library(dplyr)
library(janitor)
library(writexl)

# 01_build_master_toxin_sheet.R
# Purpose: Build one clean master toxin dataset from ALL_TOX + Sample Status sheets


# File paths -----------------------------------------------------------

input_path <- "~/Desktop/Project/Brooks_lake_2025/data_raw/toxins/Pena_2025.xlsx"

output_rds <- "~/Desktop/Project/Brooks_lake_2025/data_clean/toxins/master_tox.rds"
output_csv <- "~/Desktop/Project/Brooks_lake_2025/data_clean/toxins/master_tox.csv"


# Read Excel sheets ----------------------------------------------------

all_tox <- read_excel(input_path, sheet = "ALL_TOX") %>%
  clean_names()

sample_status_raw <- read_excel(input_path, sheet = "Sample Status") %>%
  clean_names()


# Site lookup table ----------------------------------------------------

site_lookup <- tibble::tribble(
  ~site_code, ~lake, ~site_type,
  "BKS_BL_SH_01", "brooks", "shore",
  "BKS_BL_SH_02", "brooks", "shore",
  "BKS_BL_SH_03", "brooks", "shore",
  "BKS_BL_SH_04", "brooks", "shore",
  "BKS_BL_SH_05", "brooks", "shore",
  "BKS_BL_SH_06", "brooks", "shore",
  "BKS_BL_SH_07", "brooks", "shore",
  "BKS_BL_BU_SS", "brooks", "buoy_surface",
  "BKS_BL_BU_DD", "brooks", "buoy_depth",
  "BKS_BL_BON", "tributary", "tributary",
  "BKS_BL_CREEK", "tributary", "tributary",
  "BKS_UB_NCOVE", "upper brooks", "shore",
  "BKS_UB_SH_08", "upper brooks", "shore",
  "BKS_UB_SH_09", "upper brooks", "shore",
  "BKS_UB_BU_SS", "upper brooks", "buoy_surface",
  "BKS_UB_BU_DD", "upper brooks", "buoy_depth",
  "BKS_RN_SH_10", "rainbow", "shore",
  "BKS_RN_BU_SS", "rainbow", "buoy_surface",
  "BKS_RN_BU_DD", "rainbow", "buoy_depth",
  "BKS_LJ_BU_SS", "lower jade", "buoy_surface",
  "BKS_LJ_BU_DD", "lower jade", "buoy_depth",
  "BYS_BR_SH_01", "boysen", "shore",
  "BYS_PC_SH_02", "boysen", "shore",
  "BYS_FR_SH_03", "boysen", "shore",
  "BYS_CR_SH_04", "boysen", "shore",
  "BYS_BY_BU_SS", "boysen", "buoy_surface",
  "BYS_BY_BU_DD", "boysen", "buoy_depth",
  "BKS_RN_INFLOW", "rainbow", "tributary",
  "GRAB_BLANK", "blank", "blank", 
  "Boysen Dup", "boysen", "blank",
  "SPATT_BLANK", "blank", "blank"
)


# Toxin columns --------------------------------------------------------

toxin_cols <- c(
  "rr", "yr", "lr", "la", "dm_lr", "ly",
  "nod", "lf", "wr", "atx", "hatx", "cyl"
)

microcystin_cols <- c(
  "rr", "yr", "lr", "la", "dm_lr", "ly", "lf", "wr"
)


# Build master dataset -------------------------------------------------
master_tox <- sample_status_raw %>%
  
  full_join(
    all_tox,
    by = "sample_number",
    suffix = c("_status", "_tox")
  ) %>%
  
  mutate(
    site_id = coalesce(site_id_tox, site_id_status),
    sample_date = coalesce(sample_date_tox, sample_date_status),
    method = coalesce(method_tox, method_status),
    sample_type = coalesce(sample_type_tox, sample_type_status)
  ) %>%
  
  left_join(
    site_lookup,
    by = c("site_id" = "site_code")
  ) %>%
  
  mutate(
    across(
      all_of(toxin_cols),
      ~ readr::parse_number(as.character(.x))
    ),
    
    result = if_any(
      all_of(toxin_cols),
      ~ !is.na(.x) & .x > 0
    ),
    
    sample_status = case_when(
      if_any(all_of(toxin_cols), ~ !is.na(.x)) & !result ~ "analyzed_no_detection",
      result ~ "detected",
      TRUE ~ "not_analyzed"
    ),
    
    total_mc = rowSums(
      across(all_of(microcystin_cols)),
      na.rm = TRUE
    ),
    
    toxin_richness = rowSums(
      across(all_of(toxin_cols), ~ !is.na(.x) & .x > 0),
      na.rm = TRUE
    )
  ) %>%
  
  select(
    sample_number,
    site_id,
    lake,
    site_type,
    sample_date,
    sample_type,
    method,
    result,
    sample_status,
    status,
    total_mc,
    toxin_richness,
    all_of(toxin_cols),
    notes
  )

# Save outputs ---------------------------------------------------------

saveRDS(master_tox, output_rds)

write.csv(
  master_tox,
  output_csv,
  row.names = FALSE
)


# View dataset ---------------------------------------------------------

View(master_tox)
