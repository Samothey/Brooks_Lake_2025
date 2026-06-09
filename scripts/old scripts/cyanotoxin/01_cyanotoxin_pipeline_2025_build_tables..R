
# 2025 Cyanotoxin Data ----------------------------------------------------
# Author: Samantha Peña
#
# Goals -------------------------------------------------------------------
#   - Read toxin results (WholeWater + SPATT) and sample tracking (Sample Status)
#   - Standardize fields (dates, casing, site metadata)
#   - Build:
#       1) master_sample_tracking (status_final + "missing/analyzed" logic)
#       2) toxin_results_wide_all (grab + spatt combined)
#       3) toxin_results_long_regular (regular samples only, pivoted long)


# ---- Packages ----
library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)

# ---- File paths  ----
main_toxin <- "C:/Users/spena5/Desktop/Project/Brooks_lake_2025/CODE/cyanotoxin/Cyanotoxins_2025/Pena_2025.xlsx"

out_dir   <- "C:/Users/spena5/Desktop/Project/Brooks_lake_2025/CODE/cyanotoxin/Cyanotoxins_2025"
out_master_rds <- file.path(out_dir, "master_sample_tracking.rds")
out_master_csv <- file.path(out_dir, "master_sample_tracking.csv")
out_wide_rds   <- file.path(out_dir, "toxin_results_wide_all.rds")
out_wide_csv   <- file.path(out_dir, "toxin_results_wide_all.csv")
out_long_rds   <- file.path(out_dir, "toxin_results_long_regular.rds")
out_long_csv   <- file.path(out_dir, "toxin_results_long_regular.csv")
out_grab_totals_rds   <- file.path(out_dir, "grab_class_totals.rds")
out_grab_totals_csv   <- file.path(out_dir, "grab_class_totals.csv")
out_spatt_presence_rds <- file.path(out_dir, "spatt_presence.rds")
out_spatt_presence_csv <- file.path(out_dir, "spatt_presence.csv")


# ---- Rules/constants ----
exclude_site_patterns <- c("^BKS_UB_NCOVE")  # add more as needed

site_type_levels <- c("shore", "buoy")
depth_levels     <- c("surface", "bottom")

mc_threshold <- 0.8  # EPA/ Wyoming  threshold for recreational advsiory ug/L


# filter extras
  filter_extras <- function(df) {
  df %>%
    filter(!if_any(site_id, ~ str_detect(.x, exclude_site_patterns)))
}

# derive lake/site_type/depth from site_id
add_site_fields <- function(df) {
  df %>%
    mutate(
      lake = case_when(
        str_detect(site_id, "^BKS_BL_") ~ "brooks",
        str_detect(site_id, "^BKS_UB_") ~ "upper brooks",
        str_detect(site_id, "^BKS_RN_") ~ "rainbow",
        str_detect(site_id, "^BKS_LJ_") ~ "lower jade",
        str_detect(site_id, "^BYS_")    ~ "boysen",
        str_detect(site_id, "^GRAB_")   ~ "blank",
        str_detect(site_id, "^SPATT_")  ~ "blank",
        TRUE ~ NA_character_
      ),
      site_type = case_when(
        str_detect(site_id, "_SH_") ~ "shore",
        str_detect(site_id, "_BU_") ~ "buoy",
        TRUE ~ NA_character_
      ),
      depth_category = case_when(
        str_detect(site_id, "_SS$") ~ "surface",
        str_detect(site_id, "_DD$") ~ "bottom",
        TRUE ~ NA_character_
      )
    )
}

# ---- standardize excel sheets  ----
# Function: read_toxin_sheet()
# Purpose:Reads a single sheet from the main cyanotoxin Excel workbook and applies
#   consistent cleaning and standardization steps across all datasets, (WholeWater, SPATT, Sample Status).
# Why create this function?
#   - Avoids duplicating data cleaning code for each sheet
#   - Ensures consistent formatting for joins, summaries, and plots
#   - Makes the analysis pipeline more reproducible and easier to maintain
#
# Arguments used:
#   main_toxin      = file path to the Excel workbook
#   sheet_name      = name of the sheet to read
#   method_override = "grab" or "spatt" to force method assignment,
#                     or NULL to read method from the sheet
#   status_sheet    = TRUE if the sheet is the Sample Status sheet
#
# Returns: A cleaned tibble with standardized fields and derived site metadata
#
# NOTE: This function structure was developed with assistance from ChatGPT to improve code modularity and documentation.


read_toxin_sheet <- function(main_toxin, sheet_name,
                             method_override = NULL,
                             status_sheet = FALSE) {
  
  df <- read_excel(main_toxin, sheet = sheet_name) %>%
    clean_names()
  
  # Standardize sample_number 
  if ("sample_number" %in% names(df)) {
    df <- df %>%
      mutate(sample_number = str_pad(as.character(sample_number), width = 3, pad = "0"))
  }
  
  # Standardize date column name
  if ("sample_date" %in% names(df)) {
    df <- df %>% rename(date = sample_date)
  }
  
  # Standardize core fields that are common across sheets
  df <- df %>%
    mutate(
      date = as.Date(date),
      sample_type = if ("sample_type" %in% names(df)) str_squish(str_to_lower(sample_type)) else NA_character_
    )
  
  # Apply method logic
  if (!is.null(method_override)) {
    df <- df %>% mutate(method = method_override)
  } else if ("method" %in% names(df)) {
    df <- df %>% mutate(method = str_squish(str_to_lower(method)))
  }
  
  # Status sheet extra standardization
  if (status_sheet && "status" %in% names(df)) {
    df <- df %>% mutate(status = str_squish(str_to_lower(status)))
  }
  
  # Apply shared cleaning/derivations
  df %>%
    filter_extras() %>%
    add_site_fields()
}

#read all sheets
grab_result <- read_toxin_sheet(main_toxin, "WholeWater", method_override = "grab")
spatt_result <- read_toxin_sheet(main_toxin, "SPATT", method_override = "spatt")
sample_status <- read_toxin_sheet(main_toxin, "Sample Status", status_sheet = TRUE)


# master sample tracking  -------------------------------------------------

# Keys that indicates whether a result exists for that grab or spatt sample
grab_keys <- grab_result %>%
  distinct(sample_number, site_id, date, sample_type, method) %>%
  mutate(has_result = TRUE)

spatt_keys <- spatt_result %>%
  distinct(sample_number, site_id, date, sample_type, method) %>%
  mutate(has_result = TRUE)

# Join keys onto status and define status_final
master_sample_tracking <- sample_status %>%
  left_join(grab_keys,  by = c("sample_number", "site_id", "date", "sample_type", "method")) %>%
  mutate(has_grab_result = coalesce(has_result, FALSE)) %>%
  select(-has_result) %>%
  left_join(spatt_keys, by = c("sample_number", "site_id", "date", "sample_type", "method")) %>%
  mutate(has_spatt_result = coalesce(has_result, FALSE)) %>%
  select(-has_result) %>%
  mutate(
    status_final = case_when(
      method == "grab"  & has_grab_result  ~ "analyzed",
      method == "spatt" & has_spatt_result ~ "analyzed",
      status == "shipped" & method == "grab"  & !has_grab_result  ~ "missing",
      status == "shipped" & method == "spatt" & !has_spatt_result ~ "missing",
      TRUE ~ status
    )
  )

# Quick summary (what do we have )
status_summary <- master_sample_tracking %>%
  filter(status_final %in% c("missing", "lost", "analyzed")) %>%
  count(method, status_final) %>%
  pivot_wider(names_from = status_final, values_from = n, values_fill = 0)

print(status_summary)

# Save outputs
saveRDS(master_sample_tracking, out_master_rds)
write_csv(master_sample_tracking, out_master_csv)


# Combine toxin results ------------

## STEP 1: Combine GRAB and SPATT results into one table 

toxin_results_wide_all <- bind_rows(grab_result, spatt_result) %>%
  # Ensure text fields are consistently formatted
  mutate(
    method = str_to_lower(method),
    sample_type = str_squish(str_to_lower(sample_type)) # trims spaces + lowercases
  )
# save
write_csv(toxin_results_wide_all, out_wide_csv)
saveRDS(toxin_results_wide_all, out_wide_rds)

# STEP 2: Define which columns describe the sample (metadata)

meta_cols <- c("sample_number", # unique identifier
"site_id",  #site code
"date",        #collection date
"sample_type",  # regular, dup, blank
"method", # spatt vs grab
"lake",     # lake name 
"site_type",     # shore vs buoy
"depth_category") # surafce vs bottom 

# # STEP 3: Automatically identify toxin concentration columns 
#   - Toxin measurements should be numeric
#   - Metadata columns may also be numeric, so remove them explicitly
# Result: toxin_cols` is a character vector of column names corresponding to toxins

toxin_cols <- toxin_results_wide_all %>%
  select(where(is.numeric)) %>%  # select all numeric columns
  names() %>%                  # extract their names
  setdiff(meta_cols)          # remove metadata fields if any are numeric


# STEP 4: Convert to long format (regular samples only)
# Long format structure:
#   - one row = one sample × one toxin congener
#   - required for tidy summaries, ggplot, and toxin class aggregation

toxin_results_long_regular <- toxin_results_wide_all %>%
  # Keep only "real" environmental samples
  filter(
    sample_type == "regular", # exclude blanks, duplicates, QA/QC samples
    !is.na(lake),             # drop samples without a valid lake assignment
    !lake %in% c("blank", "duplicate")
  ) %>%
  
  # Reshapes from wide -> long:
  #   LR, RR, YR, ATX, ...  -->  congener + value columns
  pivot_longer(
    cols = all_of(toxin_cols), # only pivot toxin columns
    names_to = "congener", # toxin name (e.g., "lr", "rr")
    values_to = "value" # measured concentration
  ) %>%
  # Standardize toxin names for joins and plotting
  mutate(congener = str_squish(str_to_lower(congener)))

write_csv(toxin_results_long_regular, out_long_csv)
saveRDS(toxin_results_long_regular, out_long_rds)

# Toxin key 
toxin_key <- tribble(
  ~congener, ~toxin_class,
  "dm_lr", "microcystin",
  "la",    "microcystin",
  "lf",    "microcystin",
  "lr",    "microcystin",
  "ly",    "microcystin",
  "rr",    "microcystin",
  "wr",    "microcystin",
  "yr",    "microcystin",
  "atx",   "anatoxin",
  "hatx",  "anatoxin",
  "cyl",   "cylindrospermopsin",
  "nod",   "nodularin"
)

# Save long-format toxin table
toxin_plot <- toxin_results_long_regular %>%
  left_join(toxin_key, by = "congener")

toxin_plot <- toxin_plot %>%
  mutate(
    lake = str_squish(str_to_lower(lake)),
    site_type = str_squish(str_to_lower(site_type)),
    depth_category = str_squish(str_to_lower(depth_category)),
    method = str_squish(str_to_lower(method))
  )



# STEP 5: Build toxin-class summaries used for plots
#   - GRAB: total concentration by toxin class (sum congeners)
#   - SPATT: detection by toxin class (any congener detected)

# ---- GRAB totals by toxin class -----------------------------------------
# One row per: sample_number × site_id × date × depth_category × toxin_class
grab_class_totals <- toxin_plot %>%
  filter(method == "grab") %>%                 # only GRAB samples
  filter(!is.na(toxin_class)) %>%              # drop unmapped congeners
  group_by(
    sample_number, site_id, date,              # sample identifiers
    lake, site_type, depth_category,           # location + depth context
    method, toxin_class                        # summarize within toxin class
  ) %>%
  summarise(
    total = sum(value, na.rm = TRUE),          # sum congeners in that class
    .groups = "drop"
  )

# ---- SPATT presence/absence by toxin class ------------------------------
# One row per: sample_number × site_id × date × depth_category × toxin_class
spatt_presence <- toxin_plot %>%
  filter(method == "spatt") %>%                # only SPATT samples
  filter(!is.na(toxin_class)) %>%              # drop unmapped congeners
  group_by(
    sample_number, site_id, date,
    lake, site_type, depth_category,
    method, toxin_class
  ) %>%
  summarise(
    detected = any(value > 0, na.rm = TRUE),   # TRUE if any congener detected
    .groups = "drop"
  )

# Save for plotting script
saveRDS(grab_class_totals, out_grab_totals_rds)
write_csv(grab_class_totals, out_grab_totals_csv)

saveRDS(spatt_presence, out_spatt_presence_rds)
write_csv(spatt_presence, out_spatt_presence_csv)

