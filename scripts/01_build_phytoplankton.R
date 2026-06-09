library(tidyverse)
library(lubridate)
library(janitor)
library(scales)

# ------------------------------------------------------------
# Clean phytoplankton data (removing the duplicates)
# ------------------------------------------------------------

# Load raw phytoplankton data ------------------------------------------

phyto_raw <- read_excel(
  "data_raw/phytoplankton/WYDEQ 2025 Lander Phytoplankton Heterocyst Flat Data.xlsx"
)

# Clean phytoplankton dataset ------------------------------------------

phyto_clean <- phyto_raw |>
  clean_names() |>
  select(
    sample_station_name,
    sample_date_collected,
    division,
    taxon,
    per_taxa_id,
    heterocyst_count,
    ncu_counted,
    cells_counted,
    live_cells,
    dead_cells,
    ncu_counted_l,
    cells_counted_l,
    live_cells_l,
    dead_cells_l,
    total_cells_l,
    fields_of_view_analyzed,
    volume_analyzed_l,
    analyst,
    heterocyst_density_l
  ) |>
  rename(
    lake = sample_station_name,
    date = sample_date_collected,
    taxa_id = per_taxa_id,
    heterocysts = heterocyst_count,
    total_cells = total_cells_l,
    live_cells_density = live_cells_l,
    dead_cells_density = dead_cells_l,
    volume_analyzed = volume_analyzed_l,
    heterocyst_density = heterocyst_density_l
  ) |>
  mutate(
    date = mdy(date),
    
    lake = str_to_lower(lake),
    lake = str_replace_all(lake, "\\s*-\\s*", " - "),
    
    # Identify duplicates BEFORE removing "deepest"
    sample_type = case_when(
      str_detect(lake, "dup") ~ "duplicate",
      TRUE ~ "regular"
    ),
    
    # Remove " - deepest" and " - deepest dup" from lake_site
    lake = str_remove(lake, "\\s*-\\s*deepest.*$"),
    
    division = str_to_lower(division),
    taxon = str_squish(taxon),
    analyst = str_squish(analyst),
    
    total_cells = as.numeric(total_cells),
    live_cells_density = as.numeric(live_cells_density),
    dead_cells_density = as.numeric(dead_cells_density),
    heterocyst_density = as.numeric(heterocyst_density),
    
    lake = case_when(
      str_detect(lake, "^upper brooks lake") ~ "Upper Brooks Lake",
      str_detect(lake, "^brooks lake") ~ "Brooks Lake",
      str_detect(lake, "^rainbow lake") ~ "Rainbow Lake",
      str_detect(lake, "^lower jade lake") ~ "Lower Jade Lake",
      str_detect(lake, "^upper jade lake") ~ "Upper Jade Lake",
      TRUE ~ "Other"
    )
  ) |>
  filter(sample_type == "regular")

# QA checks -------------------------------------------------------------

table(phyto_clean$sample_type)
unique(phyto_clean$lake_site)
unique(phyto_clean$lake)
nrow(phyto_clean)


# Save cleaned dataset -------------------------------------------------

dir.create(
  "~/Desktop/Project/Brooks_lake_2025/data_clean/phytoplankton",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  phyto_clean,
  "~/Desktop/Project/Brooks_lake_2025/data_clean/phytoplankton/phyto_clean.rds"
)


# ------------------------------------------------------------
# Clean phytoplankton data ( this keeps the duplicates though)
# ------------------------------------------------------------
WYDEQ_2025_Lander_Phytoplankton_Heterocyst_Flat_Data <- read_excel("data_raw/phytoplankton/WYDEQ 2025 Lander Phytoplankton Heterocyst Flat Data.xlsx")
 

phyto_cleanwdups <- WYDEQ_2025_Lander_Phytoplankton_Heterocyst_Flat_Data |>
  clean_names() |>
  select(
    sample_station_name,
    sample_date_collected,
    division,
    taxon,
    per_taxa_id,
    heterocyst_count,
    ncu_counted,
    cells_counted,
    live_cells,
    dead_cells,
    ncu_counted_l,
    cells_counted_l,
    live_cells_l,
    dead_cells_l,
    total_cells_l,
    fields_of_view_analyzed,
    volume_analyzed_l,
    analyst,
    heterocyst_density_l
  ) |>
  rename(
    lake = sample_station_name,
    date = sample_date_collected,
    taxa_id = per_taxa_id,
    heterocysts = heterocyst_count,
    total_cells = total_cells_l,
    live_cells_density = live_cells_l,
    dead_cells_density = dead_cells_l,
    volume_analyzed = volume_analyzed_l,
    heterocyst_density = heterocyst_density_l
  ) |>
  mutate(
    date = mdy(date),
    
    lake = str_to_lower(lake),
    lake = str_replace_all(lake, "\\s*-\\s*", " - "),
    lake = str_remove(lake, "\\s*-\\s*deepest$"),
    
    division = str_to_lower(division),
    taxon = str_squish(taxon),
    analyst = str_squish(analyst),
    
    total_cells = as.numeric(total_cells),
    live_cells_density = as.numeric(live_cells_density),
    dead_cells_density = as.numeric(dead_cells_density),
    heterocyst_density = as.numeric(heterocyst_density),
    
    lake = case_when(
      str_detect(lake, "^upper brooks lake") ~ "Upper Brooks Lake",
      str_detect(lake, "^brooks lake") ~ "Brooks Lake",
      str_detect(lake, "^rainbow lake") ~ "Rainbow Lake",
      str_detect(lake, "^lower jade lake") ~ "Lower Jade Lake",
      str_detect(lake, "^upper jade lake") ~ "Upper Jade Lake",
      TRUE ~ "Other"
    ),
    
    sample_type = case_when(
      str_detect(lake, "dup") ~ "duplicate",
      TRUE ~ "regular"
    )
  )

saveRDS(
  phyto_cleanwdups,
  "~/Desktop/Project/Brooks_lake_2025/data_clean/phytoplankton/phyto_clean_wdups.rds"
)