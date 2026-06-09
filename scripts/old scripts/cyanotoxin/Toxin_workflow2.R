################################################
#   2025 Cyanotoxin Data 
#   Author: Samantha Peña
#   Summary: This script contains cyanotoxin (SPATT & grab) data from the Brooks lake watershed and Boysen reservoir
#  
#   Purpose:
#      - upload cyanotoxin data (spatt, grab and sample tracker spreadsheet)
#      - 
#      -  
#.     - 
###############################################################################============= 
library(tidyverse)
library(readxl)
library(janitor)
library(stringr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(colorspace)



main_toxin <- ("~/Desktop/Project/Brooks_lake_2025/CODE/cyanotoxin/Cyanotoxins_2025/Pena_2025.xlsx")

## exclude extra samples
exclude_site_patterns <- c("^BKS_UB_NCOVE")  

filter_extras <- function(df) {
  df %>%
    filter(!if_any(site_id, ~ str_detect(.x, exclude_site_patterns)))
}

#  derived site fields, this function will ensure the site IDs are correlated with the correct lakes, site type, and depth)
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
        str_detect(site_id, "^SPATT_")   ~ "blank",
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


# ---- Sheet 1: Whole water results ----
grab_result <- read_excel(main_toxin, sheet = "WholeWater") %>%
  clean_names() %>%
  mutate(sample_number = str_pad(as.character(sample_number), width = 3, pad = "0")) %>%
  rename(date = sample_date) %>%            
  mutate(
    date = as.Date(date),
    sample_type = str_to_lower(sample_type),
    method = "grab"
  ) %>%
  filter_extras() %>%    
  add_site_fields()

# ---- Sheet 2: SPATT results (only recovered, analyzed + shipped bags will be here) ----
spatt_result <- read_excel(main_toxin, sheet = "SPATT") %>%
  clean_names() %>%
  mutate(sample_number = str_pad(as.character(sample_number), width = 3, pad = "0")) %>%
  rename(date = sample_date) %>%            
  mutate(
    date = as.Date(date),
    sample_type = str_to_lower(sample_type),
    method = "spatt"
  ) %>%
  filter_extras() %>%    
  add_site_fields()

#------- Sheet 3: All sample status( samples could either be: collected, shipped, lost)
sample_status <- read_excel(main_toxin, sheet = "Sample Status") %>%
  clean_names() %>%
  mutate(sample_number = str_pad(as.character(sample_number), width = 3, pad = "0")) %>%
  rename(date = sample_date) %>%
  mutate(
    date = as.Date(date),
    sample_type = str_to_lower(sample_type),
    status = str_to_lower(status),
    method = str_to_lower(method)
  ) %>%
  filter_extras() %>%    
  add_site_fields()


######################################################################################
############################ building a master sheet  ################################
######################################################################################

#----- A)
sample_status <- read_excel(main_toxin, sheet = "Sample Status") %>%
  clean_names() %>%
  mutate(sample_number = str_pad(as.character(sample_number), width = 3, pad = "0")) %>%
  rename(date = sample_date) %>%
  mutate(
    date = as.Date(date),
    sample_type = str_to_lower(sample_type),
    status = str_to_lower(status),
    method = str_to_lower(method)
  ) %>%
  filter_extras() %>%    
  add_site_fields()

# ---- B) Joining the results "keys", this doesnt take the toxin result, it only creates "keys" for if the sample was analyzed
grab_keys <- grab_result %>%
  distinct(sample_number, site_id, date, sample_type, method) %>%
  mutate(has_grab_result = TRUE)

spatt_keys <- spatt_result %>%
  distinct(sample_number, site_id, date, sample_type, method) %>%
  mutate(has_spatt_result = TRUE)

# ---- C) apply rules, this will help us stay organized and understand which samples (analyzed overrides; shipped + no results -> missing)  ----
sample_status_Final <- sample_status %>%
  left_join(grab_keys,  by = c("sample_number", "site_id", "date", "sample_type", "method")) %>%
  left_join(spatt_keys, by = c("sample_number", "site_id", "date", "sample_type", "method")) %>%
  mutate(
    has_grab_result  = coalesce(has_grab_result, FALSE),
    has_spatt_result = coalesce(has_spatt_result, FALSE),
    status_final = case_when(
      method == "grab"  & has_grab_result  ~ "analyzed",
      method == "spatt" & has_spatt_result ~ "analyzed",
      status == "shipped" & method == "grab"  & !has_grab_result  ~ "missing",
      status == "shipped" & method == "spatt" & !has_spatt_result ~ "missing",
      TRUE ~ status
    )
  )


## Summary of what we have
status_summary_key <- sample_status_Final %>%
  mutate(
    method = str_to_lower(method),
    status_final = str_to_lower(status_final)
  ) %>%
  filter(status_final %in% c("missing", "lost", "analyzed")) %>%
  count(method, status_final) %>%
  pivot_wider(names_from = status_final, values_from = n, values_fill = 0)

status_summary_key





#### Toxin results table wide format- all (including blanks and dups)#####
toxin_results_wide_all <- bind_rows(grab_result, spatt_result) %>%
  mutate(
    method = str_to_lower(method),
    sample_type = str_to_lower(sample_type)
  )




####################################################################################################
##  create tables with toxin results- Toxin results long (regular samples only, no dups or blanks)

# identify toxin columns from the combined wide table
meta_cols <- c("sample_number","site_id","date","sample_type","method",
               "lake","site_type","depth_category")

toxin_cols <- toxin_results_wide_all %>%
  select(where(is.numeric)) %>%
  names() %>%
  setdiff(meta_cols)

toxin_results_long_regular <- toxin_results_wide_all %>%
  filter(sample_type == "regular", !lake %in% c("blank", "duplicate"), !is.na(lake)) %>%
  pivot_longer(
    cols = all_of(toxin_cols),
    names_to = "congener",
    values_to = "value"
  )

## plotting the tables

# 1) Identify toxin (numeric) columns from wide table
meta_cols <- c("sample_number","site_id","date","sample_type","method",
               "lake","site_type","depth_category")

toxin_cols <- toxin_results_wide_all %>%
  select(where(is.numeric)) %>%
  names() %>%
  setdiff(meta_cols)

# 2) Filter to real samples and pivot long
toxin_results_long_regular <- toxin_results_wide_all %>%
  filter_extras() %>%
  mutate(
    method = str_to_lower(method),
    sample_type = str_squish(str_to_lower(sample_type)),
    lake = str_to_lower(lake),
    site_type = str_to_lower(site_type)
  ) %>%
  filter(
    sample_type == "regular",
    !is.na(lake), lake != "blank"
  ) %>%
  pivot_longer(
    cols = all_of(toxin_cols),
    names_to = "congener",
    values_to = "value"
  ) %>%
  mutate(
    congener = str_to_lower(congener)
  )
## 2.) create toxin- class key + join it

toxin_key <- tribble(
  ~congener, ~toxin_class,
  "dm_lr", "microcystin",
  "la",   "microcystin",
  "lf",   "microcystin",
  "lr",   "microcystin",
  "ly",   "microcystin",
  "rr",   "microcystin",
  "wr",   "microcystin",
  "yr",   "microcystin",
  
  "atx",  "anatoxin",
  "hatx", "anatoxin",
  
  "cyl",  "cylindrospermopsin",
  "nod",  "nodularin"
)

toxin_plot <- toxin_results_long_regular %>%
  left_join(toxin_key, by = "congener")


# GRAB totals by toxin class
grab_class_totals <- toxin_plot %>%
  filter_extras() %>%
  filter(method == "grab") %>%
  group_by(
    sample_number, site_id, date, lake,
    site_type, depth_category, method, toxin_class
  ) %>%
  summarise(
    total = sum(value, na.rm = TRUE),
    .groups = "drop"
  )

## SPATT presence/absence by toxin class
spatt_presence <- toxin_plot %>%
  filter(method == "spatt") %>%
  group_by(
    sample_number, site_id, date, lake,
    site_type, depth_category, method, toxin_class
  ) %>%
  summarise(
    detected = any(value > 0, na.rm = TRUE),
    .groups = "drop"
  )

############################################# PLOTS ###########################

############ grab sample shore vs. buoy #################

grab_class_totals_focus <- grab_class_totals %>%
  filter(
    lake != "boysen",
    toxin_class %in% c("microcystin")
  )
### violin plots 

# ---- prep data ----
df_plot <- grab_class_totals_focus %>%
  filter(!is.na(lake), lake != "blank") %>%
  filter(toxin_class == "microcystin") %>%
  filter(site_type %in% c("shore", "buoy")) %>%
  mutate(
    lake = factor(lake),
    value = total,      # <-- CHANGE if needed
    value_log1p = log1p(value)
  ) %>%
  filter(!is.na(value_log1p))

# ---- sample size labels ----
n_labs <- df_plot %>%
  group_by(lake) %>%
  summarise(n = n(), .groups = "drop")

y_top <- max(df_plot$value_log1p, na.rm = TRUE) * 1.2



# ---- violin plot ----
ggplot(df_plot, aes(x = lake, y = value_log1p, fill = lake)) +
  geom_violin(
    trim = FALSE,
    scale = "width",   # keeps violins comparable with unequal n
    color = "black",   # outline
    alpha = 0.8
  ) +
  geom_text(
    data = n_labs,
    aes(x = lake, y = y_top, label = paste0("n=", n)),
    inherit.aes = FALSE,
    vjust = -0.6,
    size = 3
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Lake",
    y = "log(1 + microcystin concentration)",
    title = "Distribution of microcystin concentrations by lake (shore and buoy samples combined)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 5.5, 18, 5.5)
  )

### violins for shore and buoy for each lake ###


# ---- prep data ----
df_plot <- grab_class_totals_focus %>%
  filter(!is.na(lake), lake != "blank", lake != "boysen") %>%
  filter(toxin_class == "microcystin") %>%
  filter(site_type %in% c("shore", "buoy")) %>%
  mutate(
    lake = factor(lake),
    site_type = factor(site_type, levels = c("shore", "buoy")),
    value = total,     
    value_log1p = log1p(value)
  ) %>%
  filter(!is.na(value_log1p))

# ---- sample sizes per lake x site_type ----
n_labs <- df_plot %>%
  group_by(lake, site_type) %>%
  summarise(n = n(), .groups = "drop")

# Put labels above the tallest violin + add padding
y_top <- max(df_plot$value_log1p, na.rm = TRUE) * 1.2

ggplot(df_plot, aes(x = lake, y = value_log1p, fill = site_type)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "black",
    alpha = 0.85,
    position = position_dodge(width = 0.85)
  ) +
  geom_text(
    data = n_labs,
    aes(x = lake, y = y_top, label = paste0("n=", n), group = site_type),
    position = position_dodge(width = 0.85),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Lake",
    y = "log(1 + microcystin concentration)",
    fill = "Site type",
    title = "Microcystin distributions by lake: shore vs buoy"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 5.5, 18, 5.5)
  )





#### box plot
mc_threshold <- 0.8  # µg/L (EPA recreational guidance)

ggplot(grab_class_totals_focus,
       aes(x = site_type, y = total)) +

  geom_boxplot(outlier.shape = NA, alpha = 0.4, fill = "grey80") +
  geom_jitter(
    aes(color = site_type),
    width = 0.15,
    alpha = 0.7,
    size = 2
  ) +
  geom_hline(
    aes(yintercept = mc_threshold, color = "Microcystin threshold"),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_grid(toxin_class ~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  scale_color_manual(
    values = c(
      "shore" = "#1f78b4",
      "buoy"  = "purple4",
      "Microcystin threshold" = "red"
    ),
    name = ""   # removes "Site type" title
  ) + labs(
    x = "Site type",
    y = "Total toxin (µg/L)",
    title = "GRAB samples: shore vs buoy",
    subtitle = "Dashed red line = EPA recreational microcystin threshold (8 µg/L)"
  ) + theme(
    legend.position = "right",
          panel.grid = element_blank()
        )





##### seperate graph for each lake #######



grab_mc_ts <- grab_class_totals %>%
  filter(lake != "boysen", lake != "upper brooks", toxin_class == "microcystin") %>%
  mutate(
    sample_date = as.Date(date),
    site_id = as.character(site_id),
    depth_category = as.character(depth_category),
    
    depth_plot = case_when(
      site_type == "shore" ~ "Shore",
      site_type == "buoy" & depth_category %in% c("surface") ~ "Buoy surface",
      site_type == "buoy" ~ "Buoy depth",
      TRUE ~ NA_character_
    ),
    
    line_id = if_else(
      site_type == "buoy",
      paste(site_id, depth_plot, sep = "_"),
      site_id
    )
  ) %>%
  filter(!is.na(depth_plot))



### graphs ( rainbow,, lower jade) 
# lakes to include in panel A
lakes_A <- c("rainbow", "lower jade")

p_A <- grab_mc_ts %>%
  filter(lake %in% lakes_A) %>%
  ggplot(aes(x = sample_date, y = total, group = line_id, color = depth_plot)) +
  geom_line(alpha = 0.8) +
  geom_point(size = 2, alpha = 0.85) +
  geom_hline(yintercept = mc_threshold, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  scale_color_manual(
    values = c(
      "Shore"        = "#0072B2",
      "Buoy surface" = "#E69F00",
      "Buoy depth"   = "#009E73"
    ),
    name = ""
  ) +
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    title = "Microcystins over time (Rainbow, Lower Jade)",
    subtitle = "Lines represent individual sites; colors indicate sampling context"
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

p_A

 ## this code doesnt work ?
lakes_b <- c("brook")

p_B <- grab_mc_ts %>%
  filter(lake %in% lakes_b) %>%
  ggplot(aes(x = sample_date, y = total, group = line_id, color = depth_plot)) +
  geom_line(alpha = 0.8) +
  geom_point(size = 2, alpha = 0.85) +
  geom_hline(yintercept = mc_threshold, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  scale_color_manual(
    values = c(
      "Shore"        = "#0072B2",
      "Buoy surface" = "#E69F00",
      "Buoy depth"   = "#009E73"
    ),
    name = ""
  ) +
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    title = "Microcystins over time (Upper brooks)",
    subtitle = "Lines represent individual sites; colors indicate sampling context"
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )
p_B




lake_target <- "brooks"   

# 1) Shore only
p_ub_shore <- grab_mc_ts %>%
  filter(lake == lake_target, site_type == "shore") %>%
  ggplot(aes(x = sample_date, y = total, group = site_id, color = site_id)) +
  geom_line(alpha = 0.85) +
  geom_point(size = 2, alpha = 0.85) +
  geom_hline(yintercept = mc_threshold, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    color = "Shore site ID",
    title = "Upper Brooks: Shore sites only",
    subtitle = "Dashed red line = EPA recreational microcystin threshold (8 µg/L)"
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )
p_ub_shore

# 2) Buoy only (surface vs depth shown by shape)
p_ub_buoy <- grab_mc_ts %>%
  filter(lake == lake_target, site_type == "buoy") %>%
  ggplot(aes(x = sample_date, y = total, group = line_id, color = site_id)) +
  geom_line(alpha = 0.85) +
  geom_point(aes(shape = depth_plot), size = 2, alpha = 0.85) +
  geom_hline(yintercept = mc_threshold, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    color = "Buoy site ID",
    shape = "Buoy depth",
    title = "Upper Brooks: Buoy sites only",
    subtitle = "Surface vs depth shown by point shape; dashed red line = 8 µg/L"
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )
p_ub_buoy

# 3) All together (shore + buoy)
p_ub_all <- grab_mc_ts %>%
  filter(lake == lake_target) %>%
  ggplot(aes(x = sample_date, y = total, group = line_id, color = site_id)) +
  geom_line(alpha = 0.85) +
  geom_point(aes(shape = depth_plot), size = 2, alpha = 0.85) +
  geom_hline(yintercept = mc_threshold, color = "red",
             linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    color = "Site ID",
    shape = "Sample type",
    title = "Upper Brooks: Shore + buoy sites",
    subtitle = "Buoy split into surface vs depth; dashed red line = 8 µg/L"
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

# Print the plots
p_ub_shore
p_ub_buoy
p_ub_all


# factor level ordering (define once)
site_type_levels <- c("shore", "buoy")
depth_levels <- c("surface", "bottom")

grab_class_totals2 <- grab_class_totals_focus %>%
  filter(!is.na(lake), lake != "blank", lake != "boysen") %>%
  filter(site_type %in% site_type_levels) %>%
  mutate(
    site_type = factor(site_type, levels = site_type_levels),
    depth_category = factor(depth_category, levels = depth_levels)
  )

spatt_presence2 <- spatt_presence %>%
  filter(!is.na(lake), lake != "blank", lake != "boysen") %>%
  filter(site_type %in% site_type_levels) %>%
  mutate(
    site_type = factor(site_type, levels = site_type_levels),
    depth_category = factor(depth_category, levels = depth_levels)
  )

## paired shore buoy plots
grab_paired_shore_buoy <- grab_class_totals2 %>%
  filter(method == "grab") %>%
  group_by(lake, date, toxin_class, site_type) %>%
  summarise(total = mean(total, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = site_type, values_from = total) %>%
  filter(!is.na(shore), !is.na(buoy)) %>%
  mutate(
    pair_id = paste(lake, date, toxin_class, sep = "_"),
    pair_id = factor(pair_id)  # ensures discrete legend
  )

# slope plots.  ***** i like this graph 
grab_paired_long <- grab_paired_shore_buoy %>%
  pivot_longer(cols = c("shore", "buoy"), names_to = "site_type", values_to = "total") %>%
  mutate(
    site_type = factor(site_type, levels = site_type_levels),
    pair_id = factor(pair_id)  # keep consistent after pivot
  )

ggplot(
  grab_paired_long,
  aes(x = site_type, y = total, group = pair_id, color = pair_id)
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2, alpha = 0.9) +
  facet_grid(toxin_class ~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = NULL,
    y = "Total toxin (log1p)",
    color = "Site ID (pair_id)",
    title = "Paired GRAB comparisons: shore vs buoy (paired by lake + date)"
  )




###surface vs bottom (Grab)


# -----------------------------
# buoy-only data
# -----------------------------
grab_buoy_depth <- grab_class_totals2 %>%
  filter(method == "grab", site_type == "buoy") %>%
  filter(depth_category %in% c("surface", "bottom"))

# -----------------------------
# exploratory box + jitter plot (unchanged)
# -----------------------------
ggplot(grab_buoy_depth,
       aes(x = depth_category, y = total)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_grid(toxin_class ~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Buoy depth",
    y = "Total toxin (log1p)",
    title = "GRAB buoy samples: surface vs bottom"
  ) +
  theme(legend.position = "none")

# -----------------------------
# paired surface–bottom by date
# -----------------------------
grab_depth_paired <- grab_buoy_depth %>%
  group_by(lake, date, toxin_class, depth_category) %>%
  summarise(total = mean(total, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = depth_category, values_from = total) %>%
  filter(!is.na(surface), !is.na(bottom)) %>%
  mutate(
    pair_id = paste(lake, date, toxin_class, sep = "_"),
    pair_id = factor(pair_id)          # ensure discrete colors
  ) %>%
  pivot_longer(
    cols = c(surface, bottom),
    names_to = "depth_category",
    values_to = "total"
  ) %>%
  mutate(
    depth_category = factor(depth_category, levels = depth_levels)
  )

# -----------------------------
# slope plot: colored by pair_id
# -----------------------------
ggplot(
  grab_depth_paired,
  aes(
    x = depth_category,
    y = total,
    group = pair_id,
    color = pair_id
  )
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2, alpha = 0.9) +
  facet_grid(toxin_class ~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = NULL,
    y = "Total toxin (log1p)",
    color = "Site ID (pair_id)",
    title = "Paired GRAB buoy comparisons: surface vs bottom (paired by date)"
  )




## detection frequency plots
## grab

grab_detection <- grab_class_totals2 %>%
  filter(method == "grab") %>%
  mutate(detected = total > 0) %>%
  group_by(lake, site_type, toxin_class) %>%
  summarise(
    n = n(),
    n_detected = sum(detected),
    detection_rate = n_detected / n,
    .groups = "drop"
  )


## spatt frequency
spatt_detection <- spatt_presence2 %>%
  group_by(lake, site_type, toxin_class) %>%
  summarise(
    n = n(),
    n_detected = sum(detected),
    detection_rate = n_detected / n,
    .groups = "drop"
  )

## time series
## grab

grab_ts <- grab_class_totals2 %>%
  filter(method == "grab") %>%
  filter(toxin_class %in% c("microcystin"))

ggplot(grab_ts,
       aes(x = date, y = total, color = site_type)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(group = interaction(site_id, site_type)), alpha = 0.25) +
  facet_grid(toxin_class ~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Date",
    y = "Total toxin (log1p)",
    title = "GRAB time series by lake"
  )

## spatt
spatt_ts <- spatt_presence2 %>%
  filter(toxin_class %in% c("microcystin")) %>%
  mutate(detected_num = as.integer(detected))

ggplot(
  spatt_ts,
  aes(x = date, y = detected_num, color = site_type)
) +
  geom_jitter(
    height = 0.06,
    width = 0,
    alpha = 0.75,
    size = 2.2
  ) +
  facet_grid(toxin_class ~ lake) +
  scale_y_continuous(
    breaks = c(0, 1),
    labels = c("Not detected", "Detected"),
    limits = c(-0.15, 1.15)
  ) +
  labs(
    x = "Date",
    y = "SPATT detection",
    title = "SPATT microcystin detections through time"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )



# ---- SPATT detection (already boolean-ish) ----
spatt_detect <- spatt_presence2 %>%
  filter(toxin_class %in% c("microcystin")) %>%
  mutate(
    method = "spatt",
    detected = as.logical(detected),              # TRUE/FALSE
    detected_num = as.integer(detected)           # 1/0
  ) %>%
  select(lake, site_id, date, toxin_class, method, detected, detected_num)

# ---- GRAB detection from concentration ----
grab_detect <- grab_class_totals2 %>%
  filter(method == "grab", toxin_class %in% c("microcystin")) %>%
  mutate(
    method = "grab",
    detected = !is.na(total) & total > 0,         # RULE: detected if total > 0
    detected_num = as.integer(detected)
  ) %>%
  select(lake, site_id, date, toxin_class, method, detected, detected_num)

# ---- Combine into one dataframe ----
detect_ts <- bind_rows(grab_detect, spatt_detect) %>%
  filter(!is.na(lake), lake != "blank", lake != "boysen") %>%
  mutate(
    method = factor(method, levels = c("grab", "spatt"))
  )




### going to plot specifically showing when they dont match #####

grab_detect <- grab_class_totals2 %>%
  filter(method == "grab", toxin_class == "microcystin") %>%
  mutate(
    method = "grab",
    detected = !is.na(total) & total > 0
  ) %>%
  select(lake, site_id, date, toxin_class, method, detected)

spatt_detect <- spatt_presence2 %>%
  filter(toxin_class == "microcystin") %>%
  mutate(
    method = "spatt",
    detected = as.logical(detected)
  ) %>%
  select(lake, site_id, date, toxin_class, method, detected)



detect_compare <- bind_rows(grab_detect, spatt_detect) %>%
  filter(!is.na(lake), lake != "blank", lake != "boysen") %>%
  distinct(lake, site_id, date, toxin_class, method, detected) %>%  # prevent duplicates
  pivot_wider(names_from = method, values_from = detected) %>%
  # keep only dates where you have both methods (true head-to-head comparison)
  filter(!is.na(grab), !is.na(spatt)) %>%
  mutate(
    mismatch_type = case_when(
      spatt & !grab ~ "SPATT detected, GRAB not",
      grab & !spatt ~ "GRAB detected, SPATT not",
      spatt & grab  ~ "Both detected",
      TRUE          ~ "Both not detected"
    ),
    mismatch_flag = mismatch_type %in% c("SPATT detected, GRAB not", "GRAB detected, SPATT not")
  )

detect_compare_long <- detect_compare %>%
  pivot_longer(cols = c(grab, spatt),
               names_to = "method",
               values_to = "detected") %>%
  mutate(
    method = factor(method, levels = c("grab", "spatt")),
    detected_num = as.integer(detected)
  )

ggplot(
  detect_compare_long,
  aes(x = date, y = detected_num)
) +
  geom_jitter(
    aes(shape = method, color = mismatch_type),
    height = 0.2, width = 0.25,
    alpha = 0.85, size = 2.6
  ) +
  facet_grid(toxin_class ~ lake) +
  scale_y_continuous(
    breaks = c(0, 1),
    labels = c("Not detected", "Detected"),
    limits = c(-0.15, 1.15)
  ) +
  # Make mismatches stand out; matches fade
  scale_color_manual(values = c(
    "SPATT detected, GRAB not" = "red3",
    "GRAB detected, SPATT not" = "dodgerblue3",
    "Both detected"            = "black",
    "Both not detected"        = "grey40"
  )) +
  labs(
    x = "Date",
    y = "Toxin detection",
    color = "Agreement",
    shape = "Method",
    title = "Microcystin detection agreement: GRAB vs SPATT (mismatches highlighted)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

detect_compare %>%
  count(mismatch_type) %>%
  filter(mismatch_type %in% c(
    "SPATT detected, GRAB not",
    "GRAB detected, SPATT not"
  ))

###
detect_compare %>%
  count(mismatch_type) %>%
  filter(mismatch_type %in% c(
    "SPATT detected, GRAB not",
    "GRAB detected, SPATT not"
  ))

detect_compare %>%
  filter(mismatch_type %in% c(
    "SPATT detected, GRAB not",
    "GRAB detected, SPATT not"
  )) %>%
  count(lake, mismatch_type) %>%
  tidyr::pivot_wider(
    names_from = mismatch_type,
    values_from = n,
    values_fill = 0
  )

detect_compare %>%
  filter(mismatch_type %in% c(
    "SPATT detected, GRAB not",
    "GRAB detected, SPATT not"
  )) %>%
  count(mismatch_type) %>%
  mutate(percent = 100 * n / sum(n))
detect_compare %>%
  count(mismatch_type) %>%
  mutate(percent = round(100 * n / sum(n), 1))



























