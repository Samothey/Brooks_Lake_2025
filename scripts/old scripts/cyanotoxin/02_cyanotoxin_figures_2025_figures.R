# 2025 Cyanotoxin Figures -------------------------------------------------
# Author: Samantha Peña
#
# Purpose ----------------------------------------------------------------
#   - Load cleaned + summarized tables produced by 01_cyanotoxin_pipeline_2025.R
#   - Generate core figures for thesis questions:
#       1) Shore vs buoy (GRAB concentrations)
#       2) Surface vs bottom (GRAB buoy vertical structure)
#       3) How patterns vary by lake
#       4) SPATT vs GRAB detection agreement
#       5) Detection frequency (GRAB + SPATT)
#
# Inputs (from pipeline) -------------------------------------------------
#   - grab_class_totals.rds
#   - spatt_presence.rds
#   - master_sample_tracking.rds (optional)
#
# Notes ------------------------------------------------------------------
#   - Keep this script "plot-only": no data cleaning, no Excel reading.
#   - Filter/scope is done via small helper functions below.

# ---- Packages ----
library(tidyverse)
library(lubridate)

# ---- Paths ----
out_dir <- "C:/Users/spena5/Desktop/Project/Brooks_lake_2025/CODE/cyanotoxin/Cyanotoxins_2025"

grab_path  <- file.path(out_dir, "grab_class_totals.rds")
spatt_path <- file.path(out_dir, "spatt_presence.rds")

# Optional: if you want status info later
# status_path <- file.path(out_dir, "master_sample_tracking.rds")

# ---- Load tables ----
grab_class_totals <- readRDS(grab_path)
spatt_presence    <- readRDS(spatt_path)
# master_sample_tracking <- readRDS(status_path)

# ---- Constants ----
mc_threshold <- 0.8  # Make sure subtitle matches this number
site_type_levels <- c("shore", "buoy")
depth_levels <- c("surface", "bottom")

#
grab_mc <- grab_class_totals %>%
  filter(
    toxin_class == "microcystin",
    !is.na(lake), lake != "blank", lake != "boysen",
    site_type %in% site_type_levels
  ) %>%
  mutate(
    lake = factor(lake),
    site_type = factor(site_type, levels = site_type_levels),
    depth_category = factor(depth_category, levels = depth_levels)
  )

spatt_mc <- spatt_presence %>%
  filter(
    toxin_class == "microcystin",
    !is.na(lake), lake != "blank", lake != "boysen",
    site_type %in% site_type_levels
  ) %>%
  mutate(
    lake = factor(lake),
    site_type = factor(site_type, levels = site_type_levels),
    depth_category = factor(depth_category, levels = depth_levels)
  )


# Graphs and plots --------------------------------------------------------



# FIGURE 1: Violin (MC by lake)  ---------------------------------


v_1 <- grab_mc %>%
  mutate(value_log1p = log1p(total)) %>%
  filter(!is.na(value_log1p))

n_labs_v1 <- v_1 %>%
  group_by(lake) %>%
  summarise(n = n(), .groups = "drop")

y_top_v1 <- max(v_1$value_log1p, na.rm = TRUE) * 1.2

p_violin_lake <- ggplot(v_1, aes(x = lake, y = value_log1p, fill = lake)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", alpha = 0.8, width = 0.8) +
  geom_text(
    data = n_labs_v1,
    aes(x = lake, y = y_top_v1, label = paste0("n=", n)),
    inherit.aes = FALSE,
    vjust = -0.6,
    size = 3
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Lake",
    y = "log(1 + microcystin concentration)",
    title = "Distribution of microcystin concentrations by lake (GRAB; shore+buoy combined)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 5.5, 18, 5.5)
  )

p_violin_lake



# FIGURE 2: Violin( MC by lake, shore vs buoy) ----------------------------


df_v2 <- grab_mc %>%
  filter(site_type %in% c("shore", "buoy")) %>%
  mutate(value_log1p = log1p(total)) %>%
  filter(!is.na(value_log1p)) 

n_labs_v2 <- df_v2 %>%
  group_by(lake, site_type) %>%
  summarise(n = n(), .groups = "drop")

y_top_v2 <- max(df_v2$value_log1p, na.rm = TRUE) * 1.2

p_violin_shore_buoy <- ggplot(df_v2, aes(x = lake, y = value_log1p, fill = site_type)) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "black",
    alpha = 0.85,
    position = position_dodge(width = 0.85)
  ) +
  geom_text(
    data = n_labs_v2,
    aes(x = lake, y = y_top_v2, label = paste0("n=", n), group = site_type),
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
    title = "Microcystin distributions by lake: shore vs buoy (GRAB)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 5.5, 18, 5.5)
  )

p_violin_shore_buoy



# FIGURE 3: Boxplot (MC shore vs buoy) ------------------------------------

p_box_shore_buoy <- ggplot(grab_mc, aes(x = site_type, y = total)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, fill = "grey80") +
  geom_jitter(aes(color = site_type), width = 0.15, alpha = 0.7, size = 2) +
  geom_hline(aes(yintercept = mc_threshold, color = "Microcystin threshold"),
linetype = "dashed", 
linewidth = 0.8
) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  scale_color_manual(
    values = c(
      "shore" = "salmon",
      "buoy"  = "green4",
      "Microcystin threshold" = "red"
    ),
    name = ""   # removes "Site type" title
  ) + labs(
    x = "Site type",
    y = "Total toxin (µg/L)",
    title = "GRAB samples: shore vs buoy",
  ) + theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

p_box_shore_buoy



# FIGURE 4: Time series — MC over time ------------------------------------
# (shore vs buoy; per-site lines)
# i need to make Brooks shore lines different colors 

grab_mc_ts <- grab_mc %>%
  mutate(
    sample_date = as.Date(date),
    depth_plot = case_when(
      site_type == "shore" ~ "Shore",
      site_type == "buoy" & depth_category == "surface" ~ "Buoy surface",
      site_type == "buoy" & depth_category == "bottom"  ~ "Buoy bottom",
      TRUE ~ NA_character_
    ),
    line_id = if_else(site_type == "buoy", paste(site_id, depth_plot, sep = "_"), site_id)
  ) %>%
  filter(!is.na(depth_plot))

p_ts <- ggplot(grab_mc_ts, aes(x = sample_date, y = total, group = line_id, color = depth_plot)) +
  geom_line(alpha = 0.8) +
  geom_point(size = 2, alpha = 0.85) +
  geom_hline(yintercept = mc_threshold, color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Sample date",
    y = "Microcystins (µg/L)",
    color = "",
    title = "Microcystins over time by lake (GRAB)",
    subtitle = "Lines represent individual sites; buoy split into surface vs bottom"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

p_ts

## Figure 5: "Paired GRAB comparisons: shore vs buoy (paired by lake + date)" 
#( i want to add in legend to this graph for each site, but also simplify legend/site_ID) 

grab_paired <- grab_mc %>%
  group_by(lake, date, toxin_class, site_type) %>%
  summarise(total = mean(total, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = site_type, values_from = total) %>%
  filter(!is.na(shore), !is.na(buoy)) %>%
  mutate(pair_id = factor(paste(lake, date, toxin_class, sep = "_"))) %>%
  pivot_longer(cols = c("shore", "buoy"), names_to = "site_type", values_to = "total") %>%
  mutate(site_type = factor(site_type, levels = site_type_levels))

p_slope_shore_buoy <- ggplot(grab_paired, aes(x = site_type, y = total, group = pair_id, color = pair_id)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2, alpha = 0.9) +
  facet_grid(toxin_class ~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = NULL,
    y = "Microcystin total (log1p)",
    title = "Paired GRAB comparisons: shore vs buoy (paired by lake + date)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p_slope_shore_buoy


# FIGURE 6: Boxplot, buoy (surface vs bottom)  ----------------------------
#need to add legend to this box plot 
grab_buoy_depth <- grab_mc %>%
  filter(site_type == "buoy", depth_category %in% c("surface", "bottom"))

p_depth_box <- ggplot(grab_buoy_depth, aes(x = depth_category, y = total)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    x = "Buoy depth",
    y = "Microcystin total (log1p)",
    title = "GRAB buoy samples: surface vs bottom"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p_depth_box


#  FIGURE 7 ---------------------------------------------------------------




# FIGURE 8: Agreement / mismatch — SPATT vs GRAB (microcystin) ------------

# Make detection tables comparable at the same grain:
# (lake, site_id, date, toxin_class) with method column

spatt_detect <- spatt_presence %>%
  scope_real_lakes(drop_boysen = TRUE) %>%
  scope_microcystin() %>%
  transmute(lake, site_id, date, toxin_class, method = "spatt", detected = as.logical(detected))

grab_detect <- grab_class_totals %>%
  scope_real_lakes(drop_boysen = TRUE) %>%
  scope_microcystin() %>%
  transmute(lake, site_id, date, toxin_class, method = "grab", detected = !is.na(total) & total > 0)

detect_compare <- bind_rows(grab_detect, spatt_detect) %>%
  distinct(lake, site_id, date, toxin_class, method, detected) %>%
  pivot_wider(names_from = method, values_from = detected) %>%
  filter(!is.na(grab), !is.na(spatt)) %>%
  mutate(
    mismatch_type = case_when(
      spatt & !grab ~ "SPATT detected, GRAB not",
      grab & !spatt ~ "GRAB detected, SPATT not",
      spatt & grab  ~ "Both detected",
      TRUE          ~ "Both not detected"
    )
  )

detect_compare_long <- detect_compare %>%
  pivot_longer(cols = c(grab, spatt), names_to = "method", values_to = "detected") %>%
  mutate(
    method = factor(method, levels = c("grab", "spatt")),
    detected_num = as.integer(detected)
  )

p_mismatch <- ggplot(detect_compare_long, aes(x = date, y = detected_num)) +
  geom_jitter(aes(shape = method, color = mismatch_type),
              height = 0.2, width = 0.25, alpha = 0.85, size = 2.6) +
  facet_grid(toxin_class ~ lake) +
  scale_y_continuous(
    breaks = c(0, 1),
    labels = c("Not detected", "Detected"),
    limits = c(-0.15, 1.15)
  ) +
  labs(
    x = "Date",
    y = "Toxin detection",
    color = "Agreement",
    shape = "Method",
    title = "Microcystin detection agreement: GRAB vs SPATT"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

p_mismatch

# Optional: mismatch counts table
detect_compare %>%
  count(mismatch_type) %>%
  arrange(desc(n)) %>%
  print()




         