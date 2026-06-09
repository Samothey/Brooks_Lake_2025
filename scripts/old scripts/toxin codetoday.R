## toxin code library(tidyverse)
library(lubridate)

# ---- Paths ----
out_dir <- "~/Desktop/Project/Brooks_lake_2025/data_clean/toxins"

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
