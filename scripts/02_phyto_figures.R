



library(tidyverse)
library(scales)

# Read data -------------------------------------------------------------

phyto_clean <- readRDS(
  "~/Desktop/Project/Brooks_lake_2025/data_clean/phytoplankton/phyto_clean.rds"
)

# Output folder ---------------------------------------------------------

fig_dir <- "~/Desktop/Project/Brooks_lake_2025/figures/phytoplankton"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ======================================================================
# 1. Total cyanobacteria abundance
# ======================================================================

cyano_total <- phyto_clean %>%
  filter(division == "cyanophyta") %>%
  group_by(lake, date) %>%
  summarise(
    cyano_cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  )

# Raw abundance: actual units ------------------------------------------

p_cyano_raw <- ggplot(
  cyano_total,
  aes(x = date, y = cyano_cells)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(
    x = NULL,
    y = expression(paste("Cyanobacteria abundance (cells ", L^{-1}, ")"))
  ) +
  theme_bw()

p_cyano_raw

ggsave(
  filename = file.path(fig_dir, "cyano_total_raw_cells_per_L.png"),
  plot = p_cyano_raw,
  width = 10,
  height = 6,
  dpi = 300
)

# Log10 abundance: transformed -----------------------------------------

p_cyano_log <- ggplot(
  cyano_total,
  aes(x = date, y = cyano_cells)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ lake) +
  scale_y_log10(labels = comma) +
  labs(
    x = NULL,
    y = expression(paste("Log"[10], " cyanobacteria abundance (cells ", L^{-1}, ")"))
  ) +
  theme_bw()

p_cyano_log

ggsave(
  filename = file.path(fig_dir, "cyano_total_log10_cells_per_L.png"),
  plot = p_cyano_log,
  width = 10,
  height = 6,
  dpi = 300
)

# ======================================================================
# 2. Relative abundance by phytoplankton division
# ======================================================================

phyto_division <- phyto_clean %>%
  mutate(
    division = case_when(
      division == "bacillariophyta" ~ "Diatoms",
      division == "chlorophyta" ~ "Green algae",
      division == "chrysophyta" ~ "Chrysophytes",
      division == "cryptophyta" ~ "Cryptophytes",
      division == "cyanophyta" ~ "Cyanobacteria",
      division %in% c("pyrrhophyta", "rotifera") ~ "Other",
      TRUE ~ division
    )
  ) %>%
  group_by(lake, date, division) %>%
  summarise(
    cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(lake, date) %>%
  mutate(
    rel_abundance = cells / sum(cells, na.rm = TRUE)
  ) %>%
  ungroup()

phyto_division_rel <- ggplot(
  phyto_division,
  aes(x = date, y = rel_abundance, fill = division)
) +
  geom_col(width = 10, alpha = 0.9) +
  facet_wrap(~ lake) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Division"
  ) +
  theme_bw()

phyto_division_rel

ggsave(
  filename = file.path(fig_dir, "phyto_division_relative_abundance.png"),
  plot = phyto_division_rel,
  width = 10,
  height = 6,
  dpi = 300
)

# ======================================================================
# 3. Cyanobacteria vs other phytoplankton
# ======================================================================

phyto_cyano <- phyto_clean %>%
  mutate(
    group = if_else(
      division == "cyanophyta",
      "Cyanobacteria",
      "Other phytoplankton"
    )
  ) %>%
  group_by(lake, date, group) %>%
  summarise(
    cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(lake, date) %>%
  mutate(
    rel_abundance = cells / sum(cells, na.rm = TRUE)
  ) %>%
  ungroup()

p_cyano_vs_other_rel <- ggplot(
  phyto_cyano,
  aes(x = date, y = rel_abundance, fill = group)
) +
  geom_area(alpha = 0.9) +
  facet_wrap(~ lake) +
  scale_y_continuous(labels = percent) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = NULL
  ) +
  theme_bw()

p_cyano_vs_other_rel

ggsave(
  filename = file.path(fig_dir, "cyano_vs_other_relative_abundance.png"),
  plot = p_cyano_vs_other_rel,
  width = 10,
  height = 6,
  dpi = 300
)

# ======================================================================
# 4. Relative abundance of cyanobacteria taxa only
# ======================================================================

cyano_taxa <- phyto_clean %>%
  filter(division == "cyanophyta") %>%
  group_by(lake, date, taxon) %>%
  summarise(
    cells = sum(total_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(lake, date) %>%
  mutate(
    rel_abundance = cells / sum(cells, na.rm = TRUE)
  ) %>%
  ungroup()

p_cyano_taxa_rel <- ggplot(
  cyano_taxa,
  aes(x = date, y = rel_abundance, fill = taxon)
) +
  geom_col(width = 10, alpha = 0.9) +
  facet_wrap(~ lake) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Relative abundance of cyanobacteria",
    fill = "Cyanobacteria taxon"
  ) +
  theme_bw()

p_cyano_taxa_rel

ggsave(
  filename = file.path(fig_dir, "cyano_taxa_relative_abundance.png"),
  plot = p_cyano_taxa_rel,
  width = 10,
  height = 6,
  dpi = 300
)

# ======================================================================
# 5. Cyanobacteria taxa abundance through time
# ======================================================================

p_cyano_taxa_cells <- ggplot(
  cyano_taxa,
  aes(x = date, y = cells, color = taxon)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ lake, scales = "free_y") +
  scale_y_log10(labels = comma) +
  labs(
    x = NULL,
    y = expression(paste("Cyanobacteria abundance (cells ", L^{-1}, ")")),
    color = "Taxon"
  ) +
  theme_bw()

p_cyano_taxa_cells

ggsave(
  filename = file.path(fig_dir, "cyano_taxa_log10_cells_per_L.png"),
  plot = p_cyano_taxa_cells,
  width = 10,
  height = 6,
  dpi = 300
)

# ======================================================================
# 6. Helpful summary table for interpretation
# ======================================================================

cyano_summary <- cyano_total %>%
  group_by(lake) %>%
  summarise(
    peak_cyano_cells = max(cyano_cells, na.rm = TRUE),
    date_peak_cyano = date[which.max(cyano_cells)],
    mean_cyano_cells = mean(cyano_cells, na.rm = TRUE),
    .groups = "drop"
  )

cyano_summary

