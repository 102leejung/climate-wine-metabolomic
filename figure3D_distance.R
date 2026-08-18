## -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

file_path <- "figure2D_PCA.xlsx"
include_temps <- c(0, 1, 2, 3, 5)
sheet_order <- as.character(include_temps)

class_levels <- c("Australia", "Europe", "North America", "South America")

base_cols <- c(
  "Australia"     = "#E9c46A",
  "Europe"        = "#2a9d8f",
  "North America" = "#2e6f95",
  "South America" = "#d1495b"
)

raw_list <- lapply(sheet_order, function(sh) {
  df <- read_excel(file_path, sheet = sh)

  stopifnot(all(c("SampleID", "Class", "PC1", "PC2") %in% colnames(df)))

  df %>%
    mutate(
      Temperature = as.numeric(sh),
      Class = str_remove(as.character(Class), "\\+\\d+°C$"),
      Class = str_squish(Class)
    )
})

pca_df <- bind_rows(raw_list) %>%
  mutate(
    Temperature = factor(Temperature, levels = include_temps),
    Class = factor(Class, levels = class_levels)
  )

class_pal <- base_cols[class_levels]

centroid_df <- pca_df %>%
  group_by(Temperature, Class) %>%
  summarise(
    cPC1 = mean(PC1, na.rm = TRUE),
    cPC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Class = as.character(Class),
    Temperature = as.numeric(as.character(Temperature))
  )

dist_df <- centroid_df %>%
  group_by(Temperature) %>%
  group_modify(~{
    tmp <- .x %>% arrange(match(Class, class_levels))

    comb <- as.data.frame(t(combn(tmp$Class, 2)), stringsAsFactors = FALSE)
    colnames(comb) <- c("Class1", "Class2")

    out <- comb %>%
      rowwise() %>%
      mutate(
        x1 = tmp$cPC1[tmp$Class == Class1],
        y1 = tmp$cPC2[tmp$Class == Class1],
        x2 = tmp$cPC1[tmp$Class == Class2],
        y2 = tmp$cPC2[tmp$Class == Class2],
        dist = sqrt((x1 - x2)^2 + (y1 - y2)^2)
      ) %>%
      ungroup()

    mean_dist <- mean(out$dist, na.rm = TRUE)

    out %>%
      mutate(
        rel_dist = dist / mean_dist,
        label = sprintf("%.2f", rel_dist)
      ) %>%
      select(Class1, Class2, rel_dist, label)
  }) %>%
  ungroup()

heat_pairs <- bind_rows(
  dist_df %>%
    transmute(
      Temperature,
      Class1 = as.character(Class1),
      Class2 = as.character(Class2),
      rel_dist,
      label
    ),
  dist_df %>%
    rename(Class1_old = Class1, Class2_old = Class2) %>%
    transmute(
      Temperature,
      Class1 = as.character(Class2_old),
      Class2 = as.character(Class1_old),
      rel_dist,
      label
    )
) %>%
  distinct(Temperature, Class1, Class2, .keep_all = TRUE)

diag_df <- expand.grid(
  Temperature = include_temps,
  Class1 = class_levels,
  Class2 = class_levels,
  stringsAsFactors = FALSE
) %>%
  filter(Class1 == Class2) %>%
  mutate(
    rel_dist = NA_real_,
    label = NA_character_
  )

heat_df <- bind_rows(heat_pairs, diag_df) %>%
  distinct(Temperature, Class1, Class2, .keep_all = TRUE)

# -----------------------------
# 하삼각만 남기기 위한 index
# -----------------------------
heat_df <- heat_df %>%
  mutate(
    x_idx = match(Class1, class_levels),
    y_idx = match(Class2, class_levels)
  )

# 하삼각 + 대각선만 유지
# 대각선은 회색으로 남기고 싶으면 <=
# 대각선도 아예 없애려면 <
heat_df_lower <- heat_df %>%
  filter(x_idx <= y_idx)

heat_limits <- range(heat_df_lower$rel_dist, na.rm = TRUE)

make_pca_panel <- function(temp_value) {
  df_pca <- pca_df %>% filter(as.numeric(as.character(Temperature)) == temp_value)
  df_cent <- centroid_df %>% filter(Temperature == temp_value)

  ggplot(df_pca, aes(PC1, PC2, color = Class, fill = Class)) +
    geom_point(size = 2.3, alpha = 0.75) +
    stat_ellipse(
      geom = "polygon",
      type = "norm",
      level = 0.68,
      alpha = 0.15,
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    geom_point(
      data = df_cent,
      aes(cPC1, cPC2),
      inherit.aes = FALSE,
      shape = 4, size = 4, stroke = 1.4, color = "black"
    ) +
    geom_text(
      data = df_cent,
      aes(cPC1, cPC2, label = Class),
      inherit.aes = FALSE,
      size = 3, fontface = "bold", vjust = -0.8, color = "black"
    ) +
    scale_color_manual(values = class_pal) +
    scale_fill_manual(values = class_pal) +
    labs(title = paste0(temp_value, "°C"), x = "PC1", y = "PC2") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    coord_equal()
}

make_heat_panel <- function(temp_value) {
  df_heat <- heat_df_lower %>% filter(Temperature == temp_value)

  ggplot(df_heat, aes(x = Class1, y = Class2, fill = rel_dist)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(
      data = df_heat %>% filter(!is.na(label)),
      aes(label = label),
      size = 3.3,
      fontface = "bold"
    ) +
    scale_x_discrete(limits = class_levels) +
    scale_y_discrete(limits = rev(class_levels)) +
   scale_fill_gradientn(
  colors = c("#F7F7F7", "#BDBDBD", "#525252"),
  limits = heat_limits,
  oob = scales::squish,
  na.value = "grey92",
  name = "Relative\nDistance"
) +
    labs(title = paste0(temp_value, "°C"), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "none"
    ) +
    coord_fixed()
}

pca_panels  <- lapply(include_temps, make_pca_panel)
heat_panels <- lapply(include_temps, make_heat_panel)

pca_grid <- wrap_plots(pca_panels, ncol = 5) +
  plot_annotation(title = "PCA by temperature")

heat_grid <- wrap_plots(heat_panels, ncol = 5) +
  plot_annotation(
    title = "Relative inter-class distance by temperature",
    subtitle = "Lower triangle only; pairwise centroid distances normalized by the mean inter-class distance within each temperature."
  )

full_plot <- (pca_grid / heat_grid) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "right")

print(full_plot)

ggsave("PCA_5panels_and_heatmap_lowertriangle.png", full_plot, width = 22, height = 10, dpi = 300)
ggsave("PCA_5panels_and_heatmap_lowertriangle.pdf", full_plot, width = 22, height = 10)
ggsave("PCA_5panels_and_heatmap_lowertriangle.svg", full_plot, width = 22, height = 10)


## -----------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# =========================================================
# 0) Input
# =========================================================
file_path <- "figure2D_PCA.xlsx"
include_temps <- c(0, 1, 2, 3, 5)
sheet_order <- as.character(include_temps)

class_levels <- c("Australia", "Europe", "North America", "South America")

pair_levels <- c(
  "Australia vs Europe",
  "Australia vs North America",
  "Australia vs South America",
  "Europe vs North America",
  "Europe vs South America",
  "North America vs South America"
)

pair_colors <- c(
  "Australia vs Europe"        = "#984EA3",
  "Australia vs North America" = "#FF7F00",
  "Australia vs South America" = "#377EB8",
  "Europe vs North America"    = "#4DAF4A",
  "Europe vs South America"    = "#A65628",
  "North America vs South America" = "#F781BF"
)

# =========================================================
# 1) Read sheets
# =========================================================
raw_list <- lapply(sheet_order, function(sh) {
  df <- read_excel(file_path, sheet = sh)

  stopifnot(all(c("SampleID", "Class", "PC1", "PC2") %in% colnames(df)))

  df %>%
    mutate(
      Temperature = as.numeric(sh),
      Class = str_remove(as.character(Class), "\\+\\d+°C$"),
      Class = str_squish(Class)
    )
})

pca_df <- bind_rows(raw_list) %>%
  mutate(
    Temperature = as.numeric(Temperature),
    Class = factor(Class, levels = class_levels)
  )

# =========================================================
# 2) Centroids
# =========================================================
centroid_df <- pca_df %>%
  group_by(Temperature, Class) %>%
  summarise(
    cPC1 = mean(PC1, na.rm = TRUE),
    cPC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Class = as.character(Class)
  )

# =========================================================
# 3) Pairwise centroid distances by temperature
# =========================================================
centroid_split <- split(centroid_df, centroid_df$Temperature)

dist_list <- lapply(names(centroid_split), function(temp_name) {
  tmp <- centroid_split[[temp_name]]
  temp_value <- as.numeric(temp_name)

  tmp <- tmp %>%
    arrange(match(Class, class_levels))

  comb <- as.data.frame(
    t(combn(tmp$Class, 2)),
    stringsAsFactors = FALSE
  )
  colnames(comb) <- c("Class1", "Class2")

  out <- comb %>%
    rowwise() %>%
    mutate(
      x1 = tmp$cPC1[tmp$Class == Class1],
      y1 = tmp$cPC2[tmp$Class == Class1],
      x2 = tmp$cPC1[tmp$Class == Class2],
      y2 = tmp$cPC2[tmp$Class == Class2],
      dist = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    ) %>%
    ungroup()

  mean_dist <- mean(out$dist, na.rm = TRUE)

  out %>%
    mutate(
      Temperature = temp_value,
      rel_dist = dist / mean_dist,
      pair = paste(Class1, "vs", Class2)
    ) %>%
    select(Temperature, Class1, Class2, pair, dist, rel_dist)
})

dist_df <- bind_rows(dist_list) %>%
  mutate(
    pair = factor(pair, levels = pair_levels)
  )

print(dist_df)

# =========================================================
# 4) Main line plot
# =========================================================
p_line <- ggplot(
  dist_df,
  aes(x = Temperature, y = rel_dist, color = pair, group = pair)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_color_manual(values = pair_pal, drop = FALSE) +
  scale_x_continuous(breaks = include_temps) +
  labs(
    title = "Temperature-dependent change in relative inter-class distance",
    subtitle = "Relative distance = pairwise centroid distance normalized by the mean inter-class distance within each temperature",
    x = "Temperature (°C)",
    y = "Relative distance",
    color = "Class pair"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_line)

# =========================================================
# 5) Faceted version
# =========================================================
p_facet <- ggplot(
  dist_df,
  aes(x = Temperature, y = rel_dist, group = pair, color = pair)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_color_manual(values = pair_pal, drop = FALSE) +
  scale_x_continuous(breaks = include_temps) +
  facet_wrap(~ pair, ncol = 3, scales = "fixed") +
  labs(
    title = "Temperature-dependent change in relative inter-class distance",
    subtitle = "Each panel shows one class pair",
    x = "Temperature (°C)",
    y = "Relative distance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p_facet)

# =========================================================
# 6) Save
# =========================================================
ggsave(
  "pairwise_relative_distance_lineplot.png",
  p_line,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  "pairwise_relative_distance_lineplot.pdf",
  p_line,
  width = 10,
  height = 6
)

ggsave(
  "pairwise_relative_distance_lineplot_faceted.png",
  p_facet,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  "pairwise_relative_distance_lineplot_faceted.pdf",
  p_facet,
  width = 12,
  height = 7
)

ggsave(
  "pairwise_relative_distance_lineplot_faceted.svg",
  p_facet,
  width = 12,
  height = 7
)


