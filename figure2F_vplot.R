# =========================================================
# Figure 1F separate plots
# 1) volcano + histogram
# 2) negative metabolites only
# 3) positive metabolites only
# =========================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# =========================================================
# 1. read files
# =========================================================
hist_df <- read.csv("figure1F_histogram_edit.csv", check.names = FALSE)
corr_df <- read.csv("figure1F_metabolite_correlation.csv", check.names = FALSE)

# remove empty / unnamed columns
corr_df <- corr_df[, !grepl("^Unnamed:|^$|^\\s*$", colnames(corr_df))]

# =========================================================
# 2. prepare volcano / histogram data
# =========================================================
hist_df <- hist_df %>%
  mutate(
    logFDR = -log10(fdr),
    volcano_group = case_when(
      fdr < 0.05 & r_value > 0 ~ "Positive",
      fdr < 0.05 & r_value < 0 ~ "Negative",
      TRUE ~ "Nonsignificant"
    ),
    hist_group = case_when(
      r_value >= 0.9 ~ "0.9 to 1.0",
      r_value >= -0.9 & r_value < -0.8 ~ "-0.9 to -0.8",
      TRUE ~ "Others"
    )
  )

# volcano labels
label_pos_volcano <- hist_df %>%
  filter(r_value >= 0.9)

label_neg_volcano <- hist_df %>%
  filter(fdr < 0.05, r_value < 0) %>%
  slice_min(order_by = r_value, n = 5)

# =========================================================
# 3. correlation lookup
# =========================================================
corr_cols <- colnames(corr_df)
corr_cols <- corr_cols[corr_cols != "Year"]

corr_lookup <- tibble(
  corr_col = corr_cols,
  Metabolite = str_remove(corr_cols, "[0-9]+$")
)

label_pos_map <- label_pos_volcano %>%
  left_join(corr_lookup, by = "Metabolite")

label_neg_map <- label_neg_volcano %>%
  left_join(corr_lookup, by = "Metabolite")

# =========================================================
# 4. long-format metabolite data
# =========================================================
long_df <- corr_df %>%
  pivot_longer(
    cols = -Year,
    names_to = "corr_col",
    values_to = "Value"
  ) %>%
  mutate(
    Metabolite = str_remove(corr_col, "[0-9]+$")
  ) %>%
  group_by(Year, Metabolite) %>%
  summarise(
    Value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Metabolite) %>%
  mutate(
    Z = as.numeric(scale(Value))
  ) %>%
  ungroup()

pos_mets <- na.omit(unique(label_pos_map$Metabolite))
neg_mets <- na.omit(unique(label_neg_map$Metabolite))

long_pos <- long_df %>%
  filter(Metabolite %in% pos_mets)

long_neg <- long_df %>%
  filter(Metabolite %in% neg_mets)

summary_pos <- long_pos %>%
  group_by(Year) %>%
  summarise(
    mean_z = mean(Z, na.rm = TRUE),
    sd_z   = sd(Z, na.rm = TRUE),
    .groups = "drop"
  )

summary_neg <- long_neg %>%
  group_by(Year) %>%
  summarise(
    mean_z = mean(Z, na.rm = TRUE),
    sd_z   = sd(Z, na.rm = TRUE),
    .groups = "drop"
  )

# =========================================================
# 5. common axis ranges
# =========================================================
year_breaks <- seq(
  floor(min(corr_df$Year, na.rm = TRUE) / 5) * 5,
  ceiling(max(corr_df$Year, na.rm = TRUE) / 5) * 5,
  by = 5
)

year_min <- min(corr_df$Year, na.rm = TRUE)
year_max <- max(corr_df$Year, na.rm = TRUE)

all_z <- c(long_pos$Z, long_neg$Z)
z_lim <- range(all_z, na.rm = TRUE)

pad <- 0.1 * diff(z_lim)
if (is.na(pad) || pad == 0) pad <- 0.2
z_lim <- c(z_lim[1] - pad, z_lim[2] + pad)

# for outside labels
left_label_space  <- 7
right_label_space <- 7

# =========================================================
# 6. label positions
#    negative: first year point
#    positive: last year point
# =========================================================
left_label_space  <- 7
right_label_space <- 7

x_label_neg <- 1995 - 1.5
x_label_pos <- year_max + 1.5

label_neg_text <- long_neg %>%
  group_by(Metabolite) %>%
  slice_min(order_by = Year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Label = str_trunc(Metabolite, 26),
    x_text = x_label_neg
  )

label_pos_text <- long_pos %>%
  group_by(Metabolite) %>%
  slice_max(order_by = Year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Label = str_trunc(Metabolite, 26),
    x_text = x_label_pos
  )

# =========================================================
# 7. themes
# =========================================================
theme_center <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9, color = "black"),
    plot.margin = margin(5, 5, 5, 5)
  )

theme_neg <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9, color = "black"),
    plot.margin = margin(5, 110, 5, 10)
  )

theme_pos <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9, color = "black"),
    plot.margin = margin(5, 10, 5, 110)
  )

# =========================================================
# 8. volcano plot
# =========================================================
volcano_plot <- ggplot(hist_df, aes(x = r_value, y = logFDR, color = volcano_group)) +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_text_repel(
    data = label_pos_volcano,
    aes(label = Metabolite),
    color = "#d73027",
    size = 2.8,
    box.padding = 0.45,
    point.padding = 0.25,
    force = 2.5,
    force_pull = 0.5,
    min.segment.length = 0,
    segment.color = "#d73027",
    segment.size = 0.3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  geom_text_repel(
    data = label_neg_volcano,
    aes(label = Metabolite),
    color = "#2c7bb6",
    size = 2.8,
    box.padding = 0.45,
    point.padding = 0.25,
    force = 2.5,
    force_pull = 0.5,
    min.segment.length = 0,
    segment.color = "#2c7bb6",
    segment.size = 0.3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Positive" = "#d73027",
      "Negative" = "#2c7bb6",
      "Nonsignificant" = "grey75"
    )
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  labs(
    x = NULL,
    y = expression(-log[10](FDR))
  ) +
  theme_center +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# =========================================================
# 9. bottom histogram
# =========================================================
hist_plot_bottom <- ggplot(hist_df, aes(x = r_value, fill = hist_group)) +
  geom_histogram(
    aes(y = -after_stat(count)),
    binwidth = 0.1,
    boundary = -1,
    color = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_vline(
    xintercept = 0.9,
    linetype = "dashed",
    linewidth = 0.4,
    color = "#d73027"
  ) +
  geom_vline(
    xintercept = -0.9,
    linetype = "dashed",
    linewidth = 0.4,
    color = "#2c7bb6"
  ) +
  scale_fill_manual(
    values = c(
      "0.9 to 1.0" = "#d73027",
      "-0.9 to -0.8" = "#2c7bb6",
      "Others" = "grey75"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  scale_y_continuous(
    labels = abs,
    breaks = seq(0, -100, by = -10)
  ) +
  labs(
    x = "Pearson's r with year",
    y = "Count"
  ) +
  theme_center

# =========================================================
# 10. center panel save object
# =========================================================
center_panel <- (volcano_plot / hist_plot_bottom) +
  plot_layout(heights = c(3.2, 1))

# =========================================================
# 11. negative plot only
#     x-axis starts at 1995
#     labels outside left, not clipped
# =========================================================
negative_plot <- ggplot() +
  geom_line(
    data = long_neg,
    aes(x = Year, y = Z, group = Metabolite),
    color = "grey60",
    alpha = 0.35,
    linewidth = 0.4
  ) +
  geom_point(
    data = long_neg,
    aes(x = Year, y = Z, group = Metabolite),
    color = "grey60",
    alpha = 0.25,
    size = 0.6
  ) +
  geom_ribbon(
    data = summary_neg,
    aes(x = Year, ymin = mean_z - sd_z, ymax = mean_z + sd_z),
    fill = "#2c7bb6",
    alpha = 0.12
  ) +
  geom_line(
    data = summary_neg,
    aes(x = Year, y = mean_z),
    color = "#2c7bb6",
    linewidth = 1.5
  ) +
  geom_point(
    data = summary_neg,
    aes(x = Year, y = mean_z),
    color = "#2c7bb6",
    size = 1.5
  ) +
  geom_segment(
    data = label_neg_text,
    aes(x = x_text + 0.2, xend = Year, y = Z, yend = Z),
    color = "grey60",
    linewidth = 0.25
  ) +
  geom_text_repel(
    data = label_neg_text,
    aes(x = x_text, y = Z, label = Label),
    color = "grey20",
    size = 2.5,
    direction = "y",
    hjust = 1,
    nudge_x = 0,
    box.padding = 0.15,
    point.padding = 0.05,
    force = 1.2,
    force_pull = 0,
    min.segment.length = 0,
    segment.color = NA,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    limits = c(1995 - left_label_space, year_max),
    breaks = year_breaks[year_breaks >= 1995 & year_breaks <= year_max]
  ) +
  scale_y_continuous(
    limits = z_lim,
    position = "right"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Year",
    y = "Z-score",
    title = "Negative metabolites"
  ) +
  theme_neg +
  theme(
    plot.title = element_text(size = 10, face = "bold", color = "#2c7bb6")
  )

# =========================================================
# 12. positive plot only
#     same size as negative
#     labels outside right, not clipped
# =========================================================
positive_plot <- ggplot() +
  geom_line(
    data = long_pos,
    aes(x = Year, y = Z, group = Metabolite),
    color = "grey60",
    alpha = 0.35,
    linewidth = 0.4
  ) +
  geom_point(
    data = long_pos,
    aes(x = Year, y = Z, group = Metabolite),
    color = "grey60",
    alpha = 0.25,
    size = 0.6
  ) +
  geom_ribbon(
    data = summary_pos,
    aes(x = Year, ymin = mean_z - sd_z, ymax = mean_z + sd_z),
    fill = "#d73027",
    alpha = 0.12
  ) +
  geom_line(
    data = summary_pos,
    aes(x = Year, y = mean_z),
    color = "#d73027",
    linewidth = 1.5
  ) +
  geom_point(
    data = summary_pos,
    aes(x = Year, y = mean_z),
    color = "#d73027",
    size = 1.5
  ) +
  geom_segment(
    data = label_pos_text,
    aes(x = Year, xend = x_text - 0.2, y = Z, yend = Z),
    color = "grey60",
    linewidth = 0.25
  ) +
  geom_text_repel(
    data = label_pos_text,
    aes(x = x_text, y = Z, label = Label),
    color = "grey20",
    size = 2.5,
    direction = "y",
    hjust = 0,
    nudge_x = 0,
    box.padding = 0.15,
    point.padding = 0.05,
    force = 1.2,
    force_pull = 0,
    min.segment.length = 0,
    segment.color = NA,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    limits = c(1995, year_max + right_label_space),
    breaks = year_breaks[year_breaks >= 1995 & year_breaks <= year_max]
  ) +
  scale_y_continuous(
    limits = z_lim,
    position = "left"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Year",
    y = "Z-score",
    title = "Positive metabolites"
  ) +
  theme_pos +
  theme(
    plot.title = element_text(size = 10, face = "bold", color = "#d73027")
  )

# draw in RStudio
center_panel
negative_plot
positive_plot

# =========================================================
# 13. save
# =========================================================
ggsave("Figure1F_centerpanel.pdf", center_panel, width = 6.5, height = 4.8)
ggsave("Figure1F_centerpanel.png", center_panel, width = 6.5, height = 4.8, dpi = 300)
ggsave("Figure1F_centerpanel.svg", center_panel, width = 6.5, height = 4.8)

ggsave("Figure1F_negative_only.pdf", negative_plot, width = 6.5, height = 4.8)
ggsave("Figure1F_negative_only.png", negative_plot, width = 6.5, height = 4.8, dpi = 300)
ggsave("Figure1F_negative_only.svg", negative_plot, width = 6.5, height = 4.8)

ggsave("Figure1F_positive_only.pdf", positive_plot, width = 6.5, height = 4.8)
ggsave("Figure1F_positive_only.png", positive_plot, width = 6.5, height = 4.8, dpi = 300)
ggsave("Figure1F_positive_only.svg", positive_plot, width = 6.5, height = 4.8)

# =========================================================
# 14. checks
# =========================================================
cat("Matched positive rows:", nrow(long_pos), "\n")
cat("Matched negative rows:", nrow(long_neg), "\n")

cat("\nPositive metabolites used:\n")
print(sort(unique(long_pos$Metabolite)))

cat("\nNegative metabolites used:\n")
print(sort(unique(long_neg$Metabolite)))

cat("\nPositive not matched:\n")
print(label_pos_map %>% filter(is.na(corr_col)) %>% pull(Metabolite))

cat("\nNegative not matched:\n")
print(label_neg_map %>% filter(is.na(corr_col)) %>% pull(Metabolite))