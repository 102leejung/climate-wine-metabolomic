library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

file_path <- "figure2A_tem,HI.xlsx"

# --------------------------------
# 1. 색상
# --------------------------------
continent_colors <- c(
  "North America" = "#2E6F95",
  "South America" = "#D1495B",
  "Europe"        = "#2A9D8F",
  "Australia"     = "#E9C46A",
  "Africa"  = "#8F77B5"
)

continent_order <- c(
  "North America",
  "South America",
  "Europe",
  "Australia",
  "Africa"
)

# --------------------------------
# 2. Temperature
# --------------------------------
temp_raw <- read_excel(file_path, sheet = "2.A Temperature Graph")

temp_long <- temp_raw %>%
  pivot_longer(
    cols = `4/10`:`10/4`,
    names_to = "month_label",
    values_to = "temperature"
  ) %>%
  mutate(
    Continent = factor(Continent, levels = continent_order),
    month_num = match(month_label, c("4/10", "5/11", "6/12", "7/1", "8/2", "9/3", "10/4"))
  )

temp_summary <- temp_long %>%
  filter(!is.na(Continent)) %>%
  group_by(Continent, month_num) %>%
  summarise(
    mean = mean(temperature, na.rm = TRUE),
    sd   = sd(temperature, na.rm = TRUE),
    ymin = mean - sd,
    ymax = mean + sd,
    .groups = "drop"
  )

# --------------------------------
# 3. HI
# --------------------------------
hi_raw <- read_excel(file_path, sheet = "2.A HI Graph")

hi_long <- hi_raw %>%
  pivot_longer(
    cols = M4_10:M10_4,
    names_to = "month_label",
    values_to = "hi"
  ) %>%
  mutate(
    Continent = factor(Continent, levels = continent_order),
    month_num = match(month_label, c("M4_10", "M5_11", "M6_12", "M7_1", "M8_2", "M9_3", "M10_4"))
  )

hi_summary <- hi_long %>%
  filter(!is.na(Continent)) %>%
  group_by(Continent, month_num) %>%
  summarise(
    mean = mean(hi, na.rm = TRUE),
    sd   = sd(hi, na.rm = TRUE),
    ymin = mean - sd,
    ymax = mean + sd,
    .groups = "drop"
  )

# --------------------------------
# 4. theme
# --------------------------------
fig_theme <- theme_classic(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

# --------------------------------
# 5. Temperature plot
# --------------------------------
p_temp <- ggplot(temp_summary, aes(x = month_num, y = mean, color = Continent, fill = Continent)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = continent_colors, breaks = continent_order) +
  scale_fill_manual(values = continent_colors, breaks = continent_order) +
  scale_x_continuous(breaks = 1:7, labels = 1:7) +
  labs(
    title = "Monthly Mean Temperature",
    x = "Month of Growing Season",
    y = "Temperature (°C)",
    color = "Continent",
    fill = "Continent"
  ) +
  fig_theme

# --------------------------------
# 6. HI plot
# --------------------------------
p_hi <- ggplot(hi_summary, aes(x = month_num, y = mean, color = Continent, fill = Continent)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = continent_colors, breaks = continent_order) +
  scale_fill_manual(values = continent_colors, breaks = continent_order) +
  scale_x_continuous(breaks = 1:7, labels = 1:7) +
  labs(
    title = "Monthly Huglin Index Accumulation",
    x = "Month of Growing Season",
    y = "Huglin index units",
    color = "Continent",
    fill = "Continent"
  ) +
  fig_theme

# --------------------------------
# 7. combine
# --------------------------------
final_plot <- p_temp + p_hi + plot_layout(ncol = 2)
final_plot

ggsave("figure2_lineplots_fixed2.png", final_plot, width = 8, height = 4, dpi = 600)
ggsave("figure2_lineplots_fixed2.pdf", final_plot, width = 8, height = 4)
ggsave("figure2_lineplots_fixed2.svg", final_plot, width = 8, height = 4)