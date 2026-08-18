# ==============================================================================
# Figure 1A: Climate variables over time
# Display: Slope and R² only
# Color:
#   Temperature, Radiation = red line + red shade
#   Precipitation, Humidity = blue line + blue shade
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Package loading
# ------------------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(patchwork)
library(svglite)

# ------------------------------------------------------------------------------
# 1. Data import
# ------------------------------------------------------------------------------
df <- read_excel("figure2A_climate.xlsx", sheet = "Sheet1")

names(df) <- c(
  "Year",
  "Temp",
  "Radiation",
  "Precipitation",
  "Humidity"
)

# ------------------------------------------------------------------------------
# 2. Slope, R², and p-value calculation
# ------------------------------------------------------------------------------
fit_info <- df %>%
  pivot_longer(
    cols = c(Temp, Radiation, Precipitation, Humidity),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(variable) %>%
  summarise(
    slope = coef(lm(value ~ Year))[2],
    r2 = summary(lm(value ~ Year))$r.squared,
    p_value = summary(lm(value ~ Year))$coefficients[2, 4],
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "< 0.001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    stat_label = paste0(
      "Slope = ", sprintf("%.4f", slope), "\n",
      "R² = ", sprintf("%.2f", r2), "\n",
      "p = ", p_label
    )
  )

print(fit_info)

# ------------------------------------------------------------------------------
# 3. Panel plotting function
# ------------------------------------------------------------------------------
make_panel <- function(data, var_name, panel_title, y_lab,
                       line_col, fill_col) {
  
  plot_df <- data %>%
    select(Year, value = all_of(var_name))
  
  label_df <- fit_info %>%
    filter(variable == var_name)
  
  ggplot(plot_df, aes(x = Year, y = value)) +
    geom_point(
      size = 2.1,
      color = "black"
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = line_col,
      fill = fill_col,
      linewidth = 0.7,
      alpha = 0.25
    ) +
    annotate(
      "text",
      x = 2006,
      y = Inf,
      label = label_df$stat_label,
      hjust = 0.5,
      vjust = 1.4,
      size = 3.2,
      family = "Arial"
    ) +
    scale_x_continuous(
      breaks = c(1990, 2000, 2010, 2020),
      limits = c(1988, 2027),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = panel_title,
      x = "Year",
      y = y_lab
    ) +
    theme_classic(base_family = "Arial") +
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(
        size = 11,
        face = "bold"
      ),
      axis.text = element_text(
        size = 10,
        color = "black"
      ),
      axis.line = element_blank(),
      axis.ticks = element_line(
        linewidth = 0.5,
        color = "black"
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.6
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# ------------------------------------------------------------------------------
# 4. Color settings
# ------------------------------------------------------------------------------

red_line <- "#C44E52"
red_fill <- "#F2B6B3"

blue_line <- "#4C78A8"
blue_fill <- "#A9CBE8"

# ------------------------------------------------------------------------------
# 5. Make each panel
# ------------------------------------------------------------------------------
p_temp <- make_panel(
  data = df,
  var_name = "Temp",
  panel_title = "Mean Temperature",
  y_lab = "Temperature (°C)",
  line_col = red_line,
  fill_col = red_fill
)

p_rad <- make_panel(
  data = df,
  var_name = "Radiation",
  panel_title = "Radiation",
  y_lab = expression(Radiation~"(W/m"^2*")"),
  line_col = red_line,
  fill_col = red_fill
)

p_prec <- make_panel(
  data = df,
  var_name = "Precipitation",
  panel_title = "Precipitation",
  y_lab = "Precipitation (mm)",
  line_col = blue_line,
  fill_col = blue_fill
)

p_hum <- make_panel(
  data = df,
  var_name = "Humidity",
  panel_title = "Humidity",
  y_lab = "Humidity (%)",
  line_col = blue_line,
  fill_col = blue_fill
)

# ------------------------------------------------------------------------------
# 6. Combine panels
# ------------------------------------------------------------------------------

# 1행 4열 배열
fig_climate <- p_temp + p_rad + p_prec + p_hum +
  plot_layout(ncol = 4)

fig_climate

# 만약 기존처럼 2행 2열로 그리고 싶으면 아래 코드 사용
# fig_climate <- (p_temp + p_rad) / (p_prec + p_hum) +
#   plot_layout(ncol = 2)

# ------------------------------------------------------------------------------
# 7. Export
# ------------------------------------------------------------------------------
ggsave(
  filename = "climate_4panel_slope_R2_colored.pdf",
  plot = fig_climate,
  width = 11.5,
  height = 3.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "climate_4panel_slope_R2_colored.png",
  plot = fig_climate,
  width = 11.5,
  height = 3.2,
  units = "in",
  dpi = 600
)

ggsave(
  filename = "climate_4panel_slope_R2_colored.svg",
  plot = fig_climate,
  width = 11.5,
  height = 3.2,
  units = "in",
  device = svglite
)