# =========================================================
# Figure 1B: Supporting wine chemistry block
# Input file: figure1B_all.csv
#
# Layout:
#   Left  : pH / Anthocyanin stacked mini plots
#           mean ± SD by vintage year
#   Right : compact forest summary for 6 chemistry traits
#
# Output:
#   figure1B_supporting_chemistry_stack_forest_cowplot.pdf
#   figure1B_supporting_chemistry_stack_forest_cowplot.svg
# =========================================================

library(tidyverse)
library(cowplot)
library(svglite)

# ---------------------------------------------------------
# 0. Turn off showtext if it is active
#    This helps keep text editable in Illustrator/PDF
# ---------------------------------------------------------
if ("showtext" %in% loadedNamespaces()) {
  showtext::showtext_auto(FALSE)
}

graphics.off()

# ---------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------
df_raw <- read.csv("figure2C_all.csv", check.names = FALSE)
colnames(df_raw)[1] <- "Year"

year_order <- c(1995, 1996, 2000, 2003, 2006, 2009, 2010, 2012, 2015, 2016, 2020, 2022)

df_long <- df_raw %>%
  pivot_longer(
    cols = -Year,
    names_to = c("Trait", "Replicate"),
    names_pattern = "(.+)_([123])$",
    values_to = "Value"
  ) %>%
  mutate(
    Year = as.numeric(Year),
    Value = as.numeric(Value),
    Year_f = factor(Year, levels = year_order)
  )

# ---------------------------------------------------------
# 2. Year color palette
# ---------------------------------------------------------
year_colors <- c(
  "1995" = "#F2D6DB",
  "1996" = "#E9C4CC",
  "2000" = "#DFA4B0",
  "2003" = "#D07A90",
  "2006" = "#C15C6F",
  "2009" = "#B23F62",
  "2010" = "#A12F57",
  "2012" = "#8F254C",
  "2015" = "#6F183C",
  "2016" = "#5D1233",
  "2020" = "#4A0B29",
  "2022" = "#33051C"
)

# ---------------------------------------------------------
# 3. Trait label table
# ---------------------------------------------------------
trait_label_df <- tibble(
  Trait = c(
    "pH",
    "TA",
    "Total Phenol (mg GAE/L)",
    "Flavonoid (mg QE/L)",
    "Tannin (mg TAE/L)",
    "Anthocyanin (mg C3GE/L)"
  ),
  Trait_label = c(
    "pH",
    "TA",
    "Total phenol",
    "Flavonoid",
    "Tannin",
    "Anthocyanin"
  ),
  Trait_label_short = c(
    "pH",
    "TA",
    "Phenol",
    "Flav.",
    "Tannin",
    "Antho."
  )
)

# ---------------------------------------------------------
# 4. Mini mean ± SD plot function
# ---------------------------------------------------------
make_mini_mean_sd <- function(trait_name, y_lab, title_lab) {
  
  df_sub <- df_long %>%
    filter(Trait == trait_name)
  
  yearly_summary <- df_sub %>%
    group_by(Year, Year_f) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot() +
    
    # Light linear trend based on all replicate values
    geom_smooth(
      data = df_sub,
      aes(x = Year, y = Value),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "gray65",
      linewidth = 0.35
    ) +
    
    # Mean ± SD
    geom_errorbar(
      data = yearly_summary,
      aes(
        x = Year,
        ymin = Mean - SD,
        ymax = Mean + SD,
        color = Year_f
      ),
      width = 0,
      linewidth = 0.35,
      alpha = 0.75
    ) +
    
    geom_point(
      data = yearly_summary,
      aes(x = Year, y = Mean, fill = Year_f),
      shape = 21,
      color = "black",
      stroke = 0.18,
      size = 1.35,
      alpha = 0.95
    ) +
    
    scale_color_manual(values = year_colors, drop = FALSE) +
    scale_fill_manual(values = year_colors, drop = FALSE) +
    
    scale_x_continuous(
      breaks = c(1995, 2005, 2015, 2022),
      limits = c(1994, 2023),
      expand = c(0.01, 0.01)
    ) +
    
    labs(
      title = title_lab,
      x = NULL,
      y = y_lab
    ) +
    
    theme_classic(base_size = 7.5, base_family = "sans") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8.2, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 6.2, color = "black"),
      axis.title.y = element_text(size = 6.8, color = "black"),
      axis.title.x = element_blank(),
      axis.line = element_line(linewidth = 0.32),
      axis.ticks = element_line(linewidth = 0.32),
      plot.margin = margin(1, 1, 1, 1)
    )
}

# ---------------------------------------------------------
# 5. Create pH and Anthocyanin panels
# ---------------------------------------------------------
p_pH <- make_mini_mean_sd(
  trait_name = "pH",
  y_lab = "pH",
  title_lab = "pH"
)

p_antho <- make_mini_mean_sd(
  trait_name = "Anthocyanin (mg C3GE/L)",
  y_lab = "Anthocyanin",
  title_lab = "Anthocyanin"
)

# ---------------------------------------------------------
# 6. Forest summary data
# ---------------------------------------------------------
forest_df <- df_long %>%
  filter(Trait %in% trait_label_df$Trait) %>%
  group_by(Trait) %>%
  summarise(
    r = cor(Year, Value, use = "complete.obs", method = "pearson"),
    p = cor.test(Year, Value, method = "pearson")$p.value,
    n = sum(!is.na(Value)),
    .groups = "drop"
  ) %>%
  mutate(
    z = atanh(r),
    se_z = 1 / sqrt(n - 3),
    ci_low = tanh(z - 1.96 * se_z),
    ci_high = tanh(z + 1.96 * se_z)
  ) %>%
  left_join(trait_label_df, by = "Trait") %>%
  mutate(
    Direction = case_when(
      p < 0.05 & r < 0 ~ "Decrease",
      p < 0.05 & r > 0 ~ "Increase",
      TRUE ~ "NS"
    ),
    Sig_label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

trait_order_short <- c("pH", "TA", "Phenol", "Flav.", "Tannin", "Antho.")

forest_df <- forest_df %>%
  mutate(
    Trait_label_short = factor(Trait_label_short, levels = rev(trait_order_short)),
    sig_x = case_when(
      r > 0 ~ pmin(r + 0.08, 0.84),
      r < 0 ~ pmax(r - 0.08, -0.84),
      TRUE ~ r
    )
  )

# ---------------------------------------------------------
# 7. Compact forest plot
# ---------------------------------------------------------
p_forest <- ggplot(forest_df, aes(y = Trait_label_short)) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "gray35"
  ) +
  
  geom_segment(
    aes(
      x = ci_low,
      xend = ci_high,
      y = Trait_label_short,
      yend = Trait_label_short,
      color = Direction
    ),
    linewidth = 0.55
  ) +
  
  geom_point(
    aes(x = r, color = Direction),
    size = 2.0
  ) +
  
  # 별표가 너무 답답하면 이 geom_text 블록만 주석 처리해도 됨
  geom_text(
    data = forest_df %>% filter(Sig_label != ""),
    aes(x = sig_x, label = Sig_label),
    size = 2.4,
    family = "sans",
    color = "black"
  ) +
  
  scale_color_manual(
    values = c(
      "Decrease" = "#2C6DA4",
      "Increase" = "#B23F4D",
      "NS" = "gray55"
    )
  ) +
  
  scale_x_continuous(
    limits = c(-0.9, 0.9),
    breaks = c(-0.8, 0, 0.8),
    expand = c(0.01, 0.01)
  ) +
  
  labs(
    title = "Chemistry summary",
    x = "r with year",
    y = NULL
  ) +
  
  theme_classic(base_size = 7.5, base_family = "sans") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 8.2, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 6.8, color = "black"),
    axis.text.x = element_text(size = 6.3, color = "black"),
    axis.title.x = element_text(size = 6.8, color = "black"),
    axis.line = element_line(linewidth = 0.32),
    axis.ticks = element_line(linewidth = 0.32),
    plot.margin = margin(1, 1.5, 1, 1)
  )

# ---------------------------------------------------------
# 8. Combine using cowplot
#    This avoids patchwork PDF saving issues
# ---------------------------------------------------------
mini_stack <- cowplot::plot_grid(
  p_pH,
  p_antho,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

final_plot <- cowplot::plot_grid(
  mini_stack,
  p_forest,
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(0.95, 1.25)
)

print(final_plot)

# ---------------------------------------------------------
# 9. Save files
# ---------------------------------------------------------
out_dir <- getwd()

pdf_file <- file.path(out_dir, "figure1B_supporting_chemistry_stack_forest_cowplot.pdf")
svg_file <- file.path(out_dir, "figure1B_supporting_chemistry_stack_forest_cowplot.svg")

# 기존 파일이 열려 있으면 삭제가 실패할 수 있음
# Illustrator / Acrobat / Chrome에서 해당 파일을 닫은 뒤 실행
if (file.exists(pdf_file)) file.remove(pdf_file)
if (file.exists(svg_file)) file.remove(svg_file)

# SVG 저장도 동일 크기
ggsave(
  filename = svg_file,
  plot = final_plot,
  width = 4.4,
  height = 2.4,
  device = svglite::svglite,
  bg = "white"
)

# PDF 저장
grDevices::pdf(
  file = pdf_file,
  width = 4.4,
  height = 2.4,
  family = "sans",
  useDingbats = FALSE,
  onefile = FALSE
)

print(final_plot)
dev.off()

# 저장 확인
cat("SVG saved to:\n", normalizePath(svg_file), "\n\n")
cat("PDF saved to:\n", normalizePath(pdf_file), "\n")