############################################
# Effect size dot plot
# - Dot: individual metabolites
# - Color: positive / negative slope
# - Transparency + size: significance
# - Black point + line: median ± IQR
# - Labels: top positive / negative metabolites
############################################

library(tidyverse)
library(readxl)
library(scales)
library(ggrepel)

############################################
# 1. Component order
############################################

custom_order <- c(
  "Glucose", "Fructose", "Sucrose", 
  "Tartaric Acid", "Malic Acid", "Citric Acid", 
  "Tannic Acid", 
  "Proline", "Arginine", "Glutamic Acid", "Alanine", 
  "Serine", "Threonine", "Valine", "Histidine", "Isoleucine"
)

############################################
# 2. Load data
############################################

wine_stats <- read_excel(
  "figure3D_bar, heatmap.xlsx",
  sheet = 1
)

############################################
# 3. Prepare data
############################################

dot_data <- wine_stats %>%
  filter(Component %in% custom_order) %>%
  filter(!is.na(slope), !is.na(p_value)) %>%
  mutate(
    Component = factor(Component, levels = custom_order),
    Trend = case_when(
      slope > 0 ~ "Positive",
      slope < 0 ~ "Negative",
      TRUE ~ NA_character_
    ),
    Significant = ifelse(
      p_value <= 0.05,
      "p ≤ 0.05",
      "NS"
    )
  ) %>%
  filter(!is.na(Trend))

############################################
# 4. Select top metabolites to label
############################################

top_pos <- dot_data %>%
  filter(slope > 0) %>%
  arrange(desc(slope)) %>%
  slice_head(n = 6)

top_neg <- dot_data %>%
  filter(slope < 0) %>%
  arrange(slope) %>%
  slice_head(n = 6)

top_labels <- bind_rows(top_pos, top_neg)

View(top_labels)

############################################
# 5. Make dot plot
############################################

dot_plot <- ggplot(
  dot_data,
  aes(
    x = Component,
    y = slope
  )
) +
  
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_jitter(
    aes(
      fill = Trend,
      alpha = Significant,
      size = Significant
    ),
    shape = 21,
    color = "black",
    stroke = 0.15,
    width = 0.25,
    height = 0
  ) +
  
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = median(x, na.rm = TRUE),
        ymin = quantile(x, 0.25, na.rm = TRUE),
        ymax = quantile(x, 0.75, na.rm = TRUE)
      )
    },
    geom = "pointrange",
    color = "black",
    linewidth = 0.5,
    size = 0.35,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    data = top_labels,
    aes(
      x = Component,
      y = slope,
      label = Metabolite
    ),
    inherit.aes = FALSE,
    size = 2.2,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2,
    segment.size = 0.25,
    min.segment.length = 0
  ) +
  
  scale_fill_manual(
    values = c(
      "Positive" = "#ff4122",
      "Negative" = "steelblue"
    )
  ) +
  
  scale_alpha_manual(
    values = c(
      "p ≤ 0.05" = 0.9,
      "NS" = 0.15
    )
  ) +
  
  scale_size_manual(
    values = c(
      "p ≤ 0.05" = 1.5,
      "NS" = 0.8
    )
  ) +
  
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 0.002)
  ) +
  
  labs(
    title = "Metabolite responses by wine component",
    x = NULL,
    y = "Effect size (slope, pseudo-log scale)",
    fill = NULL,
    alpha = NULL,
    size = NULL
  ) +
  
  theme_classic(base_size = 10) +
  
  theme(
    legend.position = "top",
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8,
      face = "bold"
    ),
    
    axis.text.y = element_text(size = 8),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

############################################
# 6. Display
############################################

print(dot_plot)

############################################
# 7. Save
############################################

ggsave(
  "Figure_effect_size_dot_plot_pseudolog_median_IQR_labels.png",
  plot = dot_plot,
  width = 7.5,
  height = 4.5,
  dpi = 300
)

ggsave(
  "Figure_effect_size_dot_plot_pseudolog_median_IQR_labels.pdf",
  plot = dot_plot,
  width = 7.5,
  height = 4.5
)

ggsave(
  "Figure_effect_size_dot_plot_pseudolog_median_IQR_labels.svg",
  plot = dot_plot,
  width = 7.5,
  height = 4.5
)