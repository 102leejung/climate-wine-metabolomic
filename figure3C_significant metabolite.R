# =========================================================
# Representative metabolites
# Main figure: representative negative + positive metabolites
# continent lines = continent colors
# median line = direction color
# range = ribbon
# =========================================================

library(tidyverse)
library(patchwork)

# ---------------------------------------------------------
# 0. Colors
# ---------------------------------------------------------

col_positive <- "red"
col_negative <- "#0072B2"

# Continent color palette
continent_cols <- c(
  "North America" = "#2E6F95",
  "South America" = "#D1495B",
  "Europe"        = "#2A9D8F",
  "Australia"     = "#E9C46A",
  "Oceania"       = "#E9C46A",
  "Africa"        = "#8F77B5",
  "Asia"          = "grey70"
)

# ---------------------------------------------------------
# 1. Selected representative metabolites
# ---------------------------------------------------------

top_neg_mets <- c(
  "(1R)-1-(2-Furyl)ethanol",
  "Kaempferol 3-gluco-xyloside",
  "Shikimic acid",
  "Trehalose"
)

top_pos_mets <- c(
  "Xylitol",
  "Sucrose",
  "isothreonic acid",
  "Fusarenon"
)

# check column names
setdiff(top_neg_mets, colnames(df))
setdiff(top_pos_mets, colnames(df))

# ---------------------------------------------------------
# 2. Plot function
# ---------------------------------------------------------

make_sig_plot <- function(metabolite_list,
                          direction_label,
                          median_color,
                          ncol_value = 2,
                          title_text = NULL,
                          show_legend = FALSE) {
  
  trend_data <- df %>%
    select(Continent, HI_7mo, all_of(metabolite_list)) %>%
    pivot_longer(
      cols = all_of(metabolite_list),
      names_to = "Metabolite",
      values_to = "Value"
    ) %>%
    mutate(Value = as.numeric(Value)) %>%
    filter(
      is.finite(HI_7mo),
      is.finite(Value),
      !is.na(Continent)
    ) %>%
    group_by(Metabolite) %>%
    mutate(Z = as.numeric(scale(Value))) %>%
    ungroup()
  
  x_grid <- seq(
    min(trend_data$HI_7mo, na.rm = TRUE),
    max(trend_data$HI_7mo, na.rm = TRUE),
    length.out = 100
  )
  
  # continent regression lines
  pred_continent <- trend_data %>%
    group_by(Metabolite, Continent) %>%
    group_modify(~{
      
      if (nrow(.x) < 3) return(tibble())
      
      fit <- lm(Z ~ HI_7mo, data = .x)
      
      tibble(
        HI_7mo = x_grid,
        Pred = predict(
          fit,
          newdata = tibble(HI_7mo = x_grid)
        )
      )
    }) %>%
    ungroup()
  
  # median + IQR range
  summary_pred <- pred_continent %>%
    group_by(Metabolite, HI_7mo) %>%
    summarise(
      median_pred = median(Pred, na.rm = TRUE),
      ymin = quantile(Pred, 0.25, na.rm = TRUE),
      ymax = quantile(Pred, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (is.null(title_text)) {
    title_text <- paste0(direction_label, " representative metabolites")
  }
  
  p <- ggplot() +
    
    # raw sample points
    geom_point(
      data = trend_data,
      aes(x = HI_7mo, y = Z),
      color = "grey80",
      size = 0.45,
      alpha = 0.18
    ) +
    
    # IQR range across continent regression lines
    geom_ribbon(
      data = summary_pred,
      aes(
        x = HI_7mo,
        ymin = ymin,
        ymax = ymax
      ),
      fill = "grey70",
      alpha = 0.25
    ) +
    
    # continent regression lines: colored by continent
    geom_line(
      data = pred_continent,
      aes(
        x = HI_7mo,
        y = Pred,
        group = Continent,
        color = Continent
      ),
      linewidth = 0.4,
      alpha = 0.6
    ) +
    
    # overall median trend line
    geom_line(
      data = summary_pred,
      aes(
        x = HI_7mo,
        y = median_pred
      ),
      color = median_color,
      linewidth = 1.5
    ) +
    
    scale_color_manual(
      values = continent_cols,
      drop = FALSE
    ) +
    
    facet_wrap(
      ~ Metabolite,
      scales = "fixed",
      ncol = ncol_value
    ) +
    
    labs(
      x = "HI_7mo",
      y = "Z-score",
      title = title_text,
      color = "Continent"
    ) +
    
    theme_classic(base_size = 11) +
    
    theme(
      legend.position = ifelse(show_legend, "bottom", "none"),
      
      legend.title = element_text(
        size = 9,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 8
      ),
      
      strip.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.5
      ),
      
      strip.text = element_text(
        size = 8.5,
        face = "bold",
        color = "black"
      ),
      
      axis.text = element_text(
        size = 8,
        color = "black"
      ),
      
      axis.title = element_text(
        size = 10.5,
        face = "bold"
      ),
      
      plot.title = element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      )
    )
  
  return(p)
}

# =========================================================
# 3. Negative metabolites
# =========================================================

# stronger decrease
neg_strong_mets <- c(
  "(1R)-1-(2-Furyl)ethanol",
  "Shikimic acid"
)

# milder decrease
neg_mild_mets <- c(
  "Kaempferol 3-gluco-xyloside",
  "Trehalose"
)

# upper row
p_neg_strong <- make_sig_plot(
  metabolite_list = neg_strong_mets,
  direction_label = "Negative",
  median_color = col_negative,
  ncol_value = 2,
  title_text = "",
  show_legend = FALSE
) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank()
  )

# lower row
p_neg_mild <- make_sig_plot(
  metabolite_list = neg_mild_mets,
  direction_label = "Negative",
  median_color = col_negative,
  ncol_value = 2,
  title_text = "",
  show_legend = FALSE
) +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank()
  )

# combine negative
p_negative_top4 <- p_neg_strong / p_neg_mild +
  plot_annotation(
    title = "Representative negative metabolites"
  ) &
  theme(
    plot.title = element_text(
      size = 13,
      face = "bold",
      hjust = 0.5
    )
  )

# =========================================================
# 4. Positive metabolites
# =========================================================

p_positive_top4 <- make_sig_plot(
  metabolite_list = top_pos_mets,
  direction_label = "Positive",
  median_color = col_positive,
  ncol_value = 2,
  title_text = "Representative positive metabolites",
  show_legend = FALSE
) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme(
    axis.title.y = element_blank()
  )

# show plots
p_negative_top4
p_positive_top4

# =========================================================
# 5. Save negative figure
# =========================================================

ggsave(
  "Figure_negative_representative_2x2.png",
  p_negative_top4,
  width = 6,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure_negative_representative_2x2.pdf",
  p_negative_top4,
  width = 6,
  height = 5.5,
  bg = "white"
)

ggsave(
  "Figure_negative_representative_2x2.svg",
  p_negative_top4,
  width = 6,
  height = 5.5,
  bg = "white"
)

# =========================================================
# 6. Save positive figure
# =========================================================

ggsave(
  "Figure_positive_representative_2x2.png",
  p_positive_top4,
  width = 6,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure_positive_representative_2x2.pdf",
  p_positive_top4,
  width = 6,
  height = 5.5,
  bg = "white"
)

ggsave(
  "Figure_positive_representative_2x2.svg",
  p_positive_top4,
  width = 6,
  height = 5.5,
  bg = "white"
)

# =========================================================
# 7. Combine: Negative | Volcano | Positive
# =========================================================

p_combined_main <- p_negative_top4 | Figure2C_main | p_positive_top4

p_combined_main <- p_combined_main +
  plot_layout(widths = c(1.05, 1.2, 1.05))

p_combined_main

# =========================================================
# 8. Save combined main figure
# =========================================================

ggsave(
  "Figure_main_negative_volcano_positive.png",
  p_combined_main,
  width = 16,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure_main_negative_volcano_positive.pdf",
  p_combined_main,
  width = 16,
  height = 5.5,
  bg = "white"
)

ggsave(
  "Figure_main_negative_volcano_positive.svg",
  p_combined_main,
  width = 16,
  height = 5.5,
  bg = "white"
)