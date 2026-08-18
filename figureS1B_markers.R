# =========================================================
# Figure 4B: Metabolite sensitivity beeswarm plot
# Input file: figure4B.xlsx
# Sheets: Fig 4B Past / Fig 4B Present / Fig 4B Future
# =========================================================

# 0. Load required libraries
install.packages("ggbeeswarm")

library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(scales)

# ---------------------------------------------------------
# 1. FILE PATH
# ---------------------------------------------------------
# figure4B.xlsx 파일이 현재 작업 폴더에 있으면 그대로 사용
file_path <- "figure4B.xlsx"

# 저장 경로: 현재 작업 폴더
out_dir <- getwd()

# ---------------------------------------------------------
# 2. LOAD AND CLEAN DATA
# ---------------------------------------------------------
clean_data <- function(sheet_name, time_label) {
  
  df <- read_excel(file_path, sheet = sheet_name)
  
  col_name <- intersect(
    names(df),
    c("Short_Name", "Metabolite", "Metabolite_Raw_Name")
  )[1]
  
  df %>%
    rename(
      Metabolite = all_of(col_name),
      Percent_Change = Percent_Change_per_C
    ) %>%
    mutate(
      Percent_Change = as.numeric(Percent_Change),
      Significance = as.character(Significance),
      Timeframe = time_label,
      Category = case_when(
        Significance == "ns" ~ "Not Significant",
        Percent_Change > 0 ~ "Significant Increase",
        Percent_Change <= 0 ~ "Significant Decrease",
        TRUE ~ "Not Significant"
      )
    ) %>%
    select(
      Metabolite,
      Percent_Change,
      Significance,
      Category,
      Timeframe
    ) %>%
    drop_na(Metabolite, Percent_Change)
}

past_df    <- clean_data("Fig 4B Past", "Past")
present_df <- clean_data("Fig 4B Present", "Present")
future_df  <- clean_data("Fig 4B Future", "Future")

# ---------------------------------------------------------
# 3. COMBINE DATA
# ---------------------------------------------------------
plot_data <- bind_rows(past_df, present_df, future_df) %>%
  mutate(
    Timeframe = factor(Timeframe, levels = c("Past", "Present", "Future")),
    Category = factor(
      Category,
      levels = c(
        "Not Significant",
        "Significant Decrease",
        "Significant Increase"
      )
    )
  ) %>%
  # Not significant 먼저 그리고, significant 점들이 위에 올라오도록 정렬
  arrange(Category)

# ---------------------------------------------------------
# 4. MAKE PLOT
# ---------------------------------------------------------
final_plot <- ggplot(
  plot_data,
  aes(
    x = "",
    y = Percent_Change,
    fill = Category,
    color = Category
  )
) +
  
  # y = 0 기준선
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.8
  ) +
  
  # Beeswarm dot
  geom_quasirandom(
    shape = 21,
    size = 3.5,
    stroke = 0.5,
    width = 0.30,
    alpha = 0.85
  ) +
  
  facet_wrap(~ Timeframe, nrow = 1) +
  
  # Dot 내부 색
  scale_fill_manual(
    values = c(
      "Not Significant" = "gray80",
      "Significant Decrease" = "#4575b4",
      "Significant Increase" = "#d73027"
    )
  ) +
  
  # Dot 테두리 색
  scale_color_manual(
    values = c(
      "Not Significant" = "gray80",
      "Significant Decrease" = "black",
      "Significant Increase" = "black"
    )
  ) +
  
  # pseudo-log y축
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 0.01, base = 10),
    breaks = c(-100, -20, -2, -0.2, 0, 0.2, 2, 20, 100, 200),
    labels = c("-100", "-20", "-2", "-0.2", "0", "0.2", "2", "20", "100", "200")
  ) +
  
  labs(
    title = "Metabolite Sensitivities Across Timelines",
    subtitle = "Y-axis scaled via pseudo-log transformation (sigma = 0.01)",
    x = NULL,
    y = "Percent Change per 1°C (%)",
    fill = "Metabolite Status",
    color = "Metabolite Status"
  ) +
  
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 16,
      color = "gray20"
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      color = "gray20"
    ),
    
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5,
      color = "gray30"
    ),
    
    axis.title.y = element_text(
      size = 13,
      face = "bold",
      color = "gray20"
    ),
    
    axis.text.y = element_text(
      size = 11,
      color = "gray20"
    ),
    
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

final_plot

# ---------------------------------------------------------
# 5. EXPORT PDF, PNG, SVG
# ---------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "Figure4B_Metabolite_Beeswarm.pdf"),
  plot = final_plot,
  width = 10,
  height = 7,
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "Figure4B_Metabolite_Beeswarm.png"),
  plot = final_plot,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "Figure4B_Metabolite_Beeswarm.svg"),
  plot = final_plot,
  width = 10,
  height = 7
)

# 화면에 출력
print(final_plot)

cat("Done! Files saved in:", out_dir, "\n")