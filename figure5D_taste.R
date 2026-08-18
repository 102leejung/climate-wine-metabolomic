# ==============================================================================
# FIGURE 4E: TASTE REMODELING SCORE + METABOLITE DISTRIBUTION
# ==============================================================================

library(tidyverse)
library(readxl)

# =========================================================
# 1. 파일 읽기
# =========================================================
df <- read_excel(
  "figure 4E.xlsx",
  sheet = "Metabolites-Taste",
  col_names = TRUE,
  .name_repair = "minimal"
)

# 엑셀 열 이름 강제 정리
names(df) <- c(
  "Metabolite",
  "Sig_Past",
  "Sig_Present",
  "Sig_Future",
  "Val_Past",
  "Val_Present",
  "Val_Future",
  "Taste",
  "Class_Taste"
)

# =========================================================
# 2. 유의성 별 개수 함수
# =========================================================
star_count <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x), 0, stringr::str_count(x, "\\*"))
}

# =========================================================
# 3. Past & Future 기준 대표 shift 계산
# =========================================================
df_calc <- df %>%
  filter(!is.na(Class_Taste) & Class_Taste != "") %>%
  mutate(
    Val_Past   = as.numeric(Val_Past),
    Val_Future = as.numeric(Val_Future),
    
    is_sig_past   = grepl("\\*", Sig_Past),
    is_sig_future = grepl("\\*", Sig_Future),
    
    sig_star_past   = star_count(Sig_Past),
    sig_star_future = star_count(Sig_Future),
    sig_star_max    = pmax(sig_star_past, sig_star_future, na.rm = TRUE),
    
    Sig_Level = case_when(
      sig_star_max == 0 ~ "ns",
      sig_star_max == 1 ~ "*",
      sig_star_max == 2 ~ "**",
      sig_star_max >= 3 ~ "***"
    )
  ) %>%
  rowwise() %>%
  mutate(
    same_direction = case_when(
      !is.na(Val_Past) & !is.na(Val_Future) ~ sign(Val_Past) == sign(Val_Future),
      TRUE ~ NA
    ),
    
    Plot_Value = case_when(
      # Past만 있는 경우
      !is.na(Val_Past) & is.na(Val_Future) ~ Val_Past,
      
      # Future만 있는 경우
      is.na(Val_Past) & !is.na(Val_Future) ~ Val_Future,
      
      # Past/Future 방향이 같은 경우 평균
      same_direction == TRUE ~ mean(c(Val_Past, Val_Future), na.rm = TRUE),
      
      # 방향이 반대이고 Past만 유의
      same_direction == FALSE & is_sig_past & !is_sig_future ~ Val_Past,
      
      # 방향이 반대이고 Future만 유의
      same_direction == FALSE & !is_sig_past & is_sig_future ~ Val_Future,
      
      # 방향이 반대이고 둘 다 유의 또는 둘 다 비유의이면 절댓값 큰 값
      same_direction == FALSE ~ if_else(
        abs(Val_Past) > abs(Val_Future),
        Val_Past,
        Val_Future
      ),
      
      TRUE ~ NA_real_
    ),
    
    Method_Used = case_when(
      is.na(same_direction) ~ "Single value",
      same_direction == TRUE ~ "Averaged",
      same_direction == FALSE ~ "Opposite direction",
      TRUE ~ "Other"
    ),
    
    Conflict = if_else(
      Method_Used == "Opposite direction",
      "Opposite Past/Future",
      "Consistent/single"
    ),
    
    Direction = if_else(Plot_Value >= 0, "Increase", "Decrease")
  ) %>%
  ungroup()

# =========================================================
# 4. plotting 데이터 준비
# =========================================================
taste_levels <- c("Bitterness", "Sourness", "Sweetness", "Umaminess")

df_plot <- df_calc %>%
  filter(!is.na(Plot_Value)) %>%
  mutate(
    Class_Taste = factor(Class_Taste, levels = taste_levels),
    Sig_Level = factor(Sig_Level, levels = c("ns", "*", "**", "***")),
    Conflict = factor(Conflict, levels = c("Consistent/single", "Opposite Past/Future"))
  )

# taste class별 summary
summary_df <- df_plot %>%
  group_by(Class_Taste) %>%
  summarise(
    n_metabolites = n(),
    n_increase = sum(Plot_Value > 0, na.rm = TRUE),
    n_decrease = sum(Plot_Value < 0, na.rm = TRUE),
    sum_shift = sum(Plot_Value, na.rm = TRUE),
    mean_shift = mean(Plot_Value, na.rm = TRUE),
    median_shift = median(Plot_Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_pos = as.numeric(factor(Class_Taste, levels = rev(taste_levels))),
    label = paste0(
      "Sum = ", sprintf("%+.1f", sum_shift), "%\n",
      "n=", n_metabolites,
      "  (+", n_increase, " / -", n_decrease, ")"
    )
  )

df_plot <- df_plot %>%
  mutate(
    y_pos = as.numeric(factor(Class_Taste, levels = rev(taste_levels)))
  )

# =========================================================
# 5. 색상 설정
# =========================================================
taste_cols <- c(
  "Bitterness" = "#C9A227",
  "Sourness"   = "#4FA78A",
  "Sweetness"  = "#D98C6A",
  "Umaminess"  = "#9E9E9E"
)

direction_cols <- c(
  "Increase" = "#D65F5F",
  "Decrease" = "#4C78A8"
)

# =========================================================
# 6. 배경 band용 데이터
# =========================================================
band_df <- summary_df %>%
  mutate(
    ymin = y_pos - 0.45,
    ymax = y_pos + 0.45
  )

# x축 범위 및 오른쪽 label 위치
x_min <- min(df_plot$Plot_Value, summary_df$mean_shift, na.rm = TRUE)
x_max <- max(df_plot$Plot_Value, summary_df$mean_shift, na.rm = TRUE)

x_label <- x_max + 18

# =========================================================
# 7. Taste remodeling score plot
# =========================================================
plot_taste_score <- ggplot() +
  
  # taste class background bands
  geom_rect(
    data = band_df,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax,
      fill = Class_Taste
    ),
    alpha = 0.10,
    color = NA
  ) +
  
  # zero 기준선
  geom_vline(
    xintercept = 0,
    color = "grey25",
    linewidth = 0.7
  ) +
  
  # reference lines
  geom_vline(
    xintercept = c(-25, 25, 50, 75),
    color = "grey88",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  
  # class-level mean shift segment
  geom_segment(
    data = summary_df,
    aes(
      x = 0,
      xend = mean_shift,
      y = y_pos,
      yend = y_pos,
      color = Class_Taste
    ),
    linewidth = 4.5,
    alpha = 0.35,
    lineend = "round"
  ) +
  
  # metabolite-level points
  geom_point(
    data = df_plot,
    aes(
      x = Plot_Value,
      y = y_pos,
      shape = Conflict,
      size = Sig_Level,
      fill = Direction
    ),
    position = position_jitter(height = 0.15, width = 0, seed = 123),
    color = "black",
    stroke = 0.45,
    alpha = 0.90
  ) +
  
  # class-level mean shift diamond
  geom_point(
    data = summary_df,
    aes(
      x = mean_shift,
      y = y_pos,
      color = Class_Taste
    ),
    shape = 18,
    size = 5.5
  ) +
  
  # right-side summary label
  geom_text(
    data = summary_df,
    aes(
      x = x_label,
      y = y_pos,
      label = label,
      color = Class_Taste
    ),
    hjust = 0,
    size = 3.7,
    fontface = "bold",
    lineheight = 0.95
  ) +
  
  # y축 taste class label
  scale_y_continuous(
    breaks = summary_df$y_pos,
    labels = summary_df$Class_Taste,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  
  scale_x_continuous(
    limits = c(x_min - 8, x_label + 35),
    breaks = c(-25, 0, 25, 50, 75),
    expand = c(0, 0)
  ) +
  
  scale_fill_manual(
    values = c(taste_cols, direction_cols),
    guide = "none"
  ) +
  
  scale_color_manual(
    values = c(taste_cols, direction_cols),
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = c(
      "Consistent/single" = 21,
      "Opposite Past/Future" = 24
    ),
    name = "Past/Future direction"
  ) +
  
  scale_size_manual(
    values = c(
      "ns" = 1.8,
      "*" = 2.4,
      "**" = 3.0,
      "***" = 3.7
    ),
    name = "Significance"
  ) +
  
  labs(
    x = "Representative shift rate (%)",
    y = NULL,
    title = "Taste profile remodeling by climate-associated metabolites",
    subtitle = "Small points indicate individual metabolites; diamonds and thick segments indicate class-level mean shifts. Right labels show cumulative shift and increase/decrease counts."
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "grey35", hjust = 0),
    
    axis.title.x = element_text(size = 11.5, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_rect(color = "grey70", linewidth = 0.7),
    
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    
    plot.margin = margin(10, 90, 10, 10)
  )

print("SUCCESS! Figure 4E taste remodeling score plot has been generated.")

plot_taste_score
# =========================================================
# 8. 저장
# =========================================================
ggsave(
  "Figure_4E_taste_remodeling_score.png",
  plot = plot_taste_score,
  width = 9.5,
  height = 5.8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure_4E_taste_remodeling_score.pdf",
  plot = plot_taste_score,
  width = 9.5,
  height = 5.8,
  device = cairo_pdf,
  bg = "white"
)

print("SUCCESS! Figure 4E taste remodeling score plot has been generated.")