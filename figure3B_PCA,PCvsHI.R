library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

#----------------------------
# 1. 데이터 불러오기
#----------------------------
file_path <- "figure2B_pca.xlsx"

df <- read_excel(file_path, sheet = 2) %>%
  mutate(.row_id = row_number())

#----------------------------
# 2. PCA 변수 설정
#----------------------------
meta_cols <- c("Country", "Continent", "HI_7mo")

pca_vars <- df %>%
  select(where(is.numeric)) %>%
  colnames()

pca_vars <- setdiff(pca_vars, c("HI_7mo", ".row_id"))

#----------------------------
# 3. Europe / South America PCA outlier 1개씩 제거
#    각 대륙 내 PCA에서 PC1-PC2 중심으로부터 가장 먼 점 제거
#----------------------------
find_one_pca_outlier <- function(data, continent_name, pca_vars) {
  
  sub_df <- data %>%
    filter(Continent == continent_name)
  
  pca_input <- sub_df %>%
    select(all_of(pca_vars)) %>%
    as.data.frame()
  
  keep <- complete.cases(pca_input)
  
  pca_input <- pca_input[keep, ]
  sub_df <- sub_df[keep, ]
  
  pca_res <- prcomp(
    pca_input,
    center = TRUE,
    scale. = TRUE
  )
  
  scores <- as.data.frame(pca_res$x[, 1:2])
  scores$.row_id <- sub_df$.row_id
  scores$Country <- sub_df$Country
  scores$Continent <- sub_df$Continent
  scores$HI_7mo <- sub_df$HI_7mo
  
  scores <- scores %>%
    mutate(
      dist_from_center = sqrt(
        (PC1 - mean(PC1, na.rm = TRUE))^2 +
          (PC2 - mean(PC2, na.rm = TRUE))^2
      )
    )
  
  scores %>%
    arrange(desc(dist_from_center)) %>%
    slice(1) %>%
    select(.row_id, Country, Continent, HI_7mo, PC1, PC2, dist_from_center)
}

outlier_rows <- bind_rows(
  find_one_pca_outlier(df, "Europe", pca_vars),
  find_one_pca_outlier(df, "South America", pca_vars)
)

print(outlier_rows)

df_plot <- df %>%
  anti_join(outlier_rows %>% select(.row_id), by = ".row_id")

#----------------------------
# 4. ellipse 함수
#----------------------------
make_axis_ellipse <- function(x, y, level = 0.95, n = 200) {
  
  cx <- mean(x, na.rm = TRUE)
  cy <- mean(y, na.rm = TRUE)
  
  sx <- sd(x, na.rm = TRUE)
  sy <- sd(y, na.rm = TRUE)
  
  r <- sqrt(qchisq(level, df = 2))
  theta <- seq(0, 2 * pi, length.out = n)
  
  data.frame(
    x = cx + r * sx * cos(theta),
    y = cy + r * sy * sin(theta)
  )
}

#----------------------------
# 5. 회귀 통계 함수
#    slope, R², p-value 표시
#----------------------------
get_lm_stats <- function(data, xvar, yvar = "HI") {
  
  fit <- lm(reformulate(xvar, yvar), data = data)
  sm <- summary(fit)
  
  slope <- coef(fit)[2]
  r2 <- sm$r.squared
  pval <- coef(sm)[2, 4]
  
  p_txt <- ifelse(
    is.na(pval),
    "NA",
    ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))
  )
  
  list(
    slope = slope,
    r2 = r2,
    p = pval,
    label = paste0(
      "Slope = ", sprintf("%.2f", slope), "\n",
      "R² = ", sprintf("%.2f", r2), "\n",
      "p = ", p_txt
    )
  )
}

#----------------------------
# 6. 공통 theme
#----------------------------
base_theme <- theme_classic() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    axis.line = element_blank(),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    axis.text = element_text(
      color = "black",
      size = 8
    ),
    axis.title = element_text(
      color = "black",
      size = 9
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 10,
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(
      size = 9,
      color = "black"
    ),
    legend.text = element_text(
      size = 8,
      color = "black"
    ),
    plot.margin = margin(4, 4, 4, 4)
  )

#----------------------------
# 7. HI 범위 설정
#----------------------------
hi_limits <- range(df_plot$HI_7mo, na.rm = TRUE)
hi_breaks <- pretty(hi_limits, n = 5)

hi_ylim <- range(df_plot$HI_7mo, na.rm = TRUE)
dy <- diff(hi_ylim)

hi_ylim <- c(
  hi_ylim[1] - dy * 0.10,
  hi_ylim[2] + dy * 0.10
)

hi_fill_scale <- scale_fill_gradientn(
  colours = c(
    "#fff5f0",
    "#fcbba1",
    "#fc9272",
    "#fb6a4a",
    "#de2d26",
    "#a50f15"
  ),
  limits = hi_limits,
  breaks = hi_breaks,
  name = "HI"
)

#----------------------------
# 8. Overall PCA 계산
#----------------------------
pca_input_all <- df_plot %>%
  select(all_of(pca_vars)) %>%
  as.data.frame()

keep_all <- complete.cases(pca_input_all)

pca_input_all <- pca_input_all[keep_all, ]
df_all <- df_plot[keep_all, ]

pca_res_all <- prcomp(
  pca_input_all,
  center = TRUE,
  scale. = TRUE
)

scores_all <- as.data.frame(pca_res_all$x[, 1:2])
scores_all$HI <- df_all$HI_7mo
scores_all$Continent <- df_all$Continent
scores_all$Country <- df_all$Country

var_expl_all <- summary(pca_res_all)$importance[2, 1:2] * 100

ellipse_all <- make_axis_ellipse(
  scores_all$PC1,
  scores_all$PC2
)

pc1_stats_all <- get_lm_stats(
  scores_all,
  "PC1",
  "HI"
)

x1_anno_all <- min(scores_all$PC1, na.rm = TRUE) +
  0.03 * diff(range(scores_all$PC1, na.rm = TRUE))

y1_anno_all <- max(scores_all$HI, na.rm = TRUE) -
  0.05 * diff(range(scores_all$HI, na.rm = TRUE))

#----------------------------
# 9. Overall PCA plot
#----------------------------
p_overall <- ggplot(
  scores_all,
  aes(PC1, PC2, fill = HI)
) +
  geom_hline(
    yintercept = 0,
    color = "grey45",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey45",
    linewidth = 0.4
  ) +
  geom_path(
    data = ellipse_all,
    aes(x, y),
    inherit.aes = FALSE,
    color = "grey30",
    linewidth = 0.4
  ) +
  geom_point(
    shape = 21,
    color = "grey25",
    stroke = 0.3,
    size = 2.8,
    alpha = 0.95
  ) +
  hi_fill_scale +
  labs(
    title = "Overall",
    x = NULL,
    y = paste0("PC2 (", round(var_expl_all[2], 1), "%)")
  ) +
  coord_equal() +
  base_theme

#----------------------------
# 10. Overall HI-PC1 plot
#----------------------------
p_overall_pc1 <- ggplot(
  scores_all,
  aes(x = PC1, y = HI)
) +
  geom_point(
    size = 2.2,
    alpha = 0.9,
    color = "grey40"
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#c62828",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = x1_anno_all,
    y = y1_anno_all,
    label = pc1_stats_all$label,
    hjust = 0,
    vjust = 1,
    size = 2.4
  ) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(var_expl_all[1], 1), "%)"),
    y = "HI"
  ) +
  coord_cartesian(ylim = hi_ylim) +
  base_theme +
  theme(
    legend.position = "none"
  )

#----------------------------
# 11. 대륙별 PCA + HI-PC1 plot 생성 함수
#----------------------------
make_continent_plots <- function(continent_name) {
  
  sub_df <- df_plot %>%
    filter(Continent == continent_name)
  
  pca_input <- sub_df %>%
    select(all_of(pca_vars)) %>%
    as.data.frame()
  
  keep <- complete.cases(pca_input)
  
  pca_input <- pca_input[keep, ]
  sub_df <- sub_df[keep, ]
  
  pca_res <- prcomp(
    pca_input,
    center = TRUE,
    scale. = TRUE
  )
  
  scores <- as.data.frame(pca_res$x[, 1:2])
  scores$HI <- sub_df$HI_7mo
  scores$Country <- sub_df$Country
  scores$Continent <- sub_df$Continent
  
  var_expl <- summary(pca_res)$importance[2, 1:2] * 100
  
  ellipse_df <- make_axis_ellipse(
    scores$PC1,
    scores$PC2
  )
  
  pc1_stats <- get_lm_stats(
    scores,
    "PC1",
    "HI"
  )
  
  x1_anno <- min(scores$PC1, na.rm = TRUE) +
    0.03 * diff(range(scores$PC1, na.rm = TRUE))
  
  y1_anno <- max(scores$HI, na.rm = TRUE) -
    0.05 * diff(range(scores$HI, na.rm = TRUE))
  
  p_pca <- ggplot(
    scores,
    aes(PC1, PC2, fill = HI)
  ) +
    geom_hline(
      yintercept = 0,
      color = "grey45",
      linewidth = 0.4
    ) +
    geom_vline(
      xintercept = 0,
      color = "grey45",
      linewidth = 0.4
    ) +
    geom_path(
      data = ellipse_df,
      aes(x, y),
      inherit.aes = FALSE,
      color = "grey30",
      linewidth = 0.4
    ) +
    geom_point(
      shape = 21,
      color = "grey25",
      stroke = 0.3,
      size = 2.8,
      alpha = 0.95
    ) +
    hi_fill_scale +
    labs(
      title = continent_name,
      x = NULL,
      y = paste0("PC2 (", round(var_expl[2], 1), "%)")
    ) +
    coord_equal() +
    base_theme
  
  p_pc1 <- ggplot(
    scores,
    aes(x = PC1, y = HI)
  ) +
    geom_point(
      size = 2.2,
      alpha = 0.9,
      color = "grey40"
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "#c62828",
      linewidth = 0.7
    ) +
    annotate(
      "text",
      x = x1_anno,
      y = y1_anno,
      label = pc1_stats$label,
      hjust = 0,
      vjust = 1,
      size = 2.4
    ) +
    labs(
      title = NULL,
      x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
      y = "HI"
    ) +
    coord_cartesian(ylim = hi_ylim) +
    base_theme +
    theme(
      legend.position = "none"
    )
  
  list(
    pca = p_pca,
    pc1 = p_pc1
  )
}

#----------------------------
# 12. 대륙별 plot 생성
#----------------------------
aus <- make_continent_plots("Australia")
eur <- make_continent_plots("Europe")
na  <- make_continent_plots("North America")
sa  <- make_continent_plots("South America")

p_aus <- aus$pca
p_eur <- eur$pca
p_na  <- na$pca
p_sa  <- sa$pca

p_aus_pc1 <- aus$pc1
p_eur_pc1 <- eur$pc1
p_na_pc1  <- na$pc1
p_sa_pc1  <- sa$pc1

#----------------------------
# 13. 2행 5열 배치
#     이전 그림 느낌을 살리기 위해 아래 plot은 보조패널 높이로 유지
#----------------------------
row1 <- p_overall | p_aus | p_eur | p_na | p_sa
row2 <- p_overall_pc1 | p_aus_pc1 | p_eur_pc1 | p_na_pc1 | p_sa_pc1

p_final_5x2_no_outlier <- (row1 / row2) +
  plot_layout(
    heights = c(1, 0.42),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

p_final_5x2_no_outlier

#----------------------------
# 14. 저장
#----------------------------
ggsave(
  filename = "Figure2B_PCA_PC1HI_5x2_EUR_SA_outlier_removed_final.png",
  plot = p_final_5x2_no_outlier,
  width = 17,
  height = 5.8,
  dpi = 300
)

ggsave(
  filename = "Figure2B_PCA_PC1HI_5x2_EUR_SA_outlier_removed_final.svg",
  plot = p_final_5x2_no_outlier,
  width = 17,
  height = 5.8
)

ggsave(
  filename = "Figure2B_PCA_PC1HI_5x2_EUR_SA_outlier_removed_final.pdf",
  plot = p_final_5x2_no_outlier,
  width = 17,
  height = 5.8
)