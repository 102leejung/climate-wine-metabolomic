library(tidyverse)
library(ggplot2)

# 1. 데이터 불러오기
df <- read.csv("figure2B_C13C12.csv")

# 2. long format 변환
df_long <- df %>%
  pivot_longer(
    cols = -Year,
    names_to = "Replicate",
    values_to = "Ratio"
  ) %>%
  drop_na()

# 3. 연도별 색상 지정
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

# 4. 연도별 평균값 계산
df_mean <- df_long %>%
  group_by(Year) %>%
  summarise(Ratio = mean(Ratio), .groups = "drop")

# 5. 회귀분석 (연도 평균값 기준)
fit <- lm(Ratio ~ Year, data = df_mean)

intercept <- coef(fit)[1]
slope <- coef(fit)[2]
r2 <- summary(fit)$r.squared

eq_text <- paste0(
  "y = ", sprintf("%.4f", slope), "x + ", sprintf("%.4f", intercept),
  "\nR² = ", sprintf("%.4f", r2)
)

# 6. 그래프 그리기
p <- ggplot(df_long, aes(x = Year, y = Ratio)) +
  
  # 회귀선 + 신뢰구간 (연도 평균값 기준)
  geom_smooth(
    data = df_mean,
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "grey80",
    alpha = 0.25,
    linewidth = 0.5
  ) +
  
  # 전체 샘플 점 표시
  geom_point(
    aes(fill = factor(Year)),
    shape = 21,
    size = 3.8,
    color = "black",
    stroke = 0.2
  ) +
  
  # 연도별 색상 적용
  scale_fill_manual(values = year_colors) +
  
  # x축 5년 간격
  scale_x_continuous(
    breaks = seq(1995, 2025, by = 5),
    limits = c(1994, 2023)
  ) +
  
  # 축 라벨
  labs(
    x = "Year",
    y = "13C/12C"
  ) +
  
  # 회귀식, R² 표시
  annotate(
    "text",
    x = 2010,
    y = max(df_long$Ratio),
    label = eq_text,
    hjust = 0,
    vjust = 1,
    size = 3.5,
    color = "black"
  ) +
  
  # 테마 설정
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.18, "cm"),
    plot.margin = margin(10, 30, 10, 10)
  )

# 7. 그래프 출력
print(p)

# 8. 저장
ggsave("figure1B_C13C12.png", p, width = 6, height = 4, dpi = 600, bg = "white")
ggsave("figure1B_C13C12.pdf", p, width = 6, height = 4, device = cairo_pdf, bg = "white")
ggsave("figure1B_C13C12.svg", p, width = 6, height = 4, bg = "white")