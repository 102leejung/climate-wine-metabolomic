## -----------------------------------------------------------------------------------------------------------
library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================================
# [Step 1] 데이터 로드 및 카테고리 태그 추가
# ==========================================
file_name <- "figure1G_enose etongue.xlsx"

# 시트별 로드 및 카테고리 구분 컬럼(Type) 생성
enose_data <- read_excel(file_name, sheet = "E-nose") %>% 
  mutate(Type = "E-nose")

etongue_data <- read_excel(file_name, sheet = "E-tongue") %>% 
  mutate(Type = "E-tongue")

# 컬럼명 통일 및 병합
df_total <- bind_rows(enose_data, etongue_data) %>%
  rename_with(~ "R_sensor", any_of(c("R_value (E-tongue)", "R_value (E-nose)")))

# 유의성 마커 및 절댓값 생성
df_total <- df_total %>%
  mutate(
    sig_label = case_when(
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    abs_R_sensor = abs(R_sensor),
    abs_R_year = abs(`R_value (year)`)
  )

# 유의미한 데이터만 필터링
sig_metabolites <- df_total %>% 
  filter(P_value < 0.05) %>% 
  pull(Metabolite) %>% 
  unique()

df_plot <- df_total %>% filter(Metabolite %in% sig_metabolites)

# ==========================================
# [Step 2] 시각화 (X축 회전 및 카테고리 분리)
# ==========================================

# 공통 테마 수정: X축 텍스트 회전 추가
my_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    # X축 텍스트를 세로로(90도) 세우고 정렬 조정
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
    strip.background = element_rect(fill = "grey95", color = NA), # 카테고리 라벨 배경
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# [Panel A] 좌측: 카테고리 분리형 센서 상관도
plot_sensor <- ggplot(df_plot, aes(x = Sensor, y = Metabolite)) +
  # 카테고리(Type)별로 구역 분리, 센서 개수에 따라 폭 자동 조절(space="free_x")
  facet_grid(. ~ Type, scales = "free_x", space = "free_x") +
  geom_point(data = filter(df_plot, sig_label == ""), color = "grey80", size = 1.5) +
  geom_point(data = filter(df_plot, sig_label != ""), 
             aes(size = abs_R_sensor, color = R_sensor), alpha = 0.85) +
  geom_text(aes(label = sig_label), vjust = 0.7, color = "white", size = 3, fontface = "bold") +
  scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026", midpoint = 0, name = "Sensor R") +
  scale_size_continuous(range = c(3, 9), guide = "none") +
  my_theme

# [Panel B] 우측: 연도 상관도
df_year <- df_plot %>% 
  select(Metabolite, `R_value (year)`, abs_R_year) %>% 
  distinct() %>%
  mutate(sig_year = ifelse(abs_R_year > 0.5, "**", ifelse(abs_R_year > 0.3, "*", ""))) 

plot_year <- ggplot(df_year, aes(x = "Vintage Year", y = Metabolite)) +
  geom_point(aes(size = abs_R_year, color = `R_value (year)`), alpha = 0.85) +
  geom_text(aes(label = sig_year), vjust = 0.7, color = "white", size = 3, fontface = "bold") +
  scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026", midpoint = 0, name = "Year R") +
  scale_size_continuous(range = c(3, 9), guide = "none") +
  my_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5) # 연도축은 회전 불필요
  )

# ==========================================
# [Step 3] 최종 병합 및 출력
# ==========================================
final_plot <- plot_sensor + plot_year + 
  plot_layout(widths = c(4, 1), guides = "collect") & 
  theme(legend.position = "bottom")

print(final_plot)

ggsave("figure1G_Combined_Correlation_Plot2.png", plot = final_plot, width = 14, height = 10, dpi = 300)


## -----------------------------------------------------------------------------------------------------------
library(readxl)
library(ggplot2)
library(dplyr)

# ==========================================
# [Step 1] 데이터 로드 및 전처리
# ==========================================
file_name <- "figure1G_enose etongue.xlsx"

enose_data <- read_excel(file_name, sheet = "E-nose") %>% 
  mutate(Type = "E-nose")

etongue_data <- read_excel(file_name, sheet = "E-tongue") %>% 
  mutate(Type = "E-tongue")

df_total <- bind_rows(enose_data, etongue_data) %>%
  rename_with(~ "R_sensor", any_of(c("R_value (E-tongue)", "R_value (E-nose)"))) %>%
  mutate(
    abs_R_year = abs(`R_value (year)`)
  )

# 유의한 대사산물만 추출
sig_metabolites <- df_total %>% 
  filter(P_value < 0.05) %>% 
  pull(Metabolite) %>% 
  unique()

df_plot <- df_total %>% 
  filter(Metabolite %in% sig_metabolites)

# ==========================================
# [Step 2] 대사산물 순서 정리
#   1. E-tongue에서 유의한 대사산물 우선
#   2. 그 안에서 abs_R_year 최대값 큰 순서
# ==========================================
metab_order <- df_plot %>%
  group_by(Metabolite) %>%
  summarise(
    etongue_sig = any(Type == "E-tongue" & P_value < 0.05),
    max_abs_r_year = max(abs_R_year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(etongue_sig), desc(max_abs_r_year), Metabolite) %>%
  pull(Metabolite)

df_plot <- df_plot %>%
  mutate(
    Metabolite = factor(Metabolite, levels = metab_order),
    Type = factor(Type, levels = c("E-nose", "E-tongue")),
    R_year_group = case_when(
      abs_R_year <= 0.40 ~ "≤0.40",
      abs_R_year <= 0.60 ~ "0.40–0.60",
      abs_R_year <= 0.80 ~ "0.60–0.80",
      TRUE ~ ">0.80"
    )
  )

# size 범례 순서 고정
df_plot$R_year_group <- factor(
  df_plot$R_year_group,
  levels = c("≤0.40", "0.40–0.60", "0.60–0.80", ">0.80")
)

# ==========================================
# [Step 3] 센서 순서 고정
# ==========================================
enose_sensor_order <- c(
  "(S)-2-Methyl-1-butanol",
  "(Z)-1,2-Dichloroethene",
  "2-Methyl-1-propanol",
  "3-Mercapto-4-methyl-2-pentanone",
  "Dill ether",
  "Ethyl hexanoate",
  "Ethyl isovalerate",
  "Methyl salicylate"
)

etongue_sensor_order <- c(
  "AHS", "ANS", "CPS", "CTS", "NMS", "PKS", "SCS"
)

df_plot <- df_plot %>%
  mutate(
    Sensor = case_when(
      Type == "E-nose" ~ as.character(factor(Sensor, levels = enose_sensor_order)),
      Type == "E-tongue" ~ as.character(factor(Sensor, levels = etongue_sensor_order)),
      TRUE ~ as.character(Sensor)
    )
  )

# ==========================================
# [Step 4] 테마
# ==========================================
my_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_text(size = 11),
    panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 10, b = 10),
    plot.margin = margin(t = 10, r = 15, b = 30, l = 10)
  )

# ==========================================
# [Step 5] 그래프
#   같은 metabolite 기준으로
#   위 = E-nose / 아래 = E-tongue
# ==========================================
final_plot <- ggplot(df_plot, aes(x = Metabolite, y = Sensor)) +
  facet_grid(Type ~ ., scales = "free_y", space = "free_y") +
  
  # 비유의 조합은 배경점처럼 연하게 표시
  geom_point(
    data = filter(df_plot, P_value >= 0.05),
    color = "grey85",
    size = 1.3
  ) +
  
  # 유의한 조합만 색 + 4단계 크기 적용
  geom_point(
  data = filter(df_plot, P_value < 0.05),
  aes(size = R_year_group, color = R_sensor),
  alpha = 0.7
) +
  
  scale_color_gradient2(
    low = "#053061",
    mid = "white",
    high = "#67001f",
    midpoint = 0,
    name = "Sensor R (Color)"
  ) +
  
  scale_size_manual(
    values = c(
      "≤0.40" = 2.5,
      "0.40–0.60" = 5,
      "0.60–0.80" = 7.5,
      ">0.80" = 10
    ),
    name = "Year |R|"
  ) +
  
  scale_x_discrete(expand = expansion(add = c(0.3, 0.3))) +
  my_theme

# ==========================================
# [Step 6] 출력 및 저장
# ==========================================
print(final_plot)

ggsave(
  "figure1G_full_ordered_size4group.png",
  plot = final_plot,
  width = 24,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "figure1G_full_ordered_size4group.svg",
  plot = final_plot,
  width = 24,
  height = 8,
  bg = "white"
)

