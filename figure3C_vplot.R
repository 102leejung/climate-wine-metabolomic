library(tidyverse)
library(patchwork)
library(ggrepel)
library(readxl)

# =========================================================
# Figure 2C
# volcano plot + mirrored histogram 
# based on figure2C_vplot.xlsx
# =========================================================

# =========================================================
# 0. Packages
# =========================================================
library(tidyverse)
library(patchwork)
library(ggrepel)
library(readxl)

# =========================================================
# 1. Read file
# =========================================================
df <- read_excel("figure2C_vplot.xlsx")

# remove duplicated column names if any
df <- df[, !duplicated(colnames(df))]

# =========================================================
# 2. Select region
# =========================================================
region_name <- "All"
# region_name <- "Europe"
# region_name <- "North America"
# region_name <- "Australia"
# region_name <- "South America"

# =========================================================
# 3. Subset data by region
# =========================================================
if (region_name != "All") {
  plot_data <- df[df$Continent == region_name, ]
} else {
  plot_data <- df
}

# optional check
dim(plot_data)
head(plot_data[, 1:min(8, ncol(plot_data))])

# =========================================================
# 4. Set metabolite columns
#    assume columns 1~5 are metadata
#    and columns 6:end are metabolite variables
# =========================================================
metabo_cols <- names(plot_data)[6:ncol(plot_data)]

# optional check
length(metabo_cols)
head(metabo_cols)

# =========================================================
# 5. Pearson correlation analysis
#    HI_7mo vs each metabolite
# =========================================================
results <- map_df(metabo_cols, function(metabo) {
  
  x <- plot_data[["HI_7mo"]]
  y <- as.numeric(plot_data[[metabo]])
  
  valid_pairs <- sum(is.finite(x) & is.finite(y))
  if (valid_pairs < 3) return(NULL)
  
  if (var(y, na.rm = TRUE) == 0) return(NULL)
  
  test <- cor.test(x, y, method = "pearson")
  
  tibble(
    Metabolite = metabo,
    r = as.numeric(test$estimate),
    p_value = test$p.value
  )
})

# optional check
dim(results)
head(results)
summary(results$r)

# =========================================================
# 6. FDR and significance
# =========================================================
results <- results %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    log_FDR = -log10(FDR),
    Significant = FDR < 0.05,
    Direction = case_when(
      Significant & r > 0 ~ "Positive",
      Significant & r < 0 ~ "Negative",
      TRUE ~ "Not Significant"
    )
  ) %>%
  arrange(FDR) %>%
  mutate(
    Label_Text = ifelse(row_number() <= 15 & Significant, Metabolite, "")
  )

# optional check
head(results, 20)
table(results$Direction)
sum(results$Significant)

# =========================================================
# 7. Colors
# =========================================================
col_positive <- "red"
col_negative <- "#0072B2"
col_neutral  <- "#B0B0B0"

# =========================================================
# 8. Common theme
# =========================================================
theme_fig2c <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9, color = "black"),
    plot.margin = margin(5, 5, 5, 5)
  )

# =========================================================
# 9. Volcano plot
# =========================================================
volcano_plot <- ggplot(results, aes(x = r, y = log_FDR, color = Direction)) +
  geom_point(
    size = 1.8,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  geom_text_repel(
    data = filter(results, Label_Text != ""),
    aes(label = Label_Text, color = Direction),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Positive" = col_positive,
      "Negative" = col_negative,
      "Not Significant" = col_neutral
    )
  ) +
  scale_x_continuous(limits = c(-1.0, 1.0)) +
  labs(
    title = paste("HI_7mo -", region_name),
    x = NULL,
    y = expression(-log[10](FDR))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(b = 0)
  )

# volcano 확인
volcano_plot

# =========================================================
# 10. Mirrored histogram
#     grey = all metabolites
#     blue = significant negative
#     red  = significant positive
# =========================================================
results_hist <- results %>%
  filter(is.finite(r), r >= -1, r <= 1)

hist_plot <- ggplot() +
  # 전체 분포 (회색)
  geom_histogram(
    data = results_hist,
    aes(x = r, y = -after_stat(count)),
    binwidth = 0.1,
    boundary = -1,
    fill = col_neutral,
    color = "black",
    linewidth = 0.2
  ) +
  
  # 유의한 negative (파랑)
  geom_histogram(
    data = results_hist %>% filter(Significant, Direction == "Negative"),
    aes(x = r, y = -after_stat(count)),
    binwidth = 0.1,
    boundary = -1,
    fill = col_negative,
    color = "black",
    linewidth = 0.2
  ) +
  
  # 유의한 positive (빨강)
  geom_histogram(
    data = results_hist %>% filter(Significant, Direction == "Positive"),
    aes(x = r, y = -after_stat(count)),
    binwidth = 0.1,
    boundary = -1,
    fill = col_positive,
    color = "black",
    linewidth = 0.2
  ) +
  
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  scale_y_continuous(labels = abs) +
  
  labs(
    x = "Pearson correlation (r)",
    y = "Count"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
  )

# histogram 확인
hist_plot

# =========================================================
# 11. Combine volcano + histogram
# =========================================================
Figure2C_main <- volcano_plot / plot_spacer() / hist_plot +
  plot_layout(heights = c(3, 0.03, 0.55))

Figure2C_main


#-----------
# v plot+histgram 파일 저장
ggsave(
  filename = "Figure2C_vplot_histogram.svg",
  plot = Figure2C_main,
  width = 7,
  height = 5,
  bg = "white"
)
ggsave(
  filename = "Figure2C_vplot_histogram.pdf",
  plot = Figure2C_main,
  width = 7,
  height = 5,
  bg = "white"
)