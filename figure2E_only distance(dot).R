# =========================================================
# Figure 1E distance dot plot
# Platform-wise min-max range + yearly mean dots
# =========================================================

# -----------------------------
# 0) Packages
# -----------------------------
packages <- c("tidyverse", "ggplot2", "svglite")

installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) install.packages(p)
}
invisible(lapply(packages, library, character.only = TRUE))

# -----------------------------
# 1) File names
# -----------------------------
file_gcms      <- "figure1E_GCMS.csv"
file_headspace <- "figure1E_headspaceGCMS.csv"
file_lcn       <- "figure1E_LCn.csv"
file_lcp       <- "figure1E_LCp.csv"

required_files <- c(file_gcms, file_headspace, file_lcn, file_lcp)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "These files were not found in the current working directory:\n",
    paste("-", missing_files, collapse = "\n")
  )
}

# -----------------------------
# 2) Year colors
# -----------------------------
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

year_levels <- names(year_colors)

# -----------------------------
# 3) Reader
# -----------------------------
read_omics_data <- function(input_file) {
  df <- read.csv(input_file, check.names = FALSE)
  
  if (ncol(df) < 3) {
    stop("File has fewer than 3 columns: ", input_file)
  }
  
  meta <- df[, 1:2, drop = FALSE]
  colnames(meta) <- c("SampleID", "Year")
  meta$SampleID <- as.character(meta$SampleID)
  meta$Year <- as.integer(as.character(meta$Year))
  
  x <- df[, -(1:2), drop = FALSE] %>%
    mutate(across(everything(), as.numeric))
  
  x <- x %>%
    mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
  
  keep_cols <- !sapply(x, function(z) all(is.na(z)))
  x <- x[, keep_cols, drop = FALSE]
  
  var_vec <- sapply(x, function(z) var(z, na.rm = TRUE))
  keep_cols2 <- !is.na(var_vec) & var_vec > 0
  x <- x[, keep_cols2, drop = FALSE]
  
  if (ncol(x) == 0) {
    stop("No usable feature columns after filtering in file: ", input_file)
  }
  
  x_mat <- as.data.frame(x)
  rownames(x_mat) <- paste0("Y", meta$Year, "_", meta$SampleID)
  
  list(meta = meta, x = x_mat)
}

# -----------------------------
# 4) PCA helper
# -----------------------------
run_pca <- function(x_mat, meta) {
  pca_res <- prcomp(x_mat, center = TRUE, scale. = TRUE)
  
  scores <- as.data.frame(pca_res$x[, 1:2, drop = FALSE]) %>%
    setNames(c("PC1", "PC2")) %>%
    bind_cols(meta)
  
  scores
}

# -----------------------------
# 5) Absolute PC1 distance from 1995
# -----------------------------
calc_pc1_distance_from_1995 <- function(input_file, platform_name) {
  dat <- read_omics_data(input_file)
  scores <- run_pca(dat$x, dat$meta)
  
  if (!1995 %in% scores$Year) {
    stop("No Year == 1995 samples found in file: ", input_file)
  }
  
  ref_pc1 <- mean(scores$PC1[scores$Year == 1995], na.rm = TRUE)
  
  raw_dist <- scores %>%
    mutate(
      Platform = platform_name,
      PC1_dist = abs(PC1 - ref_pc1)
    ) %>%
    select(SampleID, Year, Platform, PC1_dist)
  
  # 연도별 평균값 (dot용)
  mean_dist <- raw_dist %>%
    group_by(Platform, Year) %>%
    summarise(
      MeanDist = mean(PC1_dist, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 플랫폼별 전체 min/max (range용)
  range_dist <- mean_dist %>%
    group_by(Platform) %>%
    summarise(
      xmin = min(MeanDist, na.rm = TRUE),
      xmax = max(MeanDist, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(raw = raw_dist, mean = mean_dist, range = range_dist)
}

# -----------------------------
# 6) Run for 4 platforms
# -----------------------------
dist_gcms <- calc_pc1_distance_from_1995(file_gcms,      "GC–MS")
dist_hs   <- calc_pc1_distance_from_1995(file_headspace, "Headspace")
dist_lcn  <- calc_pc1_distance_from_1995(file_lcn,       "LC negative")
dist_lcp  <- calc_pc1_distance_from_1995(file_lcp,       "LC positive")

dist_mean_all <- bind_rows(
  dist_gcms$mean,
  dist_hs$mean,
  dist_lcn$mean,
  dist_lcp$mean
)

dist_range_all <- bind_rows(
  dist_gcms$range,
  dist_hs$range,
  dist_lcn$range,
  dist_lcp$range
)

platform_levels <- c(
  "Headspace",
  "GC–MS",
  "LC positive",
  "LC negative"
)

dist_mean_all$Platform <- factor(dist_mean_all$Platform, levels = platform_levels)
dist_mean_all$Year <- factor(as.character(dist_mean_all$Year), levels = year_levels)

dist_range_all$Platform <- factor(dist_range_all$Platform, levels = platform_levels)

# 라벨 위치용
range_labels <- dist_range_all %>%
  mutate(
    xmin_lab = round(xmin, 1),
    xmax_lab = round(xmax, 1)
  )

# -----------------------------
# 7) Dot plot
# -----------------------------
p_dot <- ggplot() +
  # 옅은 회색 min-max range
  geom_segment(
    data = dist_range_all,
    aes(
      x = xmin,
      xend = xmax,
      y = Platform,
      yend = Platform
    ),
    color = "grey85",
    linewidth = 5.2,
    lineend = "round"
  ) +
  
  # 연도별 평균 dot 12개
  geom_point(
    data = dist_mean_all,
    aes(
      x = MeanDist,
      y = Platform,
      fill = Year
    ),
    shape = 21,
    size = 5.2,
    color = "grey25",
    stroke = 0.4
  ) +
  
  # 최소값 라벨
  geom_text(
    data = range_labels,
    aes(
      x = xmin,
      y = Platform,
      label = xmin_lab
    ),
    hjust = 1.35,
    vjust = -0.85,
    size = 3.0,
    color = "grey40"
  ) +
  
  # 최대값 라벨
  geom_text(
    data = range_labels,
    aes(
      x = xmax,
      y = Platform,
      label = xmax_lab
    ),
    hjust = -0.35,
    vjust = -0.85,
    size = 3.0,
    color = "grey40"
  ) +
  
  scale_fill_manual(values = year_colors, drop = FALSE) +
  
  scale_x_continuous(
    breaks = seq(0, 25, by = 5),
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  
  scale_y_discrete(
    limits = rev(platform_levels)
  ) +
  
  labs(
    x = "Absolute PC1 distance from 1995",
    y = NULL
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 10) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8),
    
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.line.x = element_line(linewidth = 0.45),
    axis.ticks.x = element_line(linewidth = 0.35),
    
    # 핵심 추가
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(8, 18, 8, 18, "mm")
  )

print(p_dot)

# -----------------------------
# 8) Save
# -----------------------------
ggsave(
  "figure1E_distance_dotplot.png",
  p_dot,
  width = 6.0,
  height = 2.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  "figure1E_distance_dotplot.pdf",
  p_dot,
  width = 6.0,
  height = 2.8,
  bg = "white"
)

ggsave(
  "figure1E_distance_dotplot.svg",
  p_dot,
  width = 6.0,
  height = 2.8,
  bg = "white"
)

# -----------------------------
# 9) Export data
# -----------------------------
write.csv(dist_mean_all,  "figure1E_distance_dotplot_yearly_mean.csv", row.names = FALSE)
write.csv(dist_range_all, "figure1E_distance_dotplot_range.csv", row.names = FALSE)

# -----------------------------
# 10) Show
# -----------------------------
print(p_dot)

cat("\nDone.\n")
cat("Saved files:\n")
cat("- figure1E_distance_dotplot.png / pdf / svg\n")
cat("- figure1E_distance_dotplot_yearly_mean.csv\n")
cat("- figure1E_distance_dotplot_range.csv\n")