# ============================================================
# Figure 3A : v plot only from C sheet
# - Upper panel: v plot
# - Lower panel: count bar plot
# - negative labels: top 10 significant
# - positive labels: all significant (4)
# Input  : figure3B_vplot.xlsx
# Sheet  : "Fig C. Volcano Plot"
# Output : PNG, PDF, SVG, CSV
# ============================================================

# ------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------
library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)
library(svglite)

# ------------------------------------------------------------
# 1) File path
# ------------------------------------------------------------
file_path  <- "C:/Users/pse34/R studio/wine_science/figure3/figure3B_vplot.xlsx"
sheet_name <- "Fig C. Volcano Plot"

# 저장 폴더
out_dir <- getwd()

# ------------------------------------------------------------
# 2) Read Excel
# ------------------------------------------------------------
raw_data <- read_excel(
  file_path,
  sheet = sheet_name,
  col_names = FALSE
)

# ------------------------------------------------------------
# 3) Extract metadata and build dataframe
# ------------------------------------------------------------
years <- raw_data[1, -1] %>%
  unlist(use.names = FALSE) %>%
  as.character() %>%
  as.numeric()

sample_ids <- raw_data[2, -1] %>%
  unlist(use.names = FALSE) %>%
  as.character()

metabolites <- raw_data[-c(1, 2), 1] %>%
  unlist(use.names = FALSE) %>%
  as.character()

intensity_mat <- raw_data[-c(1, 2), -1] %>%
  mutate(across(everything(), ~ as.numeric(.x)))

df_t <- as.data.frame(t(as.matrix(intensity_mat)))
colnames(df_t) <- make.names(metabolites, unique = TRUE)

df_t <- df_t %>%
  mutate(
    Sample_ID = sample_ids,
    Year = years
  ) %>%
  relocate(Sample_ID, Year)

# ------------------------------------------------------------
# 4) Colors
# ------------------------------------------------------------
col_positive <- "#D73027"
col_negative <- "#2C7BB6"
col_neutral  <- "#BFBFBF"
col_bar_grey <- "#C9C9C9"

# ------------------------------------------------------------
# 5) Main function
# ------------------------------------------------------------
generate_vplot_only <- function(data,
                                fdr_cutoff = 0.05,
                                neg_label_n = 10,
                                x_pad_mult = 0.04) {
  
  metabo_cols <- setdiff(names(data), c("Sample_ID", "Year"))
  
  # ----------------------------------------------------------
  # 5-1) Pearson correlation
  # ----------------------------------------------------------
  results <- map_df(metabo_cols, function(metabo) {
    
    x <- data$Year
    y <- data[[metabo]]
    
    valid_idx <- is.finite(x) & is.finite(y)
    
    if (sum(valid_idx) < 3) return(NULL)
    if (var(y[valid_idx], na.rm = TRUE) == 0) return(NULL)
    
    test <- cor.test(
      x[valid_idx],
      y[valid_idx],
      method = "pearson"
    )
    
    tibble(
      Metabolite = metabo,
      r = as.numeric(test$estimate),
      p_value = test$p.value
    )
  })
  
  if (nrow(results) == 0) {
    stop("No valid metabolite data found.")
  }
  
  # ----------------------------------------------------------
  # 5-2) FDR, significance, direction
  # ----------------------------------------------------------
  results <- results %>%
    mutate(
      FDR = p.adjust(p_value, method = "BH"),
      log_FDR = -log10(FDR),
      Significant = FDR < fdr_cutoff,
      Direction = case_when(
        Significant & r > 0 ~ "Positive",
        Significant & r < 0 ~ "Negative",
        TRUE ~ "Not significant"
      )
    ) %>%
    arrange(FDR)
  
  # Inf 처리
  finite_max <- max(results$log_FDR[is.finite(results$log_FDR)], na.rm = TRUE)
  
  results <- results %>%
    mutate(
      log_FDR = ifelse(is.infinite(log_FDR), finite_max, log_FDR)
    )
  
  # ----------------------------------------------------------
  # 5-3) Label selection
  #      negative: top 10 significant
  #      positive: all significant
  # ----------------------------------------------------------
  neg_labels <- results %>%
    filter(Significant, Direction == "Negative") %>%
    arrange(FDR) %>%
    slice_head(n = neg_label_n) %>%
    pull(Metabolite)
  
  pos_labels <- results %>%
    filter(Significant, Direction == "Positive") %>%
    arrange(FDR) %>%
    pull(Metabolite)
  
  results <- results %>%
    mutate(
      Label_Text = case_when(
        Metabolite %in% neg_labels ~ Metabolite,
        Metabolite %in% pos_labels ~ Metabolite,
        TRUE ~ ""
      )
    )
  
  # 통계 결과 저장
  write.csv(
    results,
    file = file.path(out_dir, "Figure3A_vplot_statistics.csv"),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # 5-4) Upper v plot
  # ----------------------------------------------------------
  p_volcano <- ggplot(results, aes(x = r, y = log_FDR)) +
    geom_point(
      aes(color = Direction),
      alpha = 0.9,
      size = 1.8
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey65",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed",
      color = "grey65",
      linewidth = 0.4
    ) +
    
    # negative labels
    geom_text_repel(
      data = results %>% filter(Label_Text != "", Direction == "Negative"),
      aes(label = Label_Text),
      color = col_negative,
      size = 2.6,
      box.padding = 0.35,
      point.padding = 0.18,
      min.segment.length = 0,
      max.overlaps = Inf,
      segment.color = "grey55",
      segment.size = 0.22,
      direction = "both",
      show.legend = FALSE
    ) +
    
    # positive labels
    geom_text_repel(
      data = results %>% filter(Label_Text != "", Direction == "Positive"),
      aes(label = Label_Text),
      color = col_positive,
      size = 2.8,
      box.padding = 0.40,
      point.padding = 0.20,
      min.segment.length = 0,
      max.overlaps = Inf,
      segment.color = "grey55",
      segment.size = 0.22,
      direction = "both",
      show.legend = FALSE
    ) +
    
    scale_color_manual(
      values = c(
        "Positive" = col_positive,
        "Negative" = col_negative,
        "Not significant" = col_neutral
      )
    ) +
    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.5),
      expand = expansion(mult = c(x_pad_mult, x_pad_mult))
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_line(color = "black", linewidth = 0.45),
      axis.text.y  = element_text(size = 11, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      plot.margin  = margin(t = 5, r = 8, b = 0, l = 5)
    ) +
    labs(
      y = expression(-log[10](FDR))
    )
  
  # ----------------------------------------------------------
  # 5-5) Lower count plot data
  #      0.1 bin interval
  # ----------------------------------------------------------
  bin_breaks <- seq(-1, 1, by = 0.1)
  bin_width  <- 0.1
  
  hist_raw <- results %>%
    mutate(
      bin = cut(
        r,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE,
        labels = FALSE
      )
    ) %>%
    filter(!is.na(bin))
  
  bin_df <- tibble(
    bin  = 1:(length(bin_breaks) - 1),
    xmin = bin_breaks[-length(bin_breaks)],
    xmax = bin_breaks[-1]
  ) %>%
    mutate(
      xmid = (xmin + xmax) / 2
    )
  
  bin_count <- hist_raw %>%
    count(bin, name = "count")
  
  bin_sig <- hist_raw %>%
    group_by(bin) %>%
    summarise(
      sig_neg = sum(Significant & r < 0, na.rm = TRUE),
      sig_pos = sum(Significant & r > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  bin_df2 <- bin_df %>%
    left_join(bin_count, by = "bin") %>%
    left_join(bin_sig, by = "bin") %>%
    mutate(
      count   = replace_na(count, 0),
      sig_neg = replace_na(sig_neg, 0),
      sig_pos = replace_na(sig_pos, 0)
    )
  
  # ----------------------------------------------------------
  # 5-6) Color bins
  # negative: -1.0 ~ -0.9 only
  # positive: pure extreme bins only, stop before mixed bin
  # ----------------------------------------------------------
  neg_bins_to_color <- bin_df2 %>%
    filter(xmin == -1.0, xmax == -0.9) %>%
    pull(bin)
  
  pos_info <- bin_df2 %>%
    arrange(desc(xmin)) %>%
    mutate(
      pure_sig_pos = (count > 0 & sig_pos == count & sig_pos > 0)
    )
  
  pos_bins_to_color <- c()
  
  for (i in seq_len(nrow(pos_info))) {
    
    if (pos_info$xmin[i] < 0) next
    if (pos_info$count[i] == 0) next
    
    if (pos_info$pure_sig_pos[i]) {
      pos_bins_to_color <- c(pos_bins_to_color, pos_info$bin[i])
    } else {
      if (pos_info$sig_pos[i] > 0) break
    }
  }
  
  bin_df2 <- bin_df2 %>%
    mutate(
      fill_group = case_when(
        bin %in% neg_bins_to_color ~ "Negative",
        bin %in% pos_bins_to_color ~ "Positive",
        TRUE ~ "Other"
      )
    )
  
  # ----------------------------------------------------------
  # 5-7) Lower count bar plot
  # ----------------------------------------------------------
  p_hist <- ggplot(bin_df2, aes(x = xmid, y = -count)) +
    geom_hline(
      yintercept = 0,
      color = "black",
      linewidth = 0.45
    ) +
    geom_col(
      aes(fill = fill_group),
      width = bin_width * 0.92,
      color = "black",
      linewidth = 0.25,
      alpha = 0.95
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey65",
      linewidth = 0.4
    ) +
    scale_fill_manual(
      values = c(
        "Negative" = col_negative,
        "Positive" = col_positive,
        "Other"    = col_bar_grey
      )
    ) +
    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.5),
      expand = expansion(mult = c(x_pad_mult, x_pad_mult))
    ) +
    scale_y_continuous(labels = abs) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text  = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.margin = margin(t = 0, r = 8, b = 5, l = 5)
    ) +
    labs(
      x = "Pearson's r with year",
      y = "Count"
    )
  
  # ----------------------------------------------------------
  # 5-8) Combine panels
  # ----------------------------------------------------------
  final_plot <- p_volcano / plot_spacer() / p_hist +
    plot_layout(heights = c(2.30, 0.05, 1))
  
  # ----------------------------------------------------------
  # 5-9) Save
  # ----------------------------------------------------------
  ggsave(
    filename = file.path(out_dir, "Figure3B_vplot_only.png"),
    plot = final_plot,
    width = 6.0,
    height = 5.5,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(out_dir, "Figure3B_vplot_only.pdf"),
    plot = final_plot,
    width = 6.0,
    height = 5.5,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(out_dir, "Figure3B_vplot_only.svg"),
    plot = final_plot,
    width = 6.0,
    height = 5.5,
    bg = "white"
  )
  
  message("Saved:")
  message(" - Figure3B_vplot_only.png")
  message(" - Figure3B_vplot_only.pdf")
  message(" - Figure3B_vplot_only.svg")
  message(" - Figure3B_vplot_statistics.csv")
  
  return(final_plot)
}

# ------------------------------------------------------------
# 6) Run
# ------------------------------------------------------------
fig3B_vplot <- generate_vplot_only(
  data = df_t,
  fdr_cutoff = 0.05,
  neg_label_n = 10,
  x_pad_mult = 0.04
)

fig3B_vplot