# =========================================================
# Figure 1E
# PCA score plots only (1 x 4)
# Order: Headspace -> GC-MS -> LC positive -> LC negative
# =========================================================

# -----------------------------
# 0) Packages
# -----------------------------
packages <- c("tidyverse", "ggplot2", "cowplot", "grid", "svglite")

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
# 2) Colors
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
# Assumed format:
# 1st col = SampleID
# 2nd col = Year
# 3rd+ cols = features
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

  list(meta = meta, x = x_mat, raw = df)
}

# -----------------------------
# 4) PCA helper
# -----------------------------
run_pca <- function(x_mat, meta) {
  pca_res <- prcomp(x_mat, center = TRUE, scale. = TRUE)

  var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

  scores <- as.data.frame(pca_res$x[, 1:2, drop = FALSE]) %>%
    setNames(c("PC1", "PC2")) %>%
    bind_cols(meta)

  list(
    scores = scores,
    pc1_var = round(var_explained[1] * 100, 1),
    pc2_var = round(var_explained[2] * 100, 1)
  )
}

# -----------------------------
# 5) Ellipse helper
# -----------------------------
make_ellipse <- function(x, y, level = 0.95, npoints = 200) {
  center <- c(mean(x), mean(y))
  cov_mat <- cov(cbind(x, y))
  radius <- sqrt(qchisq(level, df = 2))

  angles <- seq(0, 2 * pi, length.out = npoints)
  circle <- cbind(cos(angles), sin(angles))

  eig <- eigen(cov_mat)
  ellipse <- t(center + radius * t(circle %*% diag(sqrt(eig$values)) %*% t(eig$vectors)))

  ellipse_df <- as.data.frame(ellipse)
  colnames(ellipse_df) <- c("x", "y")
  ellipse_df
}

# -----------------------------
# 6) Run PCA for each platform
# -----------------------------
run_platform_pca <- function(input_file, platform_name) {
  dat <- read_omics_data(input_file)
  res <- run_pca(dat$x, dat$meta)

  scores <- res$scores %>%
    mutate(
      Platform = platform_name,
      Year = factor(as.character(Year), levels = year_levels)
    )

  ellipse_df <- make_ellipse(scores$PC1, scores$PC2, level = 0.95)

  list(
    platform = platform_name,
    scores = scores,
    ellipse = ellipse_df,
    pc1_var = res$pc1_var,
    pc2_var = res$pc2_var
  )
}

pca_hs   <- run_platform_pca(file_headspace, "Headspace")
pca_gcms <- run_platform_pca(file_gcms, "GC–MS")
pca_lcp  <- run_platform_pca(file_lcp, "LC positive")
pca_lcn  <- run_platform_pca(file_lcn, "LC negative")

all_scores <- bind_rows(
  pca_hs$scores,
  pca_gcms$scores,
  pca_lcp$scores,
  pca_lcn$scores
)

# -----------------------------
# 7) Axis label helper
# -----------------------------
fmt_pct <- function(x) {
  if (abs(x - round(x)) < 1e-9) {
    as.character(round(x))
  } else {
    sprintf("%.1f", x)
  }
}

# -----------------------------
# 8) PCA plot function
# -----------------------------
make_pca_plot <- function(scores_df, ellipse_df, pc1_var, pc2_var,
                          title_text, xlim_use, ylim_use,
                          show_legend = FALSE) {
  ggplot(scores_df, aes(x = PC1, y = PC2)) +
    geom_path(
      data = ellipse_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      color = "grey55",
      linewidth = 0.4
    ) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_point(
      aes(fill = Year),
      shape = 21,
      size = 2.3,
      stroke = 0.4,
      color = "black"
    ) +
    scale_fill_manual(
      values = year_colors,
      breaks = year_levels,
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = xlim_use,
      breaks = pretty(xlim_use, n = 5)
    ) +
    scale_y_continuous(
      limits = ylim_use,
      breaks = pretty(ylim_use, n = 5)
    ) +
    labs(
      title = title_text,
      x = paste0("PC1(", fmt_pct(pc1_var), "%)"),
      y = paste0("PC2(", fmt_pct(pc2_var), "%)")
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.title = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 7),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      plot.margin = margin(2, 2, 2, 2, "mm"),
      aspect.ratio = 0.65
    )
}

# -----------------------------
# 9) Manual panel limits
# -----------------------------
# Headspace
xlim_hs  <- c(-8, 11)
ylim_hs  <- c(-6.5, 6.5)

# GC-MS
xlim_gcms <- c(-15, 17)
ylim_gcms <- c(-11, 11)

# LC positive
xlim_lcp <- c(-16, 16)
ylim_lcp <- c(-6.5, 6.5)

# LC negative
xlim_lcn <- c(-17, 19)
ylim_lcn <- c(-14, 14)

# -----------------------------
# 10) Build plots
# -----------------------------
p_hs <- make_pca_plot(
  pca_hs$scores, pca_hs$ellipse,
  pca_hs$pc1_var, pca_hs$pc2_var,
  "Headspace",
  xlim_use = xlim_hs,
  ylim_use = ylim_hs,
  show_legend = FALSE
)

p_gcms <- make_pca_plot(
  pca_gcms$scores, pca_gcms$ellipse,
  pca_gcms$pc1_var, pca_gcms$pc2_var,
  "GC–MS",
  xlim_use = xlim_gcms,
  ylim_use = ylim_gcms,
  show_legend = FALSE
)

p_lcp_with_legend <- make_pca_plot(
  pca_lcp$scores, pca_lcp$ellipse,
  pca_lcp$pc1_var, pca_lcp$pc2_var,
  "LC positive",
  xlim_use = xlim_lcp,
  ylim_use = ylim_lcp,
  show_legend = TRUE
)

p_lcn <- make_pca_plot(
  pca_lcn$scores, pca_lcn$ellipse,
  pca_lcn$pc1_var, pca_lcn$pc2_var,
  "LC negative",
  xlim_use = xlim_lcn,
  ylim_use = ylim_lcn,
  show_legend = FALSE
)

year_legend <- cowplot::get_legend(p_lcp_with_legend)
p_lcp <- p_lcp_with_legend + theme(legend.position = "none")

# -----------------------------
# 11) PCA row only
# -----------------------------
pca_row <- cowplot::plot_grid(
  p_hs, p_gcms, p_lcp, p_lcn,
  ncol = 4,
  align = "h"
)

p_final <- cowplot::plot_grid(
  pca_row,
  year_legend,
  ncol = 2,
  rel_widths = c(1, 0.18)
)

print(p_final)

# -----------------------------
# 12) Save
# -----------------------------
ggsave(
  "figure1E_PCA_only.svg",
  p_final,
  width = 7.2,
  height = 2.0,
  bg = "white"
)

ggsave(
  "figure1E_PCA_only.png",
  p_final,
  width = 7.2,
  height = 2.0,
  dpi = 600,
  bg = "white"
)

ggsave(
  "figure1E_PCA_only.pdf",
  p_final,
  width = 7.2,
  height = 2.0,
  bg = "white"
)

# -----------------------------
# 13) Export data
# -----------------------------
write.csv(
  all_scores,
  "figure1E_all_pca_scores.csv",
  row.names = FALSE
)

# -----------------------------
# 14) Show
# -----------------------------
print(p_final)

cat("\nDone.\n")
cat("Saved files:\n")
cat("- figure1E_PCA_only.png / pdf / svg\n")
cat("- figure1E_all_pca_scores.csv\n")