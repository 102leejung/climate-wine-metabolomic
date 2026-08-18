# ============================================================
# Figure 3 line graph
#
# Negative: significant metabolites 중 FDR 기준 상위 10개
# Positive: 유의한 positive metabolites 전체
#
# Input  : figure3A_vplot.xlsx
# Sheet  : Fig D. Line Graph
# Output : PNG, PDF, SVG, CSV
# ============================================================


# ------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------
required_packages <- c(
  "tidyverse",
  "readxl",
  "patchwork",
  "ggrepel",
  "svglite"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)
library(svglite)


# ------------------------------------------------------------
# 1) File path
# ------------------------------------------------------------
file_path <- paste0(
  "C:/Users/pse34/R studio/wine_science/",
  "figure3/figure3B_vplot.xlsx"
)

sheet_name <- "Fig D. Line Graph"

# 결과 저장 폴더
out_dir <- "C:/Users/pse34/R studio/wine_science/figure3/output"

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}

if (!file.exists(file_path)) {
  stop(
    paste0(
      "Excel 파일을 찾을 수 없습니다:\n",
      file_path
    )
  )
}

message("Input file: ", file_path)
message("Output directory: ", out_dir)


# ------------------------------------------------------------
# 2) Read Excel
# ------------------------------------------------------------
# Excel 구조
# 1행: Year
# 2행: Sample ID
# 3행 이후: metabolite name + intensity

raw_data <- read_excel(
  path = file_path,
  sheet = sheet_name,
  col_names = FALSE
)

if (nrow(raw_data) < 3) {
  stop("Excel 데이터에 최소 3개 행이 필요합니다.")
}

if (ncol(raw_data) < 2) {
  stop("Excel 데이터에 최소 2개 열이 필요합니다.")
}


# ------------------------------------------------------------
# 3) Extract data
# ------------------------------------------------------------

# 1행: 연도
years <- raw_data[1, -1] %>%
  unlist(use.names = FALSE) %>%
  as.character() %>%
  as.numeric()

# 2행: Sample ID
sample_ids <- raw_data[2, -1] %>%
  unlist(use.names = FALSE) %>%
  as.character()

# 3행 이후 첫 번째 열: 대사체 이름
metabolite_original <- raw_data[-c(1, 2), 1] %>%
  unlist(use.names = FALSE) %>%
  as.character() %>%
  trimws()

# 빈 대사체 이름 제거를 위한 위치
valid_metabolite_rows <- !is.na(metabolite_original) &
  metabolite_original != ""

metabolite_original <- metabolite_original[
  valid_metabolite_rows
]

# 3행 이후 intensity
intensity_mat <- raw_data[-c(1, 2), -1] %>%
  slice(which(valid_metabolite_rows)) %>%
  mutate(
    across(
      everything(),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# 샘플 수 확인
if (length(years) != ncol(intensity_mat)) {
  stop(
    paste0(
      "Year 개수와 intensity sample 열 개수가 다릅니다.\n",
      "Year 개수: ", length(years), "\n",
      "Sample 열 개수: ", ncol(intensity_mat)
    )
  )
}

if (length(sample_ids) != ncol(intensity_mat)) {
  stop(
    paste0(
      "Sample ID 개수와 intensity sample 열 개수가 다릅니다.\n",
      "Sample ID 개수: ", length(sample_ids), "\n",
      "Sample 열 개수: ", ncol(intensity_mat)
    )
  )
}


# ------------------------------------------------------------
# 3-1) Metabolite name mapping
# ------------------------------------------------------------
# 분석용 열 이름은 make.names()로 안전하게 만들되,
# figure에는 원래 대사체 이름을 표시

metabolite_safe <- make.names(
  metabolite_original,
  unique = TRUE
)

metabolite_name_map <- tibble(
  Metabolite = metabolite_safe,
  Display_Name = metabolite_original
)

write.csv(
  metabolite_name_map,
  file = file.path(
    out_dir,
    "Figure3_metabolite_name_mapping.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 3-2) Transpose data
# ------------------------------------------------------------
# 원래 구조: 대사체 = 행, 샘플 = 열
# 변환 구조: 샘플 = 행, 대사체 = 열

df_t <- as.data.frame(
  t(as.matrix(intensity_mat)),
  check.names = FALSE
)

colnames(df_t) <- metabolite_safe

df_t <- df_t %>%
  mutate(
    Sample_ID = sample_ids,
    Year = years
  ) %>%
  relocate(
    Sample_ID,
    Year
  ) %>%
  filter(
    !is.na(Year),
    is.finite(Year)
  )

message("Number of samples: ", nrow(df_t))
message("Number of metabolites: ", length(metabolite_safe))
message(
  "Years: ",
  paste(
    sort(unique(df_t$Year)),
    collapse = ", "
  )
)


# ------------------------------------------------------------
# 4) Colors
# ------------------------------------------------------------
col_positive <- "#D73027"
col_negative <- "#2C7BB6"
col_grey <- "grey75"


# ------------------------------------------------------------
# 5) Main plotting function
# ------------------------------------------------------------
generate_line_graph <- function(
    data,
    name_map,
    output_directory,
    fdr_cutoff = 0.05,
    neg_line_n = 10) {
  
  metabo_cols <- setdiff(
    names(data),
    c("Sample_ID", "Year")
  )
  
  if (length(metabo_cols) == 0) {
    stop("분석 가능한 대사체 열이 없습니다.")
  }
  
  
  # ----------------------------------------------------------
  # 5-1) Correlation statistics
  # ----------------------------------------------------------
  results <- map_dfr(
    metabo_cols,
    function(metabo) {
      
      x <- suppressWarnings(
        as.numeric(data$Year)
      )
      
      y <- suppressWarnings(
        as.numeric(data[[metabo]])
      )
      
      valid_idx <- is.finite(x) &
        is.finite(y)
      
      n_valid <- sum(valid_idx)
      
      if (n_valid < 3) {
        return(NULL)
      }
      
      y_variance <- var(
        y[valid_idx],
        na.rm = TRUE
      )
      
      if (
        is.na(y_variance) ||
        !is.finite(y_variance) ||
        y_variance == 0
      ) {
        return(NULL)
      }
      
      x_variance <- var(
        x[valid_idx],
        na.rm = TRUE
      )
      
      if (
        is.na(x_variance) ||
        !is.finite(x_variance) ||
        x_variance == 0
      ) {
        return(NULL)
      }
      
      test <- tryCatch(
        cor.test(
          x[valid_idx],
          y[valid_idx],
          method = "pearson"
        ),
        error = function(e) NULL
      )
      
      if (is.null(test)) {
        return(NULL)
      }
      
      tibble(
        Metabolite = metabo,
        n = n_valid,
        r = as.numeric(test$estimate),
        p_value = test$p.value
      )
    }
  )
  
  if (nrow(results) == 0) {
    stop(
      paste0(
        "유효한 상관분석 결과가 없습니다.\n",
        "Year와 intensity 값, 결측치 및 분산을 확인하세요."
      )
    )
  }
  
  results <- results %>%
    mutate(
      FDR = p.adjust(
        p_value,
        method = "BH"
      ),
      Significant = !is.na(FDR) &
        FDR < fdr_cutoff,
      Direction = case_when(
        Significant & r > 0 ~ "Positive",
        Significant & r < 0 ~ "Negative",
        TRUE ~ "Not significant"
      )
    ) %>%
    left_join(
      name_map,
      by = "Metabolite"
    ) %>%
    arrange(
      FDR,
      p_value
    )
  
  write.csv(
    results,
    file = file.path(
      output_directory,
      "Figure3_linegraph_statistics.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  
  # ----------------------------------------------------------
  # 5-2) Select metabolites
  # ----------------------------------------------------------
  # Negative: FDR 기준 significant negative 상위 10개
  # Positive: significant positive 전체
  
  top_neg <- results %>%
    filter(
      Direction == "Negative"
    ) %>%
    arrange(
      FDR,
      desc(abs(r))
    ) %>%
    slice_head(
      n = neg_line_n
    ) %>%
    pull(Metabolite)
  
  top_pos <- results %>%
    filter(
      Direction == "Positive"
    ) %>%
    arrange(
      FDR,
      desc(abs(r))
    ) %>%
    pull(Metabolite)
  
  selected_list <- bind_rows(
    results %>%
      filter(
        Metabolite %in% top_neg
      ) %>%
      arrange(FDR) %>%
      mutate(
        Group = "Negative top metabolites"
      ),
    
    results %>%
      filter(
        Metabolite %in% top_pos
      ) %>%
      arrange(FDR) %>%
      mutate(
        Group = "Positive significant"
      )
  ) %>%
    select(
      Group,
      Metabolite,
      Display_Name,
      n,
      r,
      p_value,
      FDR
    )
  
  write.csv(
    selected_list,
    file = file.path(
      output_directory,
      "Figure3_linegraph_selected_metabolites.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  message(
    "Significant negative metabolites: ",
    sum(results$Direction == "Negative")
  )
  
  message(
    "Negative metabolites selected: ",
    length(top_neg)
  )
  
  message(
    "Positive metabolites selected: ",
    length(top_pos)
  )
  
  
  # ----------------------------------------------------------
  # 5-3) Prepare trajectory data
  # ----------------------------------------------------------
  # 개별 샘플 intensity를 대사체별 Z-score로 변환한 뒤
  # 연도별 평균 Z-score 계산
  
  df_long <- data %>%
    pivot_longer(
      cols = all_of(metabo_cols),
      names_to = "Metabolite",
      values_to = "Intensity"
    ) %>%
    mutate(
      Intensity = suppressWarnings(
        as.numeric(Intensity)
      )
    ) %>%
    filter(
      is.finite(Year),
      is.finite(Intensity)
    ) %>%
    group_by(Metabolite) %>%
    mutate(
      Metabolite_Mean = mean(
        Intensity,
        na.rm = TRUE
      ),
      Metabolite_SD = sd(
        Intensity,
        na.rm = TRUE
      ),
      Z_score = case_when(
        is.finite(Metabolite_SD) &
          Metabolite_SD > 0 ~
          (Intensity - Metabolite_Mean) /
          Metabolite_SD,
        
        TRUE ~ NA_real_
      )
    ) %>%
    ungroup() %>%
    left_join(
      name_map,
      by = "Metabolite"
    )
  
  df_year_metabolite <- df_long %>%
    filter(
      is.finite(Z_score)
    ) %>%
    group_by(
      Year,
      Metabolite,
      Display_Name
    ) %>%
    summarise(
      Mean_Z = mean(
        Z_score,
        na.rm = TRUE
      ),
      SD_Z = sd(
        Z_score,
        na.rm = TRUE
      ),
      Sample_N = sum(
        is.finite(Z_score)
      ),
      SE_Z = if_else(
        Sample_N > 1,
        SD_Z / sqrt(Sample_N),
        NA_real_
      ),
      .groups = "drop"
    )
  
  if (nrow(df_year_metabolite) == 0) {
    stop("연도별 Z-score 데이터를 만들 수 없습니다.")
  }
  
  
  # ----------------------------------------------------------
  # 5-4) Function for trajectory plot
  # ----------------------------------------------------------
  create_traj_plot <- function(
    target_metabolites,
    direction = c(
      "Negative",
      "Positive"
    ),
    line_color) {
    
    direction <- match.arg(direction)
    
    if (length(target_metabolites) == 0) {
      return(
        ggplot() +
          annotate(
            geom = "text",
            x = 0,
            y = 0,
            label = paste0(
              "No significant ",
              tolower(direction),
              " metabolites"
            ),
            size = 4
          ) +
          xlim(-1, 1) +
          ylim(-1, 1) +
          theme_void() +
          labs(
            title = paste0(
              direction,
              " metabolites"
            )
          )
      )
    }
    
    df_sub <- df_year_metabolite %>%
      filter(
        Metabolite %in%
          target_metabolites
      )
    
    if (nrow(df_sub) == 0) {
      stop(
        paste0(
          direction,
          " trajectory data가 없습니다."
        )
      )
    }
    
    # 선택된 여러 대사체의 평균 trajectory
    # ribbon은 선택 대사체 간 SD
    df_summary <- df_sub %>%
      group_by(Year) %>%
      summarise(
        mean_Z = mean(
          Mean_Z,
          na.rm = TRUE
        ),
        sd_Z = sd(
          Mean_Z,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      mutate(
        sd_Z = replace_na(
          sd_Z,
          0
        )
      )
    
    if (direction == "Negative") {
      
      anchor_year <- min(
        df_sub$Year,
        na.rm = TRUE
      )
      
      df_label <- df_sub %>%
        filter(
          Year == anchor_year
        )
      
      x_nudge <- -5
      h_just <- 1
      
      plot_title <- paste0(
        "Negative metabolites (top ",
        length(target_metabolites),
        ")"
      )
      
    } else {
      
      anchor_year <- max(
        df_sub$Year,
        na.rm = TRUE
      )
      
      df_label <- df_sub %>%
        filter(
          Year == anchor_year
        )
      
      x_nudge <- 5
      h_just <- 0
      
      plot_title <- paste0(
        "Positive metabolites (n = ",
        length(target_metabolites),
        ")"
      )
    }
    
    ggplot() +
      
      # 개별 대사체 trajectory
      geom_line(
        data = df_sub,
        aes(
          x = Year,
          y = Mean_Z,
          group = Metabolite
        ),
        color = col_grey,
        alpha = 0.65,
        linewidth = 0.55
      ) +
      
      geom_point(
        data = df_sub,
        aes(
          x = Year,
          y = Mean_Z,
          group = Metabolite
        ),
        color = col_grey,
        alpha = 0.65,
        size = 1.3
      ) +
      
      # 선택 대사체 간 평균 ± SD
      geom_ribbon(
        data = df_summary,
        aes(
          x = Year,
          ymin = mean_Z - sd_Z,
          ymax = mean_Z + sd_Z
        ),
        fill = line_color,
        alpha = 0.14
      ) +
      
      geom_line(
        data = df_summary,
        aes(
          x = Year,
          y = mean_Z
        ),
        color = line_color,
        linewidth = 1.3
      ) +
      
      geom_point(
        data = df_summary,
        aes(
          x = Year,
          y = mean_Z
        ),
        color = line_color,
        size = 2.5
      ) +
      
      # 대사체 이름
      geom_text_repel(
        data = df_label,
        aes(
          x = Year,
          y = Mean_Z,
          label = Display_Name
        ),
        size = 2.5,
        color = "grey30",
        nudge_x = x_nudge,
        direction = "y",
        hjust = h_just,
        segment.size = 0.25,
        segment.color = "grey65",
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123
      ) +
      
      scale_x_continuous(
        breaks = sort(
          unique(data$Year)
        ),
        expand = expansion(
          mult = c(
            0.38,
            0.38
          )
        )
      ) +
      
      theme_classic(
        base_size = 11
      ) +
      
      theme(
        plot.title = element_text(
          color = line_color,
          size = 12,
          face = "bold",
          hjust = 0.5
        ),
        axis.text = element_text(
          color = "black",
          size = 10
        ),
        axis.title = element_text(
          color = "black",
          size = 11,
          face = "bold"
        ),
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5
        ),
        plot.margin = margin(
          t = 8,
          r = 12,
          b = 8,
          l = 12
        )
      ) +
      
      labs(
        title = plot_title,
        x = "Year",
        y = "Standardized intensity (Z-score)"
      )
  }
  
  
  # ----------------------------------------------------------
  # 5-5) Build plots
  # ----------------------------------------------------------
  p_neg <- create_traj_plot(
    target_metabolites = top_neg,
    direction = "Negative",
    line_color = col_negative
  )
  
  p_pos <- create_traj_plot(
    target_metabolites = top_pos,
    direction = "Positive",
    line_color = col_positive
  )
  
  final_plot <- (
    p_neg | p_pos
  ) +
    plot_layout(
      widths = c(1, 1)
    )
  
  # RStudio Plots 창에 표시
  print(final_plot)
  
  
  # ----------------------------------------------------------
  # 5-6) Save
  # ----------------------------------------------------------
  output_base <- paste0(
    "Figure3_linegraph_neg",
    length(top_neg),
    "_pos",
    length(top_pos)
  )
  
  png_file <- file.path(
    output_directory,
    paste0(output_base, ".png")
  )
  
  pdf_file <- file.path(
    output_directory,
    paste0(output_base, ".pdf")
  )
  
  svg_file <- file.path(
    output_directory,
    paste0(output_base, ".svg")
  )
  
  ggsave(
    filename = png_file,
    plot = final_plot,
    width = 10,
    height = 5,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = pdf_file,
    plot = final_plot,
    width = 10,
    height = 5,
    units = "in",
    bg = "white",
    device = cairo_pdf,
    limitsize = FALSE
  )
  
  ggsave(
    filename = svg_file,
    plot = final_plot,
    width = 10,
    height = 5,
    units = "in",
    bg = "white",
    device = svglite,
    limitsize = FALSE
  )
  
  message("")
  message("Figure saved successfully:")
  message("PNG: ", png_file)
  message("PDF: ", pdf_file)
  message("SVG: ", svg_file)
  message("")
  message(
    "Statistics: ",
    file.path(
      output_directory,
      "Figure3_linegraph_statistics.csv"
    )
  )
  message(
    "Selected metabolites: ",
    file.path(
      output_directory,
      "Figure3_linegraph_selected_metabolites.csv"
    )
  )
  
  return(final_plot)
}


# ------------------------------------------------------------
# 6) Run
# ------------------------------------------------------------
fig3_linegraph <- generate_line_graph(
  data = df_t,
  name_map = metabolite_name_map,
  output_directory = out_dir,
  fdr_cutoff = 0.05,
  neg_line_n = 10
)


# ------------------------------------------------------------
# 7) Explicitly display figure
# ------------------------------------------------------------
print(fig3_linegraph)