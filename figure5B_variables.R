# =========================================================
# Figure 4D: Comparative categorized heatmap
# Version: 260612
# Input file: Figure 4D Raw data_260612.xlsx
# Sheets:
#   raw_data_past
#   raw_data_present
#   raw_data_future
# =========================================================

# ---------------------------------------------------------
# 0. CLEAR WORKSPACE
# ---------------------------------------------------------
rm(list = ls())

local({
  
  # ---------------------------------------------------------
  # 1. LOAD REQUIRED LIBRARIES
  # ---------------------------------------------------------
  packages <- c(
    "tidyverse",
    "readxl",
    "showtext",
    "sysfonts",
    "scales",
    "svglite"
  )
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  
  library(tidyverse)
  library(readxl)
  library(showtext)
  library(sysfonts)
  library(scales)
  library(svglite)
  
  # ---------------------------------------------------------
  # 2. VERSION / FILE PATH SETTING
  # ---------------------------------------------------------
  version_tag <- "v260612"
  
  file_path <- "Figure4D Raw data_260612.xlsx"
  out_dir <- getwd()
  
  output_prefix <- paste0("Figure4D_", version_tag, "_Comparative_Categorized_Heatmap")
  
  # ---------------------------------------------------------
  # 3. FONT SETTING
  # ---------------------------------------------------------
  font_add(
    family = "Arial",
    regular = "C:/Windows/Fonts/arial.ttf",
    bold = "C:/Windows/Fonts/arialbd.ttf"
  )
  
  showtext_auto()
  showtext_opts(dpi = 300)
  
  # =========================================================
  # 4. PROCESS PAST DATA
  # =========================================================
  df_past_v260612 <- read_excel(
    file_path,
    sheet = "raw_data_past",
    col_names = TRUE,
    .name_repair = "unique"
  ) %>%
    as.data.frame()
  
  colnames(df_past_v260612) <- make.unique(colnames(df_past_v260612))
  
  past_factors_v260612 <- c(
    "Temperature (C)",
    "Global Radiation (J/m2)",
    "Precipitation (mm)",
    "Humidity (%)",
    "Year"
  )
  
  # Past sheet 구조:
  # 1열 = Metabolite Name
  # 2~6열 = climate variables
  # 7열부터 = metabolites
  past_metabs_v260612 <- colnames(df_past_v260612)[7:ncol(df_past_v260612)]
  
  # 첫 번째 row에 compound group 정보가 들어 있음
  past_group_dict_v260612 <- data.frame(
    Metabolite = past_metabs_v260612,
    Compound_Group = as.character(df_past_v260612[1, past_metabs_v260612]),
    stringsAsFactors = FALSE
  )
  
  # group row 제거
  df_past_v260612 <- df_past_v260612[-1, ]
  
  # 숫자형 변환
  df_past_v260612 <- df_past_v260612 %>%
    mutate(
      across(
        all_of(c(past_factors_v260612, past_metabs_v260612)),
        ~ as.numeric(.x)
      )
    )
  
  # scaling
  df_past_scaled_v260612 <- df_past_v260612 %>%
    mutate(
      across(
        all_of(c(past_factors_v260612, past_metabs_v260612)),
        ~ as.numeric(scale(.x))
      )
    )
  
  # simple 1-on-1 regressions
  past_results_v260612 <- list()
  
  for (metab in past_metabs_v260612) {
    for (f in past_factors_v260612) {
      
      form <- as.formula(paste0("`", metab, "` ~ `", f, "`"))
      
      tryCatch({
        model <- lm(form, data = df_past_scaled_v260612)
        
        coefs <- summary(model)$coefficients[2, "Estimate"]
        pvals <- summary(model)$coefficients[2, "Pr(>|t|)"]
        
        past_results_v260612[[paste(metab, f, sep = "_")]] <- data.frame(
          Metabolite = metab,
          Factor = f,
          Slope = coefs,
          P_value = pvals,
          Era = "1. Past\nClimate variables",
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    }
  }
  
  past_df_v260612 <- bind_rows(past_results_v260612) %>%
    group_by(Factor) %>%
    mutate(FDR_value = p.adjust(P_value, method = "BH")) %>%
    ungroup()
  
  # =========================================================
  # 5. PROCESS PRESENT DATA
  # =========================================================
  df_present_v260612 <- read_excel(
    file_path,
    sheet = "raw_data_present",
    col_names = TRUE,
    .name_repair = "unique"
  ) %>%
    as.data.frame()
  
  colnames(df_present_v260612) <- make.unique(colnames(df_present_v260612))
  
  regions_v260612 <- c(
    "Australia",
    "Europe",
    "North America",
    "South America",
    "Global"
  )
  
  # Present sheet 구조:
  # 1열 = Metabolite name
  # 2열 = Country
  # 3열 = Continent
  # 4열 = HI_7mo
  # 5열부터 = metabolites
  present_metabs_v260612 <- colnames(df_present_v260612)[5:ncol(df_present_v260612)]
  
  # 첫 번째 row는 compound group 정보이므로 제거
  df_present_v260612 <- df_present_v260612[-1, ]
  
  df_present_v260612 <- df_present_v260612 %>%
    mutate(
      HI_7mo = as.numeric(HI_7mo),
      across(
        all_of(present_metabs_v260612),
        ~ as.numeric(.x)
      )
    )
  
  present_results_v260612 <- list()
  
  for (region in regions_v260612) {
    
    if (region == "Global") {
      df_sub <- df_present_v260612
    } else {
      df_sub <- df_present_v260612 %>%
        filter(Continent == region)
    }
    
    if (nrow(df_sub) > 2) {
      
      df_sub_scaled <- df_sub %>%
        mutate(
          HI_7mo = as.numeric(scale(HI_7mo)),
          across(
            all_of(present_metabs_v260612),
            ~ as.numeric(scale(.x))
          )
        )
      
      for (metab in present_metabs_v260612) {
        
        tryCatch({
          model <- lm(
            as.formula(paste0("`", metab, "` ~ HI_7mo")),
            data = df_sub_scaled
          )
          
          coefs <- summary(model)$coefficients[2, "Estimate"]
          pvals <- summary(model)$coefficients[2, "Pr(>|t|)"]
          
          present_results_v260612[[paste(region, metab, sep = "_")]] <- data.frame(
            Metabolite = metab,
            Factor = region,
            Slope = coefs,
            P_value = pvals,
            Era = "2. Present\nHI_7mo vs regions",
            stringsAsFactors = FALSE
          )
        }, error = function(e) NULL)
      }
    }
  }
  
  present_df_v260612 <- bind_rows(present_results_v260612) %>%
    group_by(Factor) %>%
    mutate(FDR_value = p.adjust(P_value, method = "BH")) %>%
    ungroup()
  
  # =========================================================
  # 6. PROCESS FUTURE DATA
  # =========================================================
  df_future_raw_v260612 <- read_excel(
    file_path,
    sheet = "raw_data_future",
    col_names = TRUE,
    .name_repair = "unique"
  ) %>%
    as.data.frame()
  
  df_future_v260612 <- df_future_raw_v260612 %>%
    mutate(
      Slope = as.numeric(r_value),
      Factor = Component,
      P_value = as.numeric(p_value),
      Era = "3. Future\nPredictive models"
    ) %>%
    select(
      Metabolite,
      Factor,
      Slope,
      P_value,
      Era
    )
  
  future_df_v260612 <- df_future_v260612 %>%
    group_by(Factor) %>%
    mutate(FDR_value = p.adjust(P_value, method = "BH")) %>%
    ungroup()
  
  # =========================================================
  # 7. COMBINE DATA
  # =========================================================
  combined_df_v260612 <- bind_rows(
    past_df_v260612,
    present_df_v260612,
    future_df_v260612
  ) %>%
    mutate(
      Significance = case_when(
        P_value < 0.001 ~ "***",
        P_value < 0.01  ~ "**",
        P_value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )
  
  # Past와 Present에 공통으로 존재하는 metabolite만 기준으로 사용
  shared_metabolites_v260612 <- intersect(
    past_metabs_v260612,
    present_metabs_v260612
  )
  
  sig_metabolites_v260612 <- combined_df_v260612 %>%
    filter(Metabolite %in% shared_metabolites_v260612) %>%
    group_by(Metabolite) %>%
    filter(any(Significance != "")) %>%
    pull(Metabolite) %>%
    unique()
  
  heatmap_data_v260612 <- combined_df_v260612 %>%
    filter(Metabolite %in% sig_metabolites_v260612) %>%
    left_join(past_group_dict_v260612, by = "Metabolite")
  
  # ---------------------------------------------------------
  # 8. FACTOR ORDER
  # ---------------------------------------------------------
  factor_order_v260612 <- c(
    "Temperature (C)",
    "Global Radiation (J/m2)",
    "Precipitation (mm)",
    "Humidity (%)",
    "Year",
    "Australia",
    "Europe",
    "North America",
    "South America",
    "Global",
    "Glucose",
    "Fructose",
    "Sucrose",
    "Tartaric Acid",
    "Malic Acid",
    "Citric Acid",
    "Tannin",
    "Proline",
    "Arginine",
    "Glutamic Acid",
    "Alanine",
    "Serine",
    "Threonine",
    "Valine",
    "Histidine",
    "Isoleucine"
  )
  
  heatmap_data_v260612 <- heatmap_data_v260612 %>%
    mutate(
      Factor = factor(Factor, levels = factor_order_v260612),
      Compound_Group = ifelse(
        is.na(Compound_Group) | Compound_Group == "",
        "Uncategorized",
        Compound_Group
      )
    ) %>%
    filter(!is.na(Factor))
  
  # ---------------------------------------------------------
  # 9. METABOLITE ORDER
  # ---------------------------------------------------------
  heatmap_data_v260612 <- heatmap_data_v260612 %>%
    group_by(Compound_Group, Metabolite) %>%
    mutate(mean_abs_slope = mean(abs(Slope), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Compound_Group, desc(mean_abs_slope))
  
  metabolite_order_v260612 <- heatmap_data_v260612 %>%
    distinct(Compound_Group, Metabolite, mean_abs_slope) %>%
    arrange(Compound_Group, mean_abs_slope) %>%
    pull(Metabolite)
  
  heatmap_data_v260612 <- heatmap_data_v260612 %>%
    mutate(
      Metabolite = factor(Metabolite, levels = unique(metabolite_order_v260612))
    )
  
  # =========================================================
  # 10. DRAW HEATMAP
  # =========================================================
  heatmap_plot_v260612 <- ggplot(
    heatmap_data_v260612,
    aes(
      x = Factor,
      y = Metabolite,
      fill = Slope
    )
  ) +
    geom_tile(
      color = "gray90",
      linewidth = 0.2
    ) +
    
    geom_text(
      aes(label = Significance),
      family = "Arial",
      color = "black",
      size = 2.6,
      fontface = "bold",
      vjust = 0.75
    ) +
    
    scale_fill_gradient2(
      low = "#4F8FC9",
      mid = "white",
      high = "#F26B6B",
      midpoint = 0,
      limits = c(-1.5, 1.5),
      breaks = c(-1.5, 0, 1.5),
      labels = c("-1.5", "0", "1.5"),
      oob = scales::squish,
      name = "Standardized\nbeta slope\nor r-value"
    ) +
    
    facet_grid(
      Compound_Group ~ Era,
      scales = "free",
      space = "free"
    ) +
    
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal(base_size = 9, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      
      axis.text.x = element_text(
        family = "Arial",
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 7.5,
        face = "bold",
        color = "gray20"
      ),
      
      axis.text.y = element_text(
        family = "Arial",
        size = 5.2,
        face = "bold",
        color = "gray20"
      ),
      
      strip.text.x = element_text(
        family = "Arial",
        size = 9,
        face = "bold",
        margin = margin(b = 5)
      ),
      
      strip.text.y = element_text(
        family = "Arial",
        size = 7,
        face = "bold",
        angle = 0,
        hjust = 0
      ),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.45
      ),
      
      panel.grid = element_blank(),
      
      legend.position = "right",
      legend.title = element_text(
        family = "Arial",
        size = 7.5,
        face = "bold"
      ),
      legend.text = element_text(
        family = "Arial",
        size = 7
      ),
      
      plot.margin = margin(5, 5, 5, 5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # =========================================================
  # 11. SAVE FIGURES
  # =========================================================
  num_rows_v260612 <- length(unique(heatmap_data_v260612$Metabolite))
  num_groups_v260612 <- length(unique(heatmap_data_v260612$Compound_Group))
  
  dynamic_height_v260612 <- max(
    10,
    num_rows_v260612 * 0.16 + num_groups_v260612 * 0.45
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(output_prefix, ".pdf")),
    plot = heatmap_plot_v260612,
    width = 18,
    height = dynamic_height_v260612,
    device = cairo_pdf,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(output_prefix, ".png")),
    plot = heatmap_plot_v260612,
    width = 18,
    height = dynamic_height_v260612,
    dpi = 300,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(output_prefix, ".svg")),
    plot = heatmap_plot_v260612,
    width = 18,
    height = dynamic_height_v260612,
    limitsize = FALSE
  )
  
  # =========================================================
  # 12. EXPORT WIDE DATA
  # =========================================================
  heatmap_slopes_wide_v260612 <- heatmap_data_v260612 %>%
    select(
      Metabolite,
      Compound_Group,
      Factor,
      Slope
    ) %>%
    pivot_wider(
      names_from = Factor,
      values_from = Slope
    )
  
  write.csv(
    heatmap_slopes_wide_v260612,
    file.path(out_dir, paste0("Figure4D_", version_tag, "_3_Eras_Heatmap_Slopes.csv")),
    row.names = FALSE
  )
  
  heatmap_asterisks_wide_v260612 <- heatmap_data_v260612 %>%
    select(
      Metabolite,
      Compound_Group,
      Factor,
      Significance
    ) %>%
    pivot_wider(
      names_from = Factor,
      values_from = Significance,
      values_fill = ""
    )
  
  write.csv(
    heatmap_asterisks_wide_v260612,
    file.path(out_dir, paste0("Figure4D_", version_tag, "_3_Eras_Heatmap_Asterisks.csv")),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------
  # 13. PRINT
  # ---------------------------------------------------------
  print(heatmap_plot_v260612)
  
  cat("Done! Figure 4D files saved in:", out_dir, "\n")
  cat("Output prefix:", output_prefix, "\n")
  cat("Number of metabolites:", num_rows_v260612, "\n")
  cat("Figure height:", dynamic_height_v260612, "\n")
  
})