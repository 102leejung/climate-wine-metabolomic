## -----------------------------------------------------------------------------------------------------------
# ===================================================================
# FULL SCRIPT: METABOLITE DECOUPLING VISUALIZATION & EXPORT
# ===================================================================

# 1. Install and load necessary packages automatically
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("ggrepel")) install.packages("ggrepel", dependencies = TRUE)
library(ggplot2)
library(ggrepel)

# 2. Load the raw data
raw_data <- read.csv("raw_data_past.csv", stringsAsFactors = FALSE, header = FALSE)

# 3. Extract the metadata (Rows 3 to the end)
years <- as.numeric(raw_data[3:nrow(raw_data), 6])
temperatures <- as.numeric(raw_data[3:nrow(raw_data), 2])

# 4. Create an empty data frame for the calculated R² values
results_df <- data.frame(
  Metabolite = character(),
  Category = character(),
  Temp_R2 = numeric(),
  Year_R2 = numeric(),
  stringsAsFactors = FALSE
)

# 5. Loop through every metabolite column to calculate X and Y coordinates
for (i in 7:ncol(raw_data)) {
  
  intensities <- as.numeric(raw_data[3:nrow(raw_data), i])
  
  temp_cor <- cor(temperatures, intensities, use = "complete.obs")
  temp_r2 <- temp_cor^2
  
  year_cor <- cor(years, intensities, use = "complete.obs")
  year_r2 <- year_cor^2
  
  results_df <- rbind(results_df, data.frame(
    Metabolite = raw_data[1, i],
    Category = raw_data[2, i],
    Temp_R2 = temp_r2,
    Year_R2 = year_r2
  ))
}

# 6. Define the exact colors based on the provided legend
custom_colors <- c(
  "Amines" = "#7D3C5D",        
  "Amino Acids" = "#FFC000",   
  "Organic Acids" = "#92D050", 
  "Sugars" = "#F47C7C",        
  "Sugar Acids" = "#FF66CC",   
  "Vitamins" = "#5BC0EB",      
  "Others" = "#707B8C",        
  "Polyphenols" = "#9B59B6",   
  "Alcohol" = "#A9CCE3",       
  "Sugar Alcohols" = "#FADBD8",
  "Fatty Acids" = "#F5B041"    
)

# -------------------------------------------------------------------
# 7. CALCULATE EXACT R² THRESHOLD FOR p = 0.05
# -------------------------------------------------------------------
n_samples <- 12
p_value <- 0.05
df <- n_samples - 2

# Calculate the critical t-value for a two-tailed test
t_critical <- qt(p = p_value / 2, df = df, lower.tail = FALSE)

# Convert the t-value directly into the R² threshold (~0.33175)
exact_threshold <- (t_critical^2) / (df + t_critical^2)

# 8. Generate the plot
quadrant_plot <- ggplot(results_df, aes(x = Temp_R2, y = Year_R2, color = Category)) +
  
  geom_point(size = 3.5, alpha = 0.85) +
  
  # Use the dynamically calculated threshold for the quadrant lines
  geom_vline(xintercept = exact_threshold, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = exact_threshold, linetype = "dashed", color = "gray50") +
  
  geom_text_repel(aes(label = Metabolite), 
                  size = 3, 
                  show.legend = FALSE, 
                  max.overlaps = Inf,
                  box.padding = 0.4) +
  
  scale_color_manual(values = custom_colors) +
  
  labs(
    title = "Metabolite Decoupling: Thermal Summation vs. Bottle Aging",
    x = bquote("Temperature Correlation ("*R^2*")"),
    y = bquote("Vintage / Aging Correlation ("*R^2*")"),
    color = "Metabolite Category"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# 9. Save the plot as a high-resolution PNG file
ggsave(filename = "Metabolite_Quadrant_Map_p05.png", 
       plot = quadrant_plot, 
       width = 12,       
       height = 8,       
       dpi = 300,        
       bg = "white")     

# 10. Confirmation message in the console showing the exact math
print(paste("Plot saved! The R² threshold for p = 0.05 with n =", n_samples, "was calculated as:", round(exact_threshold, 4)))


## -----------------------------------------------------------------------------------------------------------
# ===================================================================
# FULL SCRIPT: METABOLITE DECOUPLING (PRESENT DATA)
# HI_7mo (Continuous) vs. Continent (Categorical)
# ===================================================================

if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("ggrepel")) install.packages("ggrepel", dependencies = TRUE)
library(ggplot2)
library(ggrepel)

# 1. Load the Present Data
raw_data <- read.csv("raw_data_present.csv", stringsAsFactors = FALSE, header = FALSE)

# 2. Extract Metadata (Assuming rows 3+ are data)
# Col 3 = Continent, Col 4 = HI_7mo
continents <- as.factor(raw_data[3:nrow(raw_data), 3])
hi_values <- as.numeric(raw_data[3:nrow(raw_data), 4])

# Number of samples (n) and number of continents (k)
n_samples <- length(hi_values)
k_continents <- length(levels(continents))
p_value <- 0.05

# -------------------------------------------------------------------
# 3. CALCULATE EXACT R² THRESHOLDS FOR p = 0.05
# -------------------------------------------------------------------
# X-Axis Threshold (Continuous: HI_7mo) - Uses t-distribution
df_x <- n_samples - 2
t_critical <- qt(p = p_value / 2, df = df_x, lower.tail = FALSE)
threshold_x <- (t_critical^2) / (df_x + t_critical^2)

# Y-Axis Threshold (Categorical: Continent) - Uses F-distribution (ANOVA)
df1_y <- k_continents - 1
df2_y <- n_samples - k_continents
f_critical <- qf(p = p_value, df1 = df1_y, df2 = df2_y, lower.tail = FALSE)
threshold_y <- (df1_y * f_critical) / ((df1_y * f_critical) + df2_y)

# -------------------------------------------------------------------
# 4. CALCULATE METABOLITE COORDINATES
# -------------------------------------------------------------------
results_df <- data.frame(
  Metabolite = character(),
  Category = character(),
  HI_R2 = numeric(),
  Continent_R2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through metabolites (Starting at Column 5)
for (i in 5:ncol(raw_data)) {
  
  intensities <- as.numeric(raw_data[3:nrow(raw_data), i])
  
  # X-axis: Linear correlation with HI_7mo
  hi_cor <- cor(hi_values, intensities, use = "complete.obs")
  hi_r2 <- hi_cor^2
  
  # Y-axis: Variance explained by Continent (ANOVA R²)
  cont_lm <- lm(intensities ~ continents)
  cont_r2 <- summary(cont_lm)$r.squared
  
  results_df <- rbind(results_df, data.frame(
    Metabolite = raw_data[1, i],
    Category = raw_data[2, i],
    HI_R2 = hi_r2,
    Continent_R2 = cont_r2
  ))
}

# -------------------------------------------------------------------
# 5. VISUALIZATION & EXPORT
# -------------------------------------------------------------------
custom_colors <- c(
  "Amines" = "#7D3C5D",        
  "Amino Acids" = "#FFC000",   
  "Organic Acids" = "#92D050", 
  "Sugars" = "#F47C7C",        
  "Sugar Acids" = "#FF66CC",   
  "Vitamins" = "#5BC0EB",      
  "Others" = "#707B8C",        
  "Polyphenols" = "#9B59B6",   
  "Alcohol" = "#A9CCE3",       
  "Sugar Alcohols" = "#FADBD8",
  "Fatty Acids" = "#F5B041"    
)

quadrant_plot <- ggplot(results_df, aes(x = HI_R2, y = Continent_R2, color = Category)) +
  
  geom_point(size = 3.5, alpha = 0.85) +
  
  # Draw the dynamically calculated, precise p=0.05 thresholds
  geom_vline(xintercept = threshold_x, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = threshold_y, linetype = "dashed", color = "gray50") +
  
  geom_text_repel(aes(label = Metabolite), 
                  size = 3, 
                  show.legend = FALSE, 
                  max.overlaps = Inf,
                  box.padding = 0.4) +
  
  scale_color_manual(values = custom_colors) +
  
  labs(
    title = "Metabolite Decoupling: Thermal Summation vs. Continental Origin",
    x = bquote("Huglin Index Correlation ("*R^2*")"),
    y = bquote("Continental Variance Explained ("*R^2*")"),
    color = "Metabolite Category"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Save the high-resolution output
ggsave(filename = "Metabolite_Quadrant_Map_Present_Data.png", 
       plot = quadrant_plot, 
       width = 12,       
       height = 8,       
       dpi = 300,        
       bg = "white")     

print(paste("Plot saved! X-axis threshold:", round(threshold_x, 4), "| Y-axis threshold:", round(threshold_y, 4)))


## -----------------------------------------------------------------------------------------------------------
# ===================================================================
# FULL SCRIPT: METABOLITE DECOUPLING VISUALIZATION & EXPORT (PAST DATA)
# Includes export to PNG, PDF, and SVG
# ===================================================================

# 1. Install and load necessary packages automatically
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("ggrepel")) install.packages("ggrepel", dependencies = TRUE)
if (!require("svglite")) install.packages("svglite", dependencies = TRUE) # For high-quality SVG export

library(ggplot2)
library(ggrepel)
library(svglite)

# 2. Load the raw data
raw_data <- read.csv("raw_data_past.csv", stringsAsFactors = FALSE, header = FALSE)

# 3. Extract the metadata (Rows 3 to the end)
years <- as.numeric(raw_data[3:nrow(raw_data), 6])
temperatures <- as.numeric(raw_data[3:nrow(raw_data), 2])

# 4. Create an empty data frame for the calculated R² values
results_df <- data.frame(
  Metabolite = character(),
  Category = character(),
  Temp_R2 = numeric(),
  Year_R2 = numeric(),
  stringsAsFactors = FALSE
)

# 5. Loop through every metabolite column to calculate X and Y coordinates
for (i in 7:ncol(raw_data)) {
  
  intensities <- as.numeric(raw_data[3:nrow(raw_data), i])
  
  temp_cor <- cor(temperatures, intensities, use = "complete.obs")
  temp_r2 <- temp_cor^2
  
  year_cor <- cor(years, intensities, use = "complete.obs")
  year_r2 <- year_cor^2
  
  results_df <- rbind(results_df, data.frame(
    Metabolite = raw_data[1, i],
    Category = raw_data[2, i],
    Temp_R2 = temp_r2,
    Year_R2 = year_r2
  ))
}

# 6. Define the exact colors based on the provided legend
custom_colors <- c(
  "Amines" = "#7D3C5D",        
  "Amino Acids" = "#FFC000",   
  "Organic Acids" = "#92D050", 
  "Sugars" = "#F47C7C",        
  "Sugar Acids" = "#FF66CC",   
  "Vitamins" = "#5BC0EB",      
  "Others" = "#707B8C",        
  "Polyphenols" = "#9B59B6",   
  "Alcohol" = "#A9CCE3",       
  "Sugar Alcohols" = "#FADBD8",
  "Fatty Acids" = "#F5B041"    
)

# -------------------------------------------------------------------
# 7. CALCULATE EXACT R² THRESHOLD FOR p = 0.05
# -------------------------------------------------------------------
n_samples <- 12
p_value <- 0.05
df <- n_samples - 2

# Calculate the critical t-value for a two-tailed test
t_critical <- qt(p = p_value / 2, df = df, lower.tail = FALSE)

# Convert the t-value directly into the R² threshold (~0.33175)
exact_threshold <- (t_critical^2) / (df + t_critical^2)

# 8. Generate the plot
quadrant_plot <- ggplot(results_df, aes(x = Temp_R2, y = Year_R2, color = Category)) +
  
  geom_point(size = 3.5, alpha = 0.85) +
  
  # Use the dynamically calculated threshold for the quadrant lines
  geom_vline(xintercept = exact_threshold, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = exact_threshold, linetype = "dashed", color = "gray50") +
  
  geom_text_repel(aes(label = Metabolite), 
                  size = 3, 
                  show.legend = FALSE, 
                  max.overlaps = Inf,
                  box.padding = 0.4) +
  
  scale_color_manual(values = custom_colors) +
  
  labs(
    title = "Metabolite Decoupling: Thermal Summation vs. Bottle Aging",
    x = bquote("Temperature Correlation ("*R^2*")"),
    y = bquote("Vintage / Aging Correlation ("*R^2*")"),
    color = "Metabolite Category"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# -------------------------------------------------------------------
# 9. SAVE THE PLOT IN PNG, PDF, AND SVG FORMATS
# -------------------------------------------------------------------
base_filename <- "Metabolite_Quadrant_Map_p05"

# Save PNG (300 DPI Raster image)
ggsave(filename = paste0(base_filename, ".png"), 
       plot = quadrant_plot, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Save PDF (Vector document)
ggsave(filename = paste0(base_filename, ".pdf"), 
       plot = quadrant_plot, 
       width = 12, height = 8, bg = "white")

# Save SVG (Scalable Vector Graphic)
ggsave(filename = paste0(base_filename, ".svg"), 
       plot = quadrant_plot, 
       width = 12, height = 8, device = "svg", bg = "white")

# 10. Confirmation message in the console
print(paste("Plots successfully exported to PNG, PDF, and SVG! R² threshold:", round(exact_threshold, 4)))


## -----------------------------------------------------------------------------------------------------------
# ===================================================================
# FULL SCRIPT: METABOLITE DECOUPLING (PRESENT DATA)
# HI_7mo (Continuous) vs. Continent (Categorical)
# Includes export to PNG, PDF, and SVG
# ===================================================================

if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("ggrepel")) install.packages("ggrepel", dependencies = TRUE)
if (!require("svglite")) install.packages("svglite", dependencies = TRUE) # For high-quality SVG export

library(ggplot2)
library(ggrepel)
library(svglite)

# 1. Load the Present Data
raw_data <- read.csv("raw_data_present.csv", stringsAsFactors = FALSE, header = FALSE)

# 2. Extract Metadata (Assuming rows 3+ are data)
# Col 3 = Continent, Col 4 = HI_7mo
continents <- as.factor(raw_data[3:nrow(raw_data), 3])
hi_values <- as.numeric(raw_data[3:nrow(raw_data), 4])

# Number of samples (n) and number of continents (k)
n_samples <- length(hi_values)
k_continents <- length(levels(continents))
p_value <- 0.05

# -------------------------------------------------------------------
# 3. CALCULATE EXACT R² THRESHOLDS FOR p = 0.05
# -------------------------------------------------------------------
# X-Axis Threshold (Continuous: HI_7mo) - Uses t-distribution
df_x <- n_samples - 2
t_critical <- qt(p = p_value / 2, df = df_x, lower.tail = FALSE)
threshold_x <- (t_critical^2) / (df_x + t_critical^2)

# Y-Axis Threshold (Categorical: Continent) - Uses F-distribution (ANOVA)
df1_y <- k_continents - 1
df2_y <- n_samples - k_continents
f_critical <- qf(p = p_value, df1 = df1_y, df2 = df2_y, lower.tail = FALSE)
threshold_y <- (df1_y * f_critical) / ((df1_y * f_critical) + df2_y)

# -------------------------------------------------------------------
# 4. CALCULATE METABOLITE COORDINATES
# -------------------------------------------------------------------
results_df <- data.frame(
  Metabolite = character(),
  Category = character(),
  HI_R2 = numeric(),
  Continent_R2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through metabolites (Starting at Column 5)
for (i in 5:ncol(raw_data)) {
  
  intensities <- as.numeric(raw_data[3:nrow(raw_data), i])
  
  # X-axis: Linear correlation with HI_7mo
  hi_cor <- cor(hi_values, intensities, use = "complete.obs")
  hi_r2 <- hi_cor^2
  
  # Y-axis: Variance explained by Continent (ANOVA R²)
  cont_lm <- lm(intensities ~ continents)
  cont_r2 <- summary(cont_lm)$r.squared
  
  results_df <- rbind(results_df, data.frame(
    Metabolite = raw_data[1, i],
    Category = raw_data[2, i],
    HI_R2 = hi_r2,
    Continent_R2 = cont_r2
  ))
}

# -------------------------------------------------------------------
# 5. VISUALIZATION & EXPORT
# -------------------------------------------------------------------
custom_colors <- c(
  "Amines" = "#7D3C5D",        
  "Amino Acids" = "#FFC000",   
  "Organic Acids" = "#92D050", 
  "Sugars" = "#F47C7C",        
  "Sugar Acids" = "#FF66CC",   
  "Vitamins" = "#5BC0EB",      
  "Others" = "#707B8C",        
  "Polyphenols" = "#9B59B6",   
  "Alcohol" = "#A9CCE3",       
  "Sugar Alcohols" = "#FADBD8",
  "Fatty Acids" = "#F5B041"    
)

quadrant_plot <- ggplot(results_df, aes(x = HI_R2, y = Continent_R2, color = Category)) +
  
  geom_point(size = 3.5, alpha = 0.85) +
  
  # Draw the dynamically calculated, precise p=0.05 thresholds
  geom_vline(xintercept = threshold_x, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = threshold_y, linetype = "dashed", color = "gray50") +
  
  geom_text_repel(aes(label = Metabolite), 
                  size = 3, 
                  show.legend = FALSE, 
                  max.overlaps = Inf,
                  box.padding = 0.4) +
  
  scale_color_manual(values = custom_colors) +
  
  labs(
    title = "Metabolite Decoupling: Thermal Summation vs. Continental Origin",
    x = bquote("Huglin Index Correlation ("*R^2*")"),
    y = bquote("Continental Variance Explained ("*R^2*")"),
    color = "Metabolite Category"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# -------------------------------------------------------------------
# SAVE THE PLOT IN PNG, PDF, AND SVG FORMATS
# -------------------------------------------------------------------
base_filename <- "Metabolite_Quadrant_Map_Present_Data"

# Save PNG (300 DPI Raster image)
ggsave(filename = paste0(base_filename, ".png"), 
       plot = quadrant_plot, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Save PDF (Vector document)
ggsave(filename = paste0(base_filename, ".pdf"), 
       plot = quadrant_plot, 
       width = 12, height = 8, bg = "white")

# Save SVG (Scalable Vector Graphic)
ggsave(filename = paste0(base_filename, ".svg"), 
       plot = quadrant_plot, 
       width = 12, height = 8, device = "svg", bg = "white")

print(paste("Plots successfully exported to PNG, PDF, and SVG! X-axis threshold:", round(threshold_x, 4), "| Y-axis threshold:", round(threshold_y, 4)))

