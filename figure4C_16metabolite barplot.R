library(tidyverse)
library(readxl)
library(tidyr)

# 1. Load the dataset from Excel
df <- read_excel("figure3C_16metabolite bar.xlsx")

# 2. Filter for significant metabolites
df_sig <- df %>% 
  filter(p_value < 0.05) %>%
  mutate(
    Direction = ifelse(slope > 0, "Positive", "Negative")
  )

# 3. Custom order for y-axis (Component)
custom_comp_order <- c(
  "Glucose", "Fructose", "Sucrose", 
  "Tartaric Acid", "Malic Acid", "Citric Acid", 
  "Tannic Acid", 
  "Proline", "Arginine", "Glutamic Acid", "Alanine", 
  "Serine", "Threonine", "Valine", "Histidine", "Isoleucine"
)

df_sig$Component <- factor(df_sig$Component, levels = rev(custom_comp_order))

# 4. Custom order for metabolite classes
#    -> this order will be used on BOTH positive and negative sides
#    -> first class in the vector will be closest to 0
custom_cat_order <- c(
  "Sugars",
  "Sugar Alcohols",
  "Sugar Acids",
  "Organic Acids",
  "Amino Acids",
  "Amines",
  "Alcohols",
  "Polyphenols",
  "Vitamins",
  "Others"
)

df_sig$Category <- factor(df_sig$Category, levels = custom_cat_order)

# 5. Count metabolites by Component, Category, Direction
counts <- df_sig %>%
  group_by(Component, Category, Direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  
  # Fill missing combinations with 0 to keep stack order consistent
  complete(Component, Category, Direction, fill = list(Count = 0)) %>%
  
  mutate(
    Plot_Count = ifelse(Direction == "Negative", -Count, Count),
    Direction = factor(Direction, levels = c("Negative", "Positive"))
  )

# 6. Softer, more readable color palette
cat_colors <- c(
  "Sugars"         = "#D95F5F",
  "Sugar Alcohols" = "#E6A157",
  "Sugar Acids"    = "#F2C879",
  "Organic Acids"  = "#7AA874",
  "Amino Acids"    = "#6FA8DC",
  "Amines"         = "#9E8AC9",
  "Alcohols"       = "#5C88B0",
  "Polyphenols"    = "#B07AA1",
  "Vitamins"       = "#76B7B2",
  "Others" = "#E0E0E0"
)

# 7. Plot
max_count <- counts %>%
  group_by(Component, Direction) %>%
  summarise(Total = sum(abs(Plot_Count)), .groups = "drop") %>%
  summarise(Max_Total = max(Total, na.rm = TRUE)) %>%
  pull(Max_Total)

x_lim <- ceiling(max_count / 5) * 5 + 2

fig1 <- ggplot(counts, aes(x = Plot_Count, y = Component, fill = Category)) +
  
  annotate(
    "rect",
    xmin = -x_lim, xmax = 0,
    ymin = -Inf, ymax = Inf,
    fill = "#d7ecff",
    alpha = 0.45
  ) +
  annotate(
    "rect",
    xmin = 0, xmax = x_lim,
    ymin = -Inf, ymax = Inf,
    fill = "#ffe0ea",
    alpha = 0.45
  ) +
  
  geom_col(
    position = position_stack(reverse = TRUE),
    color = "white",
    linewidth = 0.25,
    width = 0.8
  ) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.9) +
  scale_fill_manual(
    values = cat_colors,
    breaks = custom_cat_order,
    drop = FALSE
  ) +
  scale_x_continuous(
    labels = abs,
    limits = c(-x_lim, x_lim)
  ) +
  labs(
    x = "Number of significant metabolites",
    y = NULL,
    fill = "Metabolite class"
  ) +
  theme_classic(base_size = 12, base_family="arial") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y.left = element_text(color = "black"),
    panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(reverse = FALSE))

print(fig1)

ggsave("Figure3C_class_order_fixed.png", fig1, width = 9, height = 6.5, dpi = 300)
ggsave("Figure3C_class_order_fixed.svg", fig1, width = 9, height = 6.5)
ggsave(
  "Figure3C_class_order_fixed.pdf",
  plot = fig1,
  width = 9,
  height = 6.5,
  device = cairo_pdf
)