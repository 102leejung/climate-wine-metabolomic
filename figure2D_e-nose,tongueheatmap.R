library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)

# 함수: 변수 clustering 순서
get_col_order <- function(mat) {
  hc <- hclust(dist(t(mat)), method = "complete")
  colnames(mat)[hc$order]
}

# 1. 데이터 읽기
nose <- read.csv(
  "figure2D_e-nose.csv",
  row.names = 1,
  check.names = FALSE
)

tongue <- read.csv(
  "figure2D_e-tongue.csv",
  row.names = 1,
  check.names = FALSE
)

# 2. 전자코 변수 약어
colnames(nose) <- c(
  "DCE", "E-isoval", "3MMP", "E-hex",
  "2MP", "2MB", "Me-sal", "Dill"
)

# 3. z-score
nose_z <- as.matrix(scale(nose))
tongue_z <- as.matrix(scale(tongue))

# 4. 그룹 내부 clustering
nose_order <- get_col_order(nose_z)
tongue_order <- get_col_order(tongue_z)

# 5. 합치기
combined_z <- cbind(
  nose_z[, nose_order],
  tongue_z[, tongue_order]
)

# 6. x/y 반전
combined_z_t <- t(combined_z)

# 7. annotation
type <- c(
  rep("Electronic nose", length(nose_order)),
  rep("Electronic tongue", length(tongue_order))
)

ha_row <- rowAnnotation(
  Type = type,
  col = list(
    Type = c(
      "Electronic nose" = "#DD8452",
      "Electronic tongue" = "#4C72B0"
    )
  ),
  show_annotation_name = FALSE,
  width = unit(4, "mm")
)

# 8. 색상
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#313695", "white", "#A50026")
)

# 9. Heatmap 객체
ht <- Heatmap(
  combined_z_t,

  name = "Z",
  col = col_fun,

  left_annotation = ha_row,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  show_row_dend = FALSE,
  show_column_dend = FALSE,

  row_names_side = "left",

  # 연도 글씨와 변수 글씨 크기 동일
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),

  column_names_rot = 90,
  column_labels = colnames(combined_z_t),

  rect_gp = gpar(col = "white", lwd = 1),

  # 셀 정사각형 느낌 유지
  width  = unit(ncol(combined_z_t) * 5, "mm"),
  height = unit(nrow(combined_z_t) * 5, "mm"),

  use_raster = FALSE
)

# 10. 화면 출력
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# 11. PDF 저장
pdf("Figure1D_heatmap_transposed.pdf", width = 6, height = 5, useDingbats = FALSE)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

# 12. SVG 저장
svglite("Figure1D_heatmap_transposed.svg", width = 6, height = 5)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()