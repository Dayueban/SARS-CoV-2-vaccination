setwd("C:\\Users\\Maozhang He\\Desktop\\nc_hcm_code-main\\nc_hcm_code-main\\data")
rm(list = ls())
library("pheatmap")
library("dendextend")
library("RColorBrewer")
# scale before heatmap clustering
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
fec_conj <- read.csv("cG_pheatmap_fec_conjugates copy.csv")
fec_conj_num <- read.csv("cG_pheatmap_fec_conjugates_1 copy.csv")
type <- c(rep("glucuronide",13),rep("sulfate",20))
sample <- c(rep("Conv",12),rep("GF",12))
fecPhase <- fec_conj_num 
rownames(fecPhase) <- as.vector(fec_conj$X)
fecPhase_norm <- t(apply(fecPhase, 1, cal_z_score))
fecPhase_samplerow <- data.frame(sample)
fecPhase_typecol <- data.frame(type)
row.names(fecPhase_samplerow) <- colnames(fecPhase)
row.names(fecPhase_typecol) <- rownames(fecPhase)
ann_colors <- list(sample = c(Conv = "#0073C299",GF = "#EFC00099"), 
                   type = c(glucuronide = "#CCD1D1", sulfate = "#D2B4DE"))
fecPhase_pheatmap <- pheatmap(fecPhase_norm,
                              color = colorRampPalette(c("navy", "#FEF9E7", "firebrick3"))(500),
                              annotation_col = fecPhase_samplerow,
                              annotation_row = fecPhase_typecol,
                              cutree_cols = 2,
                              cutree_rows = 2,
                              border_color = "black",
                              fontsize_row = 8,
                              fontsize_col = 7,
                              annotation_colors = ann_colors
)
save_pheatmap_png <- function(x, filename, width=5000, height=3600, res = 600) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(fecPhase_pheatmap, "fecPhase_pheatmap.png")
