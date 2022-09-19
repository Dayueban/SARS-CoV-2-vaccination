rm(list = ls())
library("pheatmap")
library("dendextend")
library("RColorBrewer")
# scale before heatmap clustering
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
fec_conj <- read.csv("MSMS/msms_Intensity2.csv", row.names = 1, header = T, check.names = F)
#fec_conj_num <- read.csv("MSMS/msms_Intensity2_1.csv")
#type <- c(rep("glucuronide",13),rep("sulfate",20))
type <- fec_conj$type
table(type)
# 运行一次结果后发现M015个体和SV组聚成一个cluster，因此建议删除
fec_conj <- fec_conj[,-16]
fec_conj_num <- fec_conj[,-1]
sample <- c(rep("VB",29),rep("SV",29))
fecPhase <- fec_conj_num 
#rownames(fecPhase) <- as.vector(fec_conj$X)
fecPhase_norm <- t(apply(fecPhase, 1, cal_z_score))
fecPhase_samplerow <- data.frame(sample)
fecPhase_typecol <- data.frame(type)
row.names(fecPhase_samplerow) <- colnames(fecPhase)
row.names(fecPhase_typecol) <- rownames(fecPhase)
#ann_colors <- list(sample = c(Conv = "#0073C299",GF = "#EFC00099"), 
#                  type = c(glucuronide = "#CCD1D1", sulfate = "#D2B4DE"))
ann_colors <- list(sample = c(VB = "#E16663",SV = "#90CBD3"),
                   type = c(Amino_acid = "#cc340c", Carbohydrate = "#3f60aa",
                            Cofactors_and_Vitamins = "#f18800", Lipid = "#e4ce00",
                            Nucleotide = "#9ec417", Unknow = "#13a983",
                            Xenobiotics = "#44c1f0"))
fecPhase_pheatmap <- pheatmap(fecPhase_norm,
                              color = colorRampPalette(c("navy", "#FEF9E7", "firebrick3"))(500),
                              annotation_col = fecPhase_samplerow,
                              annotation_row = fecPhase_typecol,
                              show_colnames = FALSE,
                              cutree_cols = 2,
                              cutree_rows = 2,
                              border_color = "black",
                              fontsize_row = 4,
                              fontsize_col = 6,
                              annotation_colors = ann_colors
)
save_pheatmap_png <- function(x, filename, width=5400, height=4500, res = 600) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(fecPhase_pheatmap, "fecPhase_pheatmap20220512.png")
