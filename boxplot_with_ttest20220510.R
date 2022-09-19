# 该代码来自于文章：https://doi.org/10.1038/s41467-021-26209-8中的Figure 3c部分
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(grid)
library(cowplot)
ind_fec_neu <- read.csv("MSMS/msms_Intensity3.csv", header = T, check.names = F, row.names = 1)

#' 需要将M015个体删除
ind_fec_neu <- ind_fec_neu[-which(rownames(ind_fec_neu) == "M015"), ]
#' 确定自己想要作图的几个代谢物在第几列
#' dim(ind_fec_neu)
#' which(names(ind_fec_neu) %in% c("L-Glutamic acid", "gamma-Aminobutyric acid", "Succinic acid semialdehyde", "Succinic acid"))
nc <- c(3, 28, 103, 134)
ind_fec_neu2 <- ind_fec_neu[,nc]
label_names <- names(ind_fec_neu2)[1:4]
names(ind_fec_neu2) <- paste(rep("M", 4), seq(1,4), sep = "_")
ind_fec_neu2$group <- c(rep("VB", 29), rep("SV", 29))
Start = 1
Stop = 4
gvec <- vector("list",length = length(Start:Stop))
#gvec <- vector("list",length = length(nc))
for(i in Start:Stop){
  batch_neu <- ggplot(data = ind_fec_neu2, 
                      aes_string(x="group",y=names(ind_fec_neu2)[i], fill="group"))+   # 需要循环作图的时候需要用到aes_string参数，而不用aes参数：https://mp.weixin.qq.com/s/2tBbYJrHSxJUoZ97jg93Aw
    scale_fill_manual(values = c("#90CBD3" ,"#E16663"))+
    scale_x_discrete(limit = c("VB", "SV")) +
    geom_jitter(size = 1.5)+
    geom_boxplot(size=1.5, alpha=.6)+
    ggtitle(label_names[i])+
    xlab("")+
    ylab("Relative intensity (log10)")+
    theme_classic()+
    theme(legend.position = "none")+
    geom_signif(comparisons = list(c("VB","SV")),
                test = "wilcox.test",  # Welch's t-test
                test.args = list(alternative = "two.sided",var.equal=FALSE, paired=FALSE),
                map_signif_level = TRUE,
                textsize = 5,
                vjust = 0.3)+
    scale_y_continuous(labels=function(y) format(y,scientific=TRUE))+
    #scale_y_continuous(trans = "log10") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10, face = "bold", color = "black",vjust = -0.3),
          axis.text.y = element_text(size = 10, color = "black", face = "bold"))
  #batch_neu <- batch_neu + scale_y_continuous(trans = "log10")
  ggsave(batch_neu, file = paste0("metabo_boxplot/","plot_",colnames(ind_fec_neu)[i],".tiff"),
         units="in", width=3.6, height=4, dpi=300, compression = 'lzw')
  gvec[[i-Start+1]] <- batch_neu
}

pp <- plot_grid(plotlist = gvec, nrow = 4) # R package cowplot
ggsave(pp, file = "metabo_boxplot/merge.tiff", height = 12, width = 3, dpi = 300)

