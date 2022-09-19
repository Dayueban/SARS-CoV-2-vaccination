# 该代码来自于文章：https://doi.org/10.1038/s41467-021-26209-8中的Figure 3c部分
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(grid)
library(cowplot)
ind_fec_neu <- read.csv("blood_routine_test.csv")
Start = 2
Stop = 19
gvec <- vector("list",length = length(Start:Stop))
for(i in Start:Stop){
  batch_neu <- ggplot(data = ind_fec_neu, 
                      aes_string(x="group",y=names(ind_fec_neu)[i], fill="group"))+
    scale_fill_manual(values = c("#90CBD3" ,"#E16663"))+
    scale_x_discrete(limit = c("VB", "SV")) +
    geom_jitter(size = 1.5)+
    geom_boxplot(size=1.5, alpha=.6)+
    xlab("")+
    ylab(names(ind_fec_neu)[i])+
    theme_classic(
      base_size = 10
    )+
    theme(legend.position = "none")+
    geom_signif(comparisons = list(c("VB","SV")),
                test = "wilcox.test",  # Welch's t-test
                test.args = list(alternative = "two.sided",var.equal=FALSE, paired=TRUE),
                map_signif_level = TRUE,
                textsize = 3,
                vjust = 0.3)+
    scale_y_continuous(labels=function(y) format(y,scientific=FALSE))+
    theme(axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10, face = "bold", color = "black",vjust = -0.3),
          axis.text.y = element_text(size = 10, color = "black", face = "bold"))
  ggsave(batch_neu, file = paste0("RBI_plot/","plot_",colnames(ind_fec_neu)[i],".tiff"),
         units="in", width=3.6, height=4, dpi=300, compression = 'lzw')
  gvec[[i-Start+1]] <- batch_neu
}
pp <- plot_grid(plotlist = gvec, nrow = 5, labels = LETTERS[seq(from=1, to=18)],
                label_size = 10, hjust = -0.8) # R package cowplot
ggsave(pp, file = "RBI_plot/merge.tiff", height = 10, width = 6, dpi = 300)

