# 该代码来自于文章：https://doi.org/10.1038/s41467-021-26209-8中的Figure 3c部分
library(ggplot2)
library(ggpubr)
library(ggsignif)
ind_fec_neu <- read.csv("bp_feces_neurotransmitters.csv")
Start = 5
Stop = 23
gvec <- vector("list",length = length(Start:Stop))
for(i in Start:Stop){
  batch_neu <- ggplot(data = ind_fec_neu, 
                      aes_string(x="group",y=names(ind_fec_neu)[i], fill="group"))+
    scale_fill_manual(values = c("#0073C299","#EFC00099"))+
    geom_jitter(size = 2.7)+
    geom_boxplot(size=1.5, alpha=.6)+
    xlab("")+
    ylab("Peak intensity")+
    theme_classic(
      base_size = 30
    )+
    theme(legend.position = "none")+
    geom_signif(comparisons = list(c("CONV-R","GF")),
                test = "t.test",  # Welch's t-test
                test.args = list(alternative = "two.sided",var.equal=FALSE, paired=FALSE),
                map_signif_level = TRUE,
                textsize = 10,
                vjust = 0.67)+
    scale_y_continuous(labels=function(y) format(y,scientific=TRUE))+
    theme(axis.title.y = element_text(size = 22),
          axis.text.x = element_text(size = 18, face = "bold", color = "black",vjust = -0.3),
          axis.text.y = element_text(size = 16, color = "black", face = "bold"))
  ggsave(batch_neu, file = paste0("plot_",colnames(ind_fec_neu)[i],".tiff"),
         units="in", width=3.6, height=5.5, dpi=300, compression = 'lzw')
  gvec[[i-Start+1]] <- batch_neu
}