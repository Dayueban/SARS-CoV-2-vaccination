rm(list = ls())
#a <- read.csv("MSMS/MSMSfinal_DEM.csv", header = T, check.names = F, row.names = 1 )
a <- read.csv("metabo_diff.csv", header = T, check.names = F, row.names = 1 )
names(a)[ncol(a)] <- "FDR"
names(a)[1] <- "VIP"
#a$log2FC <- as.numeric(a$log2FC)
#logFC_cutoff <- with(a, mean(abs(log2FC)) + 2*sd(abs(log2FC)))
logFC_cutoff <- 0.25
a$change <- as.factor(ifelse(a$FDR < 0.05 & 
                               abs(a$log2_FC) > logFC_cutoff,
                             ifelse(a$log2_FC > logFC_cutoff, 'UP',
                                    'DOWN'),'NOT'))
                      

this_tile <- paste0('Cutoff for log2FC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(a[a$change == 'UP',]),
                    '\nThe number of down gene is ',nrow(a[a$change == 'DOWN',]))
library(ggplot2)
library(ggpubr)
library(ggrepel)
b = data.frame(rownames(a),a)
colnames(b)[1] <- "ID"
#metabo_list <- as.factor(c("Acetate","Propionate","Butyrate"))
g <- ggplot(data=b, aes(x=log2_FC, y=-log10(FDR), color=change)) + 
  geom_point(alpha = 0.9,aes(size = VIP)) +
  geom_point() +
  theme_classic()+
  xlab("log2 fold change") + ylab("-log10 FDR") +
  #ggtitle(this_tile) + theme(plot.title = element_text(size = 15, hjust = 0.5)) +
  scale_color_manual(values = c("#E16663","grey","#90CBD3")) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.5) + # 根据前面设定的log2 FC阈值来设定线在横轴的位置
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5) +
  xlim(-8,9) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12))
ggsave(g, filename = "volcano_differ3.png", width = 6, height = 5, dpi = 300)

#b1 <- subset(b, grepl(paste(metabo_list, collapse = "|"), rownames(b)))
#g1 <- g + geom_text_repel(
#  data = b1,
#  #data$label <- metabo_list,
#  #aes(x=log2FC, y=-log10(Pvalue),label = label),
#  aes(label = rownames(b1)),
#  position = "identity",
#  size = 5,
#  box.padding = unit(0.5,'lines'),
#  point.padding = unit(0.8,'lines'),segment.color = "black",show.legend = FALSE
#)
#print(g)
#ggsave(g1,filename = "volcano_diff3.png", width = 6, height = 8)


ggplot(b,aes(log2_FC,-log10(FDR),fill=change))+
  geom_point(shape=21,aes(size=VIP,color=change))+
  scale_fill_manual(values=c('seagreen','gray','orange'))+
  scale_color_manual(values=c('seagreen','gray','orange'))+
  geom_vline(xintercept=c(-1,1),lty=2,col="gray30",lwd=0.6) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="gray30",lwd=0.6)+
  theme_bw(base_rect_size = 1)+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        plot.title = element_text(family = 'regular',hjust = 0.5),
        legend.position = c(0.5, 1),
        legend.justification = c(0.5, 1),
        legend.key.height = unit(0.5,'cm'),
        legend.background = element_rect(fill = NULL, colour = "black",size = 0.5))+
  xlim(-4,4)+
  guides(size=F,color=F)+
  ylab('-log10 (FDR)')+xlab('log2 (Fold Change)')


