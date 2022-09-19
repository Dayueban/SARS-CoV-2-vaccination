rm(list = ls())
library(plyr)
library(psych)
#install.packages("ggcorrplot")
library(ggcorrplot)
library(pheatmap)
library(reshape2)
library(tidyverse)
library(circlize)
df1 <- read.csv("cytokine_Ig_BRT.csv",header = T, row.names = 1, check.names = F)
df2 <- read.csv("MSMS/msms_Intensity2_t.csv", header = T, row.names = 1, check.names = F)
df2_1 <- df2[, which(names(df2) %in% c("gamma-Aminobutyric acid", "Indole", "L-Glutamic acid",
                                       "Succinic acid", "Succinic acid semialdehyde"))]
df1 <- df1[rownames(df2), ] # 确保两个表格行名是一样的
# Spearman correlation analysis
corr_df <- corr.test(df1, df2_1, method = "spearman", adjust = "BH",alpha = .05)
corr_df_cor <- corr_df$r
corr_df_p <- corr_df$p

# Reset rownames
corr_df_cor <- data.frame(row=rownames(corr_df_cor),corr_df_cor,check.names = F) # create a column called "row" 
rownames(corr_df_cor) <- NULL
corr_df_p <- data.frame(row=rownames(corr_df_p),corr_df_p,check.names = F) # create a column called "row" 
rownames(corr_df_p) <- NULL
# Melt
nbar.m <- melt(corr_df_cor)
nbap.m <- melt(corr_df_p)
# Classify (you can classify differently for nbar and for nbap also)         
nbar.m$value2<-cut(nbar.m$value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE, label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # the label for the legend
nbap.m$value2<-cut(nbap.m$value,breaks=c(-Inf, 0.001, 0.01, 0.05),label=c("***", "** ", "*  ")) 
nbar.m<-cbind.data.frame(nbar.m,nbap.m$value,nbap.m$value2) # adding the p value and its cut to the first dataset of R coefficients
names(nbar.m)[5]<-paste("valuep") # change the column names of the dataframe 
names(nbar.m)[6]<-paste("signif.")
nbar.m <- na.omit(nbar.m)

pa <- ggplot(nbar.m, aes(row, variable)) +
  geom_tile(aes(fill=value),colour="white") +
  #scale_fill_brewer(palette = "RdYlGn",name="Correlation")# RColorBrewer package
  scale_fill_gradient2(low="blue", high="red", guide="colorbar",name="correlation") +
  theme_classic() +
  theme(axis.text.x=element_text(face="bold",angle=45,color="black",vjust = 0.95,hjust = 0.95,size=13),
        axis.text.y=element_text(face = "bold",size = 12, color = "black"),
        axis.title.y=element_text(size = 20),axis.title.x=element_text(size = 20))+
  labs(title="Cecum_negative") + 
  xlab("Serum immune factors") +
  ylab("Serum metabolites") +
  theme(plot.title = element_text(size = 15,color = "black")) +
  theme(axis.line.x=element_line(colour = "black"),axis.line.y=element_line(colour = "black"),
        legend.key.height=unit(0.8,"cm"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=12),
        legend.title=element_text(size=15))
# Adding the significance level stars using geom_text 
pp<- pa +
  geom_text(aes(label=signif.),size=6,na.rm=TRUE)
ggsave("metabo_immune_assoc.tiff",width = 10,height = 6)


#' part two: chord diagram plot
# 添加一列数值1，分别计算相同代谢物的个数，用于后续排序
df3 <- nbar.m
write.csv(df3, file = "metabo_assoc_clinical_param.csv")
# remove column which the |correlation coefficients| less than 0.3
#df4 <- df3[-which(abs(df3$value) < 0.3) ,]
df4 <- df3[,1:3]
names(df4) <- c("immune", "metabo", "Coeff")
df4$value <- rep(1, nrow(df4))
# 根据列名taxo将相同细菌凑一起，然后根据value计算个数

data_total <- df4 %>%
  group_by(immune) %>% 
  transmute(Total=sum(value))
# 根据列名metabo将相同代谢物凑一起，然后根据value计算个数
data_total2 <- df4 %>%
  group_by(metabo) %>% 
  transmute(Total2=sum(value))
# 将data_total中的Total列合并到df4数据框中
df4 <- cbind(df4, data_total$Total, data_total2$Total2)
names(df4)[c(ncol(df4)-1, ncol(df4))] <- c("Total", "Total2")

# change the first column by paste metabo names with numbers in column "Total"
df4$taxa <- paste(df4$taxa, "\t(",df4$Total, ")",sep = "")
df4$metabo <- paste(df4$metabo, "\t(", df4$Total2,")", sep = "")
# sort df4 based on the numbers of links
df4 <- df4[order(df4$Total2, decreasing = T),]
# remove column 'value' and 'Total'
df4$value <- NULL
df4$Total <- NULL
df4$Total2 <- NULL
## set colours for segments
#df4 <- df4[-which(abs(df4$Corr_score) < 5), ] # only |Corr_score| > 3 were used for exhibition
circos.clear()
# 第一种图是将links根据正负相关显示的
#pdf(file = paste(temp_nm2,"_otus_metabo_circlize.pdf",sep = ""), height = 10, width = 10)
tiff(file = "immune_metabo_circlize20220511.tiff", height = 2400, width = 2400,
     res = 300)
#circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
#Put horizontally or vertically symmetric
# gap.after是为了将Otu和代谢物之间分一个gap
circos.par(start.degree = 0, gap.degree = 1, points.overflow.warning = FALSE,
           gap.after = c(rep(1.5, length(unique(df4$immune))-1), 5, 
                         rep(1.5, length(unique(df4$metabo))-1), 5),
           circle.margin = c(0.1,0.7,0.1,0.5)) #controls the margins on the left, right, bottom and top sides of the circle
#par(mar = rep(-2, 4))
#par(mar = c(-8,-8,0,0))

# 给代谢物上色，而Otu全部给予灰色
grid.col <- setNames(c(topo.colors(length(unique(df4$immune))),rep("#BEBEBE", length(unique(df4$metabo)))), 
                     c(unique(df4$immune), unique(df4$metabo)))
# now, plot the image with rotated labels
chordDiagram(df4,  
             preAllocateTracks = 1, 
             annotationTrack = "grid",
             annotationTrackHeight = c(0.03, 0.1),
             grid.col = grid.col,
             col = ifelse(df4$Coeff > 0, "#F76F72", "#A8D1DF"), # link为正负相关可以用两种颜色表示
             link.sort = TRUE, link.decreasing = FALSE)
#directional = 1, 
#direction.type = c("diffHeight", "arrows")
#link.arr.type = "big.arrow")

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.4), cex = 1.2)
}, bg.border = NA)
title(paste("Metabolites assocaited with" ,"clinical param", sep = " "))
legend("topleft", pch = 15, bty = "n",col = c("#F76F72", "#A8D1DF"), title = "Connection type",
       legend = c("Positive","Negative"), cex = 1.2)
dev.off()

# 图形后期通过PS修正一下


# scatter plot between γ-aminobutyric acid and IgG
df3 <- cbind(df1$IgM, df2$Taurine)
df3 <- as.data.frame(df3)
names(df3) <- c("IgM", "Taurine")
rownames(df3) <- rownames(df2)
df3_1 <- scale(df3, center = TRUE, scale = TRUE)
df3_1 <- as.data.frame(df3_1)
library(ggpubr)
ggscatter(df3_1, x = "IgM", y = "Taurine", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", cor.coef.size = 6,
          add.params = list(color = "#F76F72", fill = "lightgray"),
          xlab = "IgM", ylab = "Taurine")
ggsave("IgM_correlated_Taurine.tiff", height = 6, width = 6)
