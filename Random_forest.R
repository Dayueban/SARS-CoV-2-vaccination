rm(list=ls())
library(ggplot2)
library(readxl)
library(randomForest)
library(rfPermute)
library(ggthemes)
library(tidyverse)
# ----------- #
#  fecal VIP  #
# ----------- #
# 注意：msms_Intensity2_t.csv这个表根据前面Figure 3C结果，需要删除到M015号个体
vip_fecCID <- read.csv("MSMS/msms_Intensity2_t.csv", header = T, 
                       row.names = 1, check.names = F)
#vip_rownames <- c("Conv_M1","Conv_M2","Conv_M3","Conv_M4","Conv_M5","Conv_M6",
#                 "Conv_F1","Conv_F2","Conv_F3","Conv_F4","Conv_F5","Conv_F6",
#                  "GF_F1","GF_F2","GF_F3","GF_F4","GF_F5","GF_F6",
#                  "GF_M1","GF_M2","GF_M3","GF_M4","GF_M5","GF_M6")
#rownames(vip_fecCID) <- vip_rownames
metabo_names <- names(vip_fecCID)
names(vip_fecCID) <- paste(rep("M", ncol(vip_fecCID)), seq(1,100), sep = "_")
head(vip_fecCID,3)
str(vip_fecCID)
library("randomForest","rfPermute")
set.seed(2020)
Group <- c(rep("VB",29),rep("SV",29))
vip_fecCID <- cbind(vip_fecCID, Group)

rf_out_fecCID <- randomForest(as.factor(Group) ~ ., data = vip_fecCID)
plot(rf_out_fecCID)
rf_out_fecCID
# Extracts variable importance (Mean Decrease in Gini Index)
# Sorts by variable importance and relevels factors to match ordering
var_imp_fec <- tibble(variable=setdiff(colnames(vip_fecCID), "Group"),
                      importance=as.vector(importance(rf_out_fecCID)))
var_imp_fec <- arrange(var_imp_fec, desc(importance))
var_imp_fec$variable <- factor(var_imp_fec$variable, levels=var_imp_fec$variable)
pseu_names_index <- match(var_imp_fec$variable, names(vip_fecCID)) # 找代谢物的位置，下一步返回去匹配具体的代谢物名称
var_imp_fec$metabo <- metabo_names[pseu_names_index]
write.csv(var_imp_fec,"var_imp_fec.csv")
# 将输出的'var_imp_fec.csv'文件重新整理成下面要读取的形式
#var_imp_fec_top50 <- read_excel("var_imp_serum_top50.xlsx")
var_imp_fec_top50 <- read.csv("var_imp_fec.csv", header = T)
str(var_imp_fec_top50)
# nlogp <- -log10(var_imp_fec_top50$pvalue)
# var_imp_fec_top50$nlogp <- nlogp
vip_fec_bar <- ggplot(data=var_imp_fec_top50,aes(x=reorder(metabo,importance),y=importance,fill=updown))+
  geom_bar(stat = "identity",
           color = "white",
           width = 0.6,
           size = 0)+
  scale_fill_manual(values=c("#90CBD3","#E16663"))+
  labs(title = "Metabolites importance",
       x="serum metabolite variables \n top 40 out of 100",
       y="Mean decrease in Gini index")+
  coord_flip()+
  theme_linedraw()+
  theme(axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_blank())+
  theme(legend.background = element_rect(
    fill = "white",
    color = "black",
    size = 0.5),
    legend.position = c(0.6,0.12),
    legend.text = element_text(face = "bold",size = 11),
    legend.title = element_text(face = "bold",size = 12))+
  labs(fill = "Serum metabolites")
vip_fec_bar
ggsave("vip_serum20220112.tiff", units="in", width=5, height=6, dpi=300, compression = 'lzw')







rm(list = ls())
library(randomForest)
library(ggplot2)
library(caret)
library(tidyverse)
library(e1071)
df <- read.csv("MSMS/msms_Intensity3.csv", header = T, row.names = 1, check.names = F)
#' dim(df)
metabo_names <- names(df)
names(df) <- paste(rep("M", ncol(df)), seq(1, ncol(df)), sep = "_")
#' 删除M015个体
df <- df[-which(rownames(df) == "M015"), ]

Group <- c(rep("VB",29), rep("SV", 29))
outcome = Group
outcome <- sub("VB",0,outcome) # 将分类变量转为0,1变量
outcome <- sub("SV",1,outcome)
# merge the data without classification with outcome
X <- cbind(df,outcome)
# 按照一定比例选择，训练集70%样本，验证集30%的样本
ind <- sample(x=2,size = nrow(X),replace = TRUE,prob = c(0.67,0.33))
ind.train <- X[ind == 1,]
ind.test <- X[ind == 2,]

set.seed(123)

# 1.在随机森林算法的函数randomForest()中有两个非常重要的参数，
# 而这两个参数又将影响模型的准确性
# 2.它们分别是mtry和ntree。一般对mtry的选择是逐一尝试，直到找
# 到比较理想的值，ntree的选择可通过图形大致判断模型内误差稳定时的值。

# 选取randomforest函数中的mtry节点值，一般可默认，mtry指定节点中用于二叉树的变量个数
# 默认情况下为数据集变量个数的二次方根（分类模型）或三分之一（预测模型）

n <- length(names(ind.train))
set.seed(100)
errRate <- NULL
for(i in 1:(n-1)){
  mtry_fit <- randomForest(as.factor(outcome) ~ ., data = ind.train, mtry = i)
  err <- mean(mtry_fit$err.rate)
  errRate[i] <- err
}
# 选择平均误差最小的mtry
m <- which.min(errRate)
#print(m)
# m为12

# 之后选择ntree值，ntree指定随机森林所包含的决策树数目，默认为500,
set.seed(123)
ntree_fit <- randomForest(as.factor(outcome) ~ ., data = ind.train, mtry = 26, ntree = 1000)
plot(ntree_fit) # ntree到600以后就基本不变

# 根据以上结果，默认情况下的mtry效果更好，所以以mtry=12,ntree=500为参数构建随机森林模型。
rf.train <- randomForest(as.factor(outcome) ~ ., data = ind.train, mtry = 85,
                         importance = TRUE,proximity=TRUE,ntree = 500)
print(rf.train)
# OOB estimate of  error rate: 0%
#summary(rf)
#' 交叉验证选择features
# set.seed(123)
#' rfcv是随机森林交叉验证函数：Random Forest Cross Validation
# result <- rfcv(X[,-ncol(X)],X$outcome,cv.fold = 10)
# result$error.cv #' 查看错误率表，21时错误率最低，为最佳模型
#' 绘制验证结果
# with(result,plot(n.var,error.cv,log="x",type = "o",lwd=2))
# 使用replicate进行多次交叉验证，可选
result <- replicate(5, rfcv(ind.train[,-ncol(ind.train)],as.factor(ind.train$outcome),cv.fold = 20), simplify=FALSE)
error.cv <- sapply(result, "[[", "error.cv")
error.cv <- cbind(rowMeans(error.cv),error.cv)
n.var = rownames(error.cv) %>% as.numeric()
error.cv = error.cv[,2:6]
colnames(error.cv) = paste('err',1:5,sep='.')
err.mean = apply(error.cv,1,mean)
allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
# number of features selected
optimal = 5

# 表格输出
write.table(allerr, file = "metabo_rfcv.txt", sep = "\t", quote = F, row.names = T, col.names = T)
# the pre-setted parameters used for ploting afterwards
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))

p = ggplot() + 
  geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.5), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black') + 
  geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype="dashed") + 
  #  geom_hline(yintercept = min(allerr$err.mean), colour='black', lwd=0.36, linetype="dashed") + 
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(1, 3, 5, 10, 20, 40, 60)) + # , max(allerr$num)
  labs(title=paste('Training set (n = ', dim(ind.train)[1],')', sep = ''), 
       x='Number of Metabolites', y='Cross-validation error rate') + 
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep="")) + 
  main_theme
ggsave(p, file = "vaccine_metabo_rfcv.pdf", width = 100, height = 60, unit = 'mm')
p

###--------------------------------------######
###----       将最佳的那6个metabo        ######

imp <- as.data.frame(rf.train$importance)
# 随机森林中的特征重要性表示在该特征上拆分的所有节点的基尼不纯度减少的总和
# 减少的越多，说明该特征在分类过程中起的作用越大。
imp <- imp[order(imp[,1],decreasing = T),]
#head(imp)
#' 输出重要性排序的表格
#write.table(imp,file = "importance_class_20180719.txt",quote = F,sep = '\t',row.names = T,col.names = T)
#' 简单可视化
varImpPlot(rf.train,main = "Top 5-feature OTU importance",n.var = 6,bg=par("bg"),
           color = par("fg"),gcolor=par("fg"),lcolor="gray")

#' 美化feature贡献度柱状图
#' 分析选择top5效果最好
imp = head(imp, n=5)
write.table(imp,file = "importance_metabo_selected20211227.txt",quote = F,sep = '\t',row.names = T,col.names = T)
#' 反向排序X轴，让柱状图从上往下画
imp = imp[order(imp[,3],decreasing = F),]
imp$name <- rownames(imp)
imp$name <- factor(imp$name,levels = imp$name)
OTU_attribution <- read.csv("taxonomy/Species/species_attribution.csv",header = T,
                            check.names = F)
rownames(OTU_attribution) <- OTU_attribution$ID
OTU_attribution <- OTU_attribution[rownames(imp),]
imp$phylum <- OTU_attribution$Phylum
imp2 <- imp
rownames(imp2) <- paste(rownames(imp2),OTU_attribution$Species,sep = " ")
imp2$name <- rownames(imp2)
imp2$name <- factor(imp2$name,levels = imp2$name)
#' load ggplot2 package
library(ggplot2)
p <- ggplot(imp2,aes(x=name,y=MeanDecreaseAccuracy,fill=phylum)) +
  geom_bar(stat = "identity") + coord_flip() 
p <- p+theme_bw()
p <- p+theme(legend.position = c(0.85,0.15))
p <- p + theme(plot.background=element_blank(),panel.background=element_blank())
p <- p + theme(axis.text.x=element_text(face="bold",color="black",size=14),
               axis.text.y=element_text(face = "bold",color="black",size = 14),
               axis.title.y=element_text(size = 16,face = "bold"),axis.title.x=element_text(size = 16))
p <- p + scale_fill_manual(values= "blue4")
p <- p + ylim(0,0.05)
p + theme(legend.background = element_blank())
p <- p + theme(
  legend.key.height=unit(0.6,"cm"),
  legend.key.width=unit(0.6,"cm"),
  legend.text=element_text(lineheight=0.6,face="bold",size=14),
  legend.title=element_text(size=16,face = "bold"))
p <- p+labs(x="Species rank")
ggsave(paste("taxonomy/Species/rf_imp_feature" , ".tiff", sep=""), p, width = 12, height = 4)

###---------------------------------###
###---     热图展示其丰度      -----###
# ' 筛选个feature 展示
sub_abu = X[,rownames(imp)]
names(sub_abu) <- c("Pseudomonas veronii", "Pseudomonas fragi", "Sphingomonas yabuuchiae")
#' transposition the data
sub_abu <- as.data.frame(t(sub_abu))
#' table(X$outcome)
#'  0   1 
#' 108 185
#sub_abu$class <- X$outcome
#sub_abu$class <- sub("0","Gilts",sub_abu$class)
#sub_abu$class <- sub("1","Entire_boars",sub_abu$class)
sub_abu_matrix <- as.matrix(sub_abu)
library(pheatmap)
annotation_col <- data.frame(Group=factor(rep(c("Mild","Health"),c(10,10))))
rownames(annotation_col) <- colnames(sub_abu_matrix)
ann_colors <- list(Group=c(Mild="#FB8072",Health="#80B1D3"))
#plot.new()
tiff("species_rf_pheatmap.tiff",width = 2600,height = 800,res = 300)
pheatmap(sub_abu_matrix,treeheight_col = 10,annotation_colors = ann_colors,cluster_cols = FALSE,
         cluster_rows = FALSE,show_colnames = FALSE,annotation_col=annotation_col,
         gaps_col = c(10,10),border_color = "black",cellwidth = 18,cellheight = 20)
dev.off()

#' 用选择好的3个OTU重新构建randomforest模型
ind.train2 <- cbind(ind.train[,rownames(imp)], ind.train$outcome)
names(ind.train2)[ncol(ind.train2)] <- "class"
rf.train2 <- randomForest(as.factor(class) ~ ., data = ind.train2, importance = TRUE,
                          proximity = TRUE, ntree = )
print(rf.train2)
# OOB estimate of  error rate: 2.5%
# 预测响应变量
ind.train2$predicted.response <- predict(rf.train2, ind.train2)

# confusionMatrix function from caret package can be used for 
# creating confusion matrix based on actual response 
# variable and predicted value(构建混淆矩阵建立真实响应变量和预测响应变量值)
confusionMatrix(data=ind.train2$predicted.response, reference=as.factor(ind.train2$class))

# 测试训练集
ind.test <- cbind(ind.test[,rownames(imp)], ind.test$outcome)
names(ind.test)[ncol(ind.test)] <- "class"
# 预测训练集的响应变量
ind.test$predicted.reponse <- predict(rf.train2, ind.test)
confusionMatrix(data=ind.test$predicted.reponse, reference=as.factor(ind.test$class))

#rf.test <- randomForest(class ~ ., data = ind.test,
#                        importance = TRUE,proximity=TRUE,ntree = 1000)
#rf.test.pre <- predict(rf.test,type="prob")
#p.train<-rf.test.pre[,2]
########ROC in train######
library(pROC)
predicted.reponse <- predict(rf.train2, ind.test, type = "prob")
roc(ind.test$class,predicted.reponse[,2])
###################
roc1 <- roc(ind.test$class,
            predicted.reponse[,2],
            partial.auc.correct=TRUE,
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
#####################
tiff("vaccine_metabo_pROC.tiff",width=2200,height=1800,res=300)
roc1 <- roc(ind.test$class, predicted.reponse[,2],
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,font=2, font.lab=2,
            plot=TRUE, percent=roc1$percent,col="#F781BF",cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
plot(roc1,col="black",add=T)
legend("bottomright",cex=1.4,c(paste("AUC=",round(roc1$ci[2],2)*100,"%"),
                               paste("95% CI:",round(roc1$ci[1],2)*100,"%-",round(roc1$ci[3],2),"%")))
dev.off()

















