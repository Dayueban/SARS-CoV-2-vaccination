# Description: this script was arranged for metabolomics data analysis
# under the guidelines of article "A Large-Scale, Multicenter Serum Metabolite Biomarker Identification
# Study for the Early Detection of Hepatocellular Carcinoma"
# the mock data was from the R package statTarget
rm(list = ls())
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("statTarget")
library(statTarget)
library(impute)
datapath <- system.file('extdata',package = 'statTarget')
file <- paste(datapath,'data_example.csv',sep = '/')
example1 <- read.csv(file)

#' correct for phenotype data for further analysis
pheno <- read.table("pheno.csv",header = T,sep = ',',row.names = 1)

a <- pheno[,3:ncol(pheno)]
sex <- pheno[,1]
batch <- as.factor(pheno[,2])
residual <- matrix(0,nrow(a),ncol(a),dimnames = list(rownames(a),colnames(a)))
for (i in 1:ncol(a)){
  fit <- lm(a[,i] ~ sex + batch, na.action = na.exclude)
  residual[,i] <- residuals(fit)
}
colnames(residual) <- NULL
write.table(residual,file = "DLY_sexbatch_correct.txt",quote = FALSE,row.names= FALSE,col.names = FALSE)



## the second example data
example2 <- read.csv(choose.files(),check.names = F,header = T,row.names = 1)
# to impute the value NA by KNN method in the package impute
group_label <- as.factor(example2[,1])
example2[,1] <- NULL
example3 <- as.data.frame(t(example2))
example4 <- impute.knn(as.matrix(example3),k=9)
example5 <- as.data.frame(t(example4$data))

## PCAj无监督分析用于评估整体代谢物在组间的变化
library("FactoMineR")
library("factoextra")
# 数据间进行标度化处理
Z_score <- function(x){
  x <- (x - mean(x))/sd(x)
}
example5 <- apply(example5, 1, Z_score)
example5 <- as.data.frame(t(example5))
data.pca <- PCA(example5,graph = FALSE)
Group <- group_label
data <- cbind(example5, Group)
colnames(data)[ncol(data)] <- "class"
# Scree plot to determine the number of principal components
tiff("meta_scree.tiff",width=2000,height=2000,res=300)
fviz_eig(data.pca, addlabels = TRUE, ylim = c(0, 100))
dev.off()
# individual plot
tiff("meta_PCA.tiff",width=2000,height=2000,res=300)
ind.p <- fviz_pca_ind(data.pca, geom = "point", col.ind = data$class)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "metabolomic sets",
              #caption = "Source: factoextra",
              xlab = "PC1(53.2%)", ylab = "PC2(7.3%)",
              legend.title = "Group", legend.position = "top",
              ggtheme = theme_gray(), palette = "npg"
)
dev.off()

## pls-da分析
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ropls")
library(ropls)
bb <- read.csv("MSMS/meta_all.csv", header = T, row.names = 1, check.names = F)
Group <- c(rep("VB",30), rep("SV", 29))
bb_oplsda <- opls(bb,Group,predI = 1,orthoI = NA)
# calculate the vip value
VipVn <- getVipVn(opls(bb,Group,predI = 1,orthoI = NA,plot = FALSE))
VipVn <- getVipVn(bb_oplsda)
# extract the variable name which the VIP value > 1
VipVn_great_than_1 <- names(VipVn[VipVn > 1])
# univariate statistics to extract variables differred between groups
pvaVn <- apply(bb,2,function(feaVn) wilcox.test(feaVn ~ Group,paired=FALSE)$p.value)
pvaVn_adj <- p.adjust(pvaVn,method = "BH",n=length(pvaVn))
# extract the adjusted p-value less than 0.05
pvaVn_adj_less_0.05 <- names(pvaVn_adj[pvaVn_adj < 0.05])

nc <- ncol(bb)
FC <- rep(0,nc)
log2_FC <- rep(0,nc)
for (i in 1:ncol(bb)) {
  FC[i] <- (mean(as.numeric(bb[1:30,i]))+0.001)/ (mean(as.numeric(bb[31:nrow(bb),i]))+0.001)
  log2_FC[i] <- log2((mean(as.numeric(bb[1:30,i]))+0.001)/ (mean(as.numeric(bb[31:nrow(bb),i]))+0.001))
}
final_bb <- cbind(VipVn, pvaVn, FC, log2_FC, pvaVn_adj)
write.csv(final_bb, "metabo_diff.csv")



# use intersect function to obtain the variables both VIP value great than 1
# and p value less than 0.05
final_variable <- intersect(VipVn_great_than_1,pvaVn_adj_less_0.05)

# 热图展示13个final_variable的丰度
heatmap_df <- example4$data
heatmap_df <- as.data.frame(t(heatmap_df))
heatmap_df <- heatmap_df[,final_variable]
heatmap_df_log10 <- log10(heatmap_df)
heatmap_df_log10 <- as.matrix(t(heatmap_df_log10))
library(pheatmap)
annotation_col <- data.frame(Group=factor(rep(c("N","T"),c(25,61))))
rownames(annotation_col) <- colnames(heatmap_df_log10)
ann_colors <- list(Group=c(N="#4DAF4A",T="#984EA3"))
#plot.new()
tiff("metabo_biomarker_pheatmap.tiff",width = 3000,height = 700,res = 300)
pheatmap(heatmap_df_log10,treeheight_col = 20,annotation_colors = ann_colors,cluster_cols = FALSE,
         cluster_rows = TRUE,show_colnames = FALSE,annotation_col=annotation_col,
         gaps_col = c(25,25),border_color = "NA",cellwidth = 4.7,cellheight = 10)
dev.off()


## 建立模型
## first transform the Group 'N' and 'T' into '0' and '1'
example6 <- example5[,final_variable]
outcome = Group
outcome <- sub("N",0,outcome) # 将分类变量转为0,1变量
outcome <- sub("T",1,outcome)
# merge the data without classification with outcome
data_logi <- cbind(example6,outcome)
# 按照一定比例选择，训练集70%样本，验证集30%的样本
index <- sample(x=2,size = nrow(data_logi),replace = TRUE,prob = c(0.7,0.3))
traindata <- data_logi[index==1,]
testdata <- data_logi[index==2,]
# 先对训练集数据构建logit模型
logit_lm <- glm(outcome ~.,family = binomial,data=traindata)
# summary(logit_lm)
logit_step <- step(logit_lm,direction = "backward")
# summary(logit_step)
# 模型的显著性检验，不止变量要显著，模型也要显著
anova(object = logit_step,test = "Chisq")
# 模型对样本外数据（验证集）的预测
prob <- predict(object = logit_step,newdata=testdata,type="response")
prediction <- ifelse(prob >=0.5, 1,0)
prediction <- factor(prediction,levels = c(1,0),ordered = TRUE)
f <- table(testdata$outcome,prediction)
f
agreement <- as.vector(prediction) == testdata$outcome
table(agreement)
prop.table(table(agreement))

# ROC curve and AUC value calculation
library(pROC)
roc(testdata$outcome,prob)
roc1 <- roc(testdata$outcome,
            prob,
            percent=TRUE,
            partial.auc.correct=TRUE,
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
tiff("example_variable_pROC.tiff",width=2000,height=2000,res=300)
roc1 <- roc(testdata$outcome,prob,
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            plot=TRUE, percent=roc1$percent,col="#F781BF",cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
sens.ci <- ci.se(roc1, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col="#FF69B4")
plot(sens.ci, type="bars")
plot(roc1,col="black",add=T)
legend("bottomright",cex=1.2,c(paste("AUC=",round(roc1$ci[2],2),"%"),
                               paste("95% CI:",round(roc1$ci[1],2),"%-",round(roc1$ci[3],2),"%")))
dev.off()

## 使用随机森林进行潜在biomarker的选择
library(randomForest)

# 首先确定下后面要用到的参数mtry和ntree
n <- length(names(traindata))
set.seed(100)
names(traindata)[1:ncol(traindata)-1] <- paste("M",1:13,sep = "_")
# 选择mtry参数
for(i in 1:(n-1)){
  mtry_fit <- randomForest(outcome ~ ., data = traindata, mtry = i)
  err <- mean(mtry_fit$err.rate)
  print(err)  # 当mtry=2的时候误差最小
}

# 选择ntree参数
set.seed(100)
ntree_fit <- randomForest(outcome ~., data = traindata,mtry=2,ntree=1000)
plot(ntree_fit) # 在400左右模型内误差基本稳定，默认值为500，可以选默认值

# 变量重要性选择
set.seed(100)
#' rfcv是随机森林交叉验证函数：Random Forest Cross Validation
result <- rfcv(traindata[,-ncol(traindata)],traindata$outcome,cv.fold = 10)
result$error.cv #' 查看错误率表，21时错误率最低，为最佳模型
#' 绘制验证结果
with(result,plot(n.var,error.cv,log="x",type = "o",lwd=2)) # 结果是选择13个变量误差最小

# 构建模型
set.seed(100)
rf <- randomForest(outcome ~ ., data = traindata, mtry=2,ntree=400,importance=TRUE)
rf

# 验证并预测
names(testdata)[1:ncol(testdata)-1] <- paste("M",1:13,sep = "_")
pred1 <- predict(rf, testdata)
# 查看准确性
test_matrix <- table(pred1, testdata$outcome)
sum(diag(test_matrix))/sum(test_matrix)
# 0.8148
