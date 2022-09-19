## 建立模型
## first transform the Group 'N' and 'T' into '0' and '1'
rm(list = ls())
library(caret)
library(car)
df <- read.csv("MSMS/msms_Intensity2_tt.csv", header = T, row.names = 1, check.names = F)

set.seed(123)
#Z_score <- function(x){
#  x <- (x - mean(x))/sd(x)
#}
#df_z <- apply(df, 2, Z_score)
#df_z <- as.data.frame(df_z)
df_z <- log10(df)

Group <- c(rep("VB",29), rep("SV", 29))
outcome = Group
outcome <- sub("VB",0,outcome) # 将分类变量转为0,1变量
outcome <- sub("SV",1,outcome)
# merge the data without classification with outcome
data_logi <- cbind(df_z,outcome)
#data_logi <- cbind(df,outcome)
# 按照一定比例选择，训练集70%样本，验证集30%的样本
index <- sample(x=2,size = nrow(data_logi),replace = TRUE,prob = c(0.7,0.3))
traindata <- data_logi[index==1,]
testdata <- data_logi[index==2,]
# 先对训练集数据构建logit模型
# logit_lm <- glm(as.factor(outcome) ~ M1+M2+M3+M4+M5+M6+M7, 
#                family = binomial(), data = traindata)
#logit_lm <- glm(as.factor(outcome) ~ M1+M2+M3+M4, 
#                family =  binomial(link = "logit"), data = data_logi)
logit_glm <- glm(as.factor(outcome) ~ ., 
                family =  binomial(link = "logit"), data = traindata)
#logit_lm <- lm(outcome ~ ., data = data_logi)
#vif(logit_glm)
# summary(logit_glm)
logit_step <- step(logit_glm,direction = "backward")
# summary(logit_step)
# 模型的显著性检验，不止变量要显著，模型也要显著
anova(object = logit_step,test = "Chisq")

# Step:  AIC=6
# as.factor(outcome) ~ M3 + M11
metabo_index <- which(names(df_z) %in% c("gamma-Aminobutyric acid", "Indole"))
df_z1 <- df_z[, metabo_index]
df_z1 <- cbind(df_z1, outcome)
df_z1$outcome <- as.factor(df_z1$outcome)
#traindata1 <- df_z1[index==1,]
#testdata1 <- df_z1[index==2,]
#model_1 <- glm(as.factor(outcome) ~., data = traindata1, family = binomial)
#model_pre <- predict(model_1, type = "response", newdata = testdata1)

#merge_df <- cbind(as.numeric(testdata1$outcome), as.numeric(model_pre))

#' 合并四种代谢物进行回归模型预测
folds <- createFolds(y = df_z1[, ncol(df_z1)], k = 58)
max=0
min=0
fc<-as.numeric()
mod_pre<-as.numeric()
#' metabo_diff_1
for(i in 1:58){
  fold_test <- df_z1[folds[[i]],] # 其中一份拿来验证
  fold_train <- df_z1[-folds[[i]],] # 其余4份作为训练
  model <- glm(as.factor(outcome) ~ ., data=fold_train, family = binomial)
  model_pre <- predict(model, type = "response", newdata=fold_test)
  fc <- append(fc, fold_test$outcome)
  mod_pre <- append(mod_pre, as.numeric(model_pre))
}
merge_df<-cbind(fc,as.numeric(mod_pre))



#' 对结果进行ROC曲线的绘制
library(pROC)
roc_s <- roc(merge_df[,1],merge_df[,2], auc=TRUE, ci=TRUE)
mycoords <- coords(roc_s, "best", "threshold", transpose=TRUE, best.method="youden")
my.threshold = mycoords[[1]]
tiff("ROC-merge-metabolite20220907.tiff", width = 1600, height = 1500, res = 300)
mycol <- c("blue","darkmagenta","darkred","goldenrod2")
x<-plot.roc(merge_df[,1],merge_df[,2],
            smooth=F,
            lwd=2,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=T,
            main="",
            print.thres = my.threshold,
            print.thres.pch=20,
            print.thres.col="blue",
            print.thres.cex=1.5,
            identity.col="#57A3D7",
            identity.lwd=1.5,
            xlab="False positive rate",
            ylab="True positive rate",
            col=mycol[3])
legend.name <- c(paste("Combined:","AUC=","0.96",sep=" "),
                 paste("CI:", "0.89-1.0", sep = " "),
                 paste("Sensitivity:", "96.55%", sep = " "),
                 paste("Speciticity:", "96.55%", sep = " "))
legend("bottomright", 
       legend=legend.name,
       col = mycol[1:4],
       lwd = 2,
       bty="n")
dev.off()


#每个代谢物分别计算AUC值并作图
# GAMA
single_one1 <- df_z1[,c(2,3)]
folds <- createFolds(y = single_one1[,ncol(single_one1)], k=58)
max=0
min=0
fc<-as.numeric()
mod_pre<-as.numeric()
for(i in 1:58){
  fold_test<-single_one1[folds[[i]],] # 其中一份拿来验证
  fold_train<-single_one1[-folds[[i]],] # 其余4份作为训练
  model<-glm(as.factor(outcome) ~ ., data=fold_train, family = binomial)
  model_pre<-predict(model, type = "response", newdata=fold_test)
  fc <- append(fc, fold_test$outcome)
  mod_pre <- append(mod_pre, as.numeric(model_pre))
}
merge_df<-cbind(fc,as.numeric(mod_pre))

#' 对结果进行ROC曲线的绘制
library(pROC)
roc_s <- roc(merge_df[,1],merge_df[,2], auc=TRUE, ci=TRUE)
mycoords <- coords(roc_s, "best", "threshold", transpose=TRUE, best.method="youden")
my.threshold = mycoords[[1]]
tiff("ROC-Indole-metabolite20220907.tiff", width = 1600, height = 1500, res = 300)
mycol <- c("blue","darkmagenta","darkred","goldenrod2")
x<-plot.roc(merge_df[,1],merge_df[,2],
            smooth=F,
            lwd=2,
            ylim=c(0,1),
            xlim=c(1,0),
            legacy.axes=T,
            main="",
            print.thres = my.threshold,
            print.thres.pch=20,
            print.thres.col="blue",
            print.thres.cex=1.5,
            identity.col="#57A3D7",
            identity.lwd=1.5,
            xlab="False positive rate",
            ylab="True positive rate",
            col=mycol[3])
legend.name <- c(paste("Indole:","AUC=","0.84",sep=" "),
                 paste("CI:", "0.73-0.95", sep = " "),
                 paste("Sensitivity:", "89.66%", sep = " "),
                 paste("Speciticity:", "72.89%", sep = " "))
legend("bottomright", 
       legend=legend.name,
       col = mycol[1:4],
       lwd = 2,
       bty="n")
dev.off()




# 使用glm函数出现如下问题：Warning messages:
# 1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred
# 出现上述问题往往是因为自变量已经能够很好的区分因变量了，此时需要使用套索回归等增加了惩罚方式的方法

library(glmnet)
library(plotmo)
data_lasso <- cbind(df_z,outcome)
#data_lasso <- cbind(df,outcome)
##change data in matrix form
x=as.matrix(data_lasso[,-ncol(data_lasso)])
class=as.matrix(data_lasso[,ncol(data_lasso)])

##build a model, set family as multinomial for multinomial logistic regression
lasso.mod=glmnet(x,class,family="binomial", alpha=0, nlambda = 50)
print(lasso.mod)
head(coef(lasso.mod, s=c(lasso.mod$lambda[50],0.009)))
plot_glmnet(lasso.mod, label=5,nresponse = 2) 
## 10 fold CV
cvfit=cv.glmnet(x, class, family="binomial", alpha = 1, nlambda = 1000)
plot(cvfit)
# 两条虚线分别指示了两个特殊的λ值:
c(cvfit$lambda.min,cvfit$lambda.1se) 

fit <- glmnet(x, class, alpha = 1, lambda=cvfit$lambda.1se)
fit <- glmnet(x, class, alpha = 1, lambda=cvfit$lambda.min)
head(fit$beta)
choose_metabo=rownames(fit$beta)[as.numeric(fit$beta)!=0]
# [1] "M1"  "M6"  "M15" "M19" "M27" "M28" "M31" "M33" "M34" "M36" "M37" "M38" "M39" "M45"
# [15] "M49" "M53" "M54" "M55" "M57" "M65" "M72" "M73" "M74" "M77" "M84" "M89" "M93"

keep_metabo <- data_lasso[, choose_metabo]
X <- cbind(keep_metabo, outcome)
summary(glm(as.factor(outcome) ~ ., family = binomial(),data = X))
summary(lm(outcome ~ ., data = X))


