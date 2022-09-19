####基于R语言的SparCC分析
rm(list=ls())
#install.packages("robCompositions")
#library(robCompositions)
library(vegan)
library(igraph)
library(qgraph)
library(RColorBrewer)
library(readxl)
library(ggraph)
library(graphlayouts)
library(reshape)
#library(devtools)
#install_github("kasperwelbers/semnet")
library(semnet)
#install.packages('devtools')
#devtools::install_github('zdk123/SpiecEasi')
#加载 SpiecEasi 包
library(SpiecEasi)

#temp_nm1 <- "../otu/ileum/"
#temp_nm2 <- "../otu/cecum/"
#temp_nm3 <- "../otu/fecal/"
#temp_nm4 <- "cecum" # 不同部位时该参数需要改动
# load working datasheet
options (stringsAsFactors = FALSE)
#work_df <- read.csv("Ileum_OTU.csv",header = T,row.names = 1,check.names = F)
a <- read.csv("MSMS/msms_Intensity2_t.csv",
              header = T,row.names = 1,check.names = F)
#overlap_index <- intersect(intersect(rownames(a), rownames(b)), rownames(c))
#meta_df <- b #a,b,c三个数据框根据所需要的肠道部位来改动
#otu <- meta_df
# centered log ratio transformation
# otu.clr <- cenLR(otu, exp(10))
# otu <- decostand(otu,"log")
#otu <- log10(otu)
# load library


#sub_otu <- read.delim("clipboard", header = FALSE)  # micro_asso_metabo_overlap_f_g_s_otu.xlsx，这里就是根据要分析的代谢物来选择，比如在盲肠部位与代谢物相关的Otus是那些，拷贝过来
#keep <- intersect(colnames(otu), sub_otu$V1)        # 和所选择的相应部位总Otus取交集
#C_H.otu <- otu[,keep]                               # 留下有交集的那些Otus，及其丰度表


# sparcc.otu<-sparcc(otu)##Sparcc计算
sparcc.metabo <- sparcc(a, iter = 20, inner_iter = 10, th = 0.1)
cor.metabo <- sparcc.metabo$Cor
#diag(cor.otu) <- 0   #将相关矩阵中对角线中的值（代表了自相关）转为 0
cor.metabo <- as.data.frame(cor.metabo)
rownames(cor.metabo) <- colnames(a)
colnames(cor.metabo) <- colnames(a)
#dir.create("network")
temp_nm1 <- "./network/"
write.table(cor.metabo, paste(temp_nm1, 'sparcc0.txt', sep = ""), sep = '\t', col.names = NA, quote = FALSE)

#通过 100 次自举抽样获取随机矩阵
set.seed(123)
n = 100


for (i in 1:n) {
  amgut1.filt.boot <- sample(a, replace = TRUE)  #bootstrap
  amgut1.filt.sparcc_boot <- sparcc(amgut1.filt.boot, iter = 20, inner_iter = 10, th = 0.1)  #sparcc 参数设置和上文保持一致
  sparcc_boot <- amgut1.filt.sparcc_boot$Cor
  colnames(sparcc_boot) <- colnames(amgut1.filt.boot)
  rownames(sparcc_boot) <- colnames(amgut1.filt.boot)
  write.table(sparcc_boot, paste(temp_nm1,'sparcc_boot', i, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)  #输出随机值的稀疏相关性矩阵
}


#基于上述观测值的稀疏相关性矩阵以及 100 次 bootstrap 的结果，获取稀疏相关性的伪 p 值
p <- cor.metabo
p[p!=0] <- 0

for (i in 1:n) {
  p_boot <- read.delim(paste(temp_nm1,'sparcc_boot', i, '.txt', sep = ''), sep = '\t', row.names = 1)
  p[abs(p_boot)>=abs(cor.metabo)] <- p[abs(p_boot)>=abs(cor.metabo)] + 1
}

p <- p / n
write.table(p, paste(temp_nm1,'pvals.two_sided.txt', sep = ""), sep = '\t', col.names = NA, quote = FALSE)


# 合并相关矩阵和伪p值矩阵，构建网络，以及可视化等
#观测值的相关矩阵
cor_sparcc <- read.delim(paste(temp_nm1, 'sparcc0.txt', sep = ""), row.names = 1, sep = '\t', check.names = FALSE)

#伪 p 值矩阵
pvals <- read.delim(paste(temp_nm1, 'pvals.two_sided.txt', sep = ""), row.names = 1, sep = '\t', check.names = FALSE)


#保留 |相关性|>0.3 且 p<0.05的值（不满足条件的相关性均赋值为 0）
cor_sparcc[abs(cor_sparcc)<=0.3 | pvals>=0.05] <- 0
pvals[abs(cor_sparcc)<=0.3 | pvals>=0.05] <- 1

cor_sparcc2 <- data.frame(row = rownames(cor_sparcc), cor_sparcc, check.names = F)
pvals2 <- data.frame(row = rownames(pvals), pvals, check.names = F)
nbar.m <- melt(cor_sparcc2)
nbap.m <- melt(pvals2)
# 将nbap.m的P值列加到nbar.m中去
nbar.m2 <- cbind(nbar.m, nbap.m$value)
names(nbar.m2)[3:4] <- c("Coef", "Pvalue")
# 删除相关系数为0的行
nbar.m2 <- nbar.m2[-which(nbar.m2$Coef == 0), ]
nbar.m2 <- nbar.m2[which(nbar.m2$Pvalue != 1), ]
dim(nbar.m2)
write.csv(nbar.m2, "association_metabolites20220906.csv")


#将相关矩阵中对角线中的值（代表了自相关）转为 0
diag(cor_sparcc) <- 0

#输出邻接矩阵类型的网络文件
write.table(cor_sparcc, paste(temp_nm1, 'neetwork.adj.txt', sep = ""), col.names = NA, sep = '\t', quote = FALSE)

## 开始用igraph构建网络
occor = cor_sparcc
# 确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
# 粪便将相关系数的阈值设置为0.3
#occor[abs(occor) < 0.3] = 0

# 构建igraph对象
igraph = graph_from_adjacency_matrix(as.matrix(occor),mode="undirected",weighted=TRUE,diag=FALSE)
igraph

# remove isolated nodes，即去掉和所有otu均无相关性的otu 可省略，前期矩阵已处理过
bad.vs = V(igraph)[degree(igraph) == 0]
igraph = delete.vertices(igraph, bad.vs)
igraph

# 将igraph weight属性赋值到igraph.weight
igraph.weight = E(igraph)$weight

# 做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
E(igraph)$weight = NA

# 简单出图
# 设定随机种子数，后续出图都从同一随机种子数出发，保证前后出图形状相对应
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=FALSE,margin=c(0,0,0,0))

#2.按相关类型设置边颜色
# 如果构建网络时，weighted=NULL,此步骤不能统计
sum(igraph.weight>0)# number of postive correlation
sum(igraph.weight<0)# number of negative correlation

# set edge color, postive correlation 设定为red, negative correlation设定为blue
E.color = igraph.weight
E.color = ifelse(E.color>0, "#D65DB1",ifelse(E.color<0, "#2C73D2","grey"))
E(igraph)$color = as.character(E.color)

# 改变edge颜色后出图
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=FALSE,margin=c(0,0,0,0))
#3.按相关性设置边宽度
# 可以设定edge的宽 度set edge width，例如将相关系数与edge width关联
#E(igraph)$width = abs(igraph.weight)*20
E(igraph)$Coef. = abs(igraph.weight)

metabo_pro <- read.csv("Table_1.csv", header = T, row.names = 1, check.names = F)
igraph.size = metabo_pro[V(igraph)$name,9] # 筛选对应OTU属性
#colnames(igraph.size) <- "abundance"
#igraph.size1 = log(igraph.size*100) # 原始数据是什么，为什么*100再取e对数
V(igraph)$VIP = igraph.size

#' 2.按模块着色
# 模块性 modularity
fc = cluster_fast_greedy(igraph,weights =NULL)# cluster_walktrap cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity = modularity(igraph,membership(fc))
# 按照模块为节点配色
comps = membership(fc)
# 每个模块otus数目的多少
module_size = sizes(fc)
module_size
#colbar = rainbow(max(comps))
#colbar = c("gold","deeppink","deepskyblue","yellow","brown","pink")
#colbar = c("gold","deeppink","deepskyblue", "yellow","brown","pink","indianred1")
# display.brewer.all(type = "qual") 列出所有调色板，供选择
colbar = brewer.pal(7,"Dark2")
modbar = c(paste("module",seq(1,7,1),sep = " "))
V(igraph)$color = colbar[comps]
V(igraph)$modulname = modbar[comps]
modulcolor = colbar[as.numeric(as.factor(V(igraph)$modulname))]

#' 3.按照组间上调和下调进行配色
V(igraph)$Group = metabo_pro[V(igraph)$name,10]
colbar2 = c("#E16663", "#90CBD3")
my_color = colbar2[as.numeric(as.factor(V(igraph)$Group))]

#' 4.把每个
#' 
#' 
#' 
#' 
#' 
# 用ggplot2画网络图参考来源：https://mr.schochastics.net/material/netvizr/
p3 <- ggraph(igraph, layout = "auto") +
  geom_edge_link0(aes(edge_width = Coef.), edge_colour = E(igraph)$color) +
  geom_node_point(aes(fill = Group, size = VIP, shape = Group)) +
  #geom_node_point(aes(fill = modulname, size = size),shape = 21) +
  geom_node_text(aes(label = name, size = VIP), family = "serif",
                 repel = TRUE) +
  scale_shape_manual(values = c(21, 22)) + 
  scale_fill_manual(values = my_color) +
  scale_edge_width(range = c(1, 2)) +
  scale_size(range = c(2, 6)) +
  theme_graph() +
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size = 15))
ggsave("metabo_network20220511.png", width = 16, height = 8, dpi = 300)

p4 <- ggraph(igraph, layout = "stress") +
  geom_edge_link0(aes(edge_width = Coef.), edge_colour = E(igraph)$color) +
  geom_node_point(aes(fill = modulname, size = VIP), shape = 21) +
  #geom_node_point(aes(fill = modulname, size = size),shape = 21) +
  geom_node_text(aes(label = name, size = VIP), family = "serif",
                 repel = TRUE) +
  #scale_shape_manual(values = c(21, 22)) + 
  #scale_fill_manual(values = modulcolor) +
  scale_edge_width(range = c(1, 2)) +
  scale_size(range = c(2, 6)) +
  theme_graph() +
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size = 15))
ggsave("metabo_network_module20220511.png", width = 16, height = 8, dpi = 300)
library(cowplot)
p5 <- plot_grid(p3,p4,labels = c("A", "B"), nrow = 2)
ggsave("metabo_network_merge0906.tiff",plot = p5, width = 12, height = 12, dpi = 300)






# note: we choose the hub Otu(s) based on at least two parameters from 'degree', 'closeness'
# or 'betweenness'.
# degree
mostcloseness_igraph <- names(closeness(igraph))[which(as.vector(closeness(igraph)) == max(closeness(igraph)))]
mostcloseness_igraph

mostdegree_igraph <- names(degree(igraph))[which(as.vector(degree(igraph)) == max(degree(igraph)))]
mostdegree_igraph

mostbetweenness_igraph <- names(betweenness(igraph))[which(as.vector(betweenness(igraph)) == max(betweenness(igraph)))]
mostbetweenness_igraph
# obtain the hub Otu 
most_hub <- intersect(intersect(mostcloseness_igraph, mostdegree_igraph),
                      mostbetweenness_igraph) # 有时候三个条件不能同时满足，那就看满足任意两个条件的核心Otu是哪个
most_hub  
hub_score(igraph)$vector
# 获得hub metabolites的名字
names(hub_score(igraph)$vector)[as.vector(hub_score(igraph)$vector) == 1.0]
# D-Galactose
as.vector(comps[names(comps) == "D-Galactose"]) # 快速确定核心Otu是哪个模块的

#radian.rescale <- function(x, start=0, direction=1) {
#  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
#  c.rotate(ggplot2::rescale(x, c(0, 2 * pi), range(x)))
#}
#n <- length(V(igraph)$name)
#lab.locs <- radian.rescale(x=n, direction=1, start=0)

tiff("metabolite_network_module.tiff",width = 4000,height = 4000, res = 300)
#set.seed(100)
#par(mar=c(2,2,2,4))
e <- get.edgelist(igraph, names = FALSE)
#l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(subgraph_sub)) #https://stackoverflow.com/questions/39290909/igraph-resolving-tight-overlapping-nodes
#l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(subgraph_sub))
ll <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(igraph),
                                        area=8*(vcount(igraph)^2),repulse.rad=(vcount(igraph)^3.2))
#l <- layout_with_fr(igraph)
#ll <- norm_coords(ll, ymin = -1,ymax = 1,xmin = -1,xmax = 1)
plot(igraph,
     rescale=TRUE,
     layout=ll, 
     vertex.label.dist=0.1, 
     #vertex.label=NA,
     #vertex.color = my_color,
     vertex.label.color="black", 
     vertex.label.degree=0, 
     edge.curved=FALSE,
     vertex.shape=ifelse(V(igraph)$name %in% "D-Galactose", "square","circle"),
     vertex.label.cex=ifelse(V(igraph)$name %in% "D-Galactose", 1.6, 1.2),
     vertex.label.font=1, 
     edge.width=E(igraph)$width, 
     vertex.frame.color=NA,
     margin=c(0,0,0,0))
title("Co-occurrence network", cex.main=1.8)
library(stringr)
text(0.5, 1.1, labels = "Hub metabolites: D-Galactose", cex=1.5)
legend("topleft",legend = levels(as.factor(unique(V(igraph)$modulname))), 
       fill=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462"),bty = "n",cex=1.4,title = "Module names")
#legend("topleft",legend = levels(as.factor(unique(V(igraph)$updown))), 
#       col = colbar2,bty = "n", pch = 20, cex=1.4, text.col = colbar2)


#legend("topleft",legend = paste("module", seq(1,11,1), sep = ""), 
#       fill=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
#              "#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5"),bty = "n",cex=1.2,title = "Module names")
#E.color = ifelse(E.color>0, "lightcoral",ifelse(E.color<0, "lightgreen","grey"))
#legend(-1.2, -1.1,legend = c("module1"), col = c("lightcoral"), lty = 1, lwd = 3, 
#      cex = 1.4,title = "Orientation Of module to Indosyl sulfate", bty = "n")
dev.off()





################### step2 取核心Otu所在的模块 ###################
###                                                          ####
#################################################################
# 将某个模块抽出来
module_sub <- fc[2]$`2` # 将属于模块XX的otu选出来
subgraph_sub <- induced_subgraph(igraph, module_sub) # 提取属于module1的网络节点图

fc2 = cluster_fast_greedy(subgraph_sub,weights =NULL)# cluster_walktrap cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity2 = modularity(subgraph_sub,membership(fc2))
# 按照模块为节点配色
comps2 = membership(fc2)
# 每个模块otus数目的多少
module_size2 = sizes(fc2)
module_size2
# display.brewer.all(type = "qual") 列出所有调色板，供选择
colbar2 = brewer.pal(3,"Set2")
modbar2 = c(paste("sub_module",seq(1,2,1),sep = ""))
V(subgraph_sub)$color = colbar2[comps2]
V(subgraph_sub)$modulname = modbar2[comps2]

#V(subgraph_sub)$color <- ifelse(V(subgraph_sub)$name == 'Otu36', "gold", "lightgoldenrod1")
#V(subgraph_sub)$color <- ifelse(V(subgraph_sub)$name %in% c('Otu38', 'Otu238'), "deepskyblue", "lightblue")
#V(subgraph_sub)$size <- ifelse(V(subgraph_sub)$name == 'Otu238', 5, 2)

tiff(paste(temp_nm5, "_otu_network_3-Indole carboxylic acid glucuronide_sub_module2.tiff", sep = ""),width = 4000,height = 4000,res = 300)
#set.seed(123)
#par(mar=c(2,2,2,4))
#l <- layout_with_fr(subgraph_sub)
#l <- layout_in_circle(subgraph_sub)
#l <- norm_coords(l, ymin = -1,ymax = 1,xmin = -1,xmax = 1)
e <- get.edgelist(subgraph_sub, names = FALSE)
#l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(subgraph_sub)) #https://stackoverflow.com/questions/39290909/igraph-resolving-tight-overlapping-nodes
#l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(subgraph_sub))
ll <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(subgraph_sub),
                                        area=10*(vcount(subgraph_sub)^2),repulse.rad=(vcount(subgraph_sub)^3.1))
#ll <- as.matrix(ll)
plot(subgraph_sub,rescale=T,layout=ll, vertex.label.cex=ifelse(V(subgraph_sub)$name %in% c('Otu238'), 1.6,1.2), 
     vertex.label.dist=0.5, vertex.shape=ifelse(V(subgraph_sub)$name %in% c('Otu238'), "square","circle"),
     vertex.frame.color=FALSE,vertex.label.color="black",
     vertex.label.degree=-0.6,edge.curved=FALSE,
     vertex.label.font=1)
text(0.5, 1.1, labels = "Hub Otu: Otu238 (Prevotella copri_N)", cex=1.5)
title("Co-occurrence network with 3-Indole carboxylic acid glucuronide under Module2", cex.main=1.8)
legend("topleft",legend = levels(as.factor(unique(V(subgraph_sub)$modulname))),
       fill=c("#66C2A5","#FC8D62"),bty = "n",cex=1.4,title = "Sub_module names")
legend(-1.2, -1,legend = c("sub_module1 (Otu768 excluded)"), col = c("lightcoral"), lty = 1, lwd = 3, 
       cex = 1.4,title = "Orientation Of module to LCA-GC", bty = "n")
dev.off()




################## step3 可以将网络属性储存到相应文件中 ############################
#####                                                                         ######
################# 用于其它可视化软件的输入文件          ############################

#这种转换模式下，默认的边权重代表了 sparcc 计算的相关性（存在负值）
#由于边权重通常为正值，因此最好取个绝对值，相关性重新复制一列作为记录
E(igraph)$sparcc <- E(igraph)$weight
E(igraph)$weight <- abs(E(igraph)$weight)

#再转为其它类型的网络文件，例如
#边列表，包括节点相关性的信息、边的权重系数等
edge <- data.frame(as_edgelist(igraph))

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(igraph)$weight,
  sparcc = E(igraph)$sparcc
)
head(edge_list)

write.table(edge_list, paste(temp_nm5, 'network.edge_list.txt',sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表，对应边列表，记录节点属性，例如
node_list <- data.frame(
  id = V(igraph)$name,
  label = V(igraph)$name,
  degree = degree(igraph),
  betweenness = betweenness(igraph),
  closness = closeness(igraph)
  
)
head(node_list)

write.table(node_list, paste(temp_nm5, 'network.node_list.txt',sep = ""), sep = '\t', row.names = FALSE, quote = FALSE)

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(igraph, paste(temp_nm5, 'network.graphml', sep = ""), format = 'graphml')

