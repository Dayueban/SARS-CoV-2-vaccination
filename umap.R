# install UMAP package

install.packages("umap")
library(umap)

# Example
head(iris, 3)
iris.data = iris[, grep("Sepal|Petal", colnames(iris))]
iris.labels = iris[, "Species"]
# Let’s load the umap package and apply the UMAP transformation
#library(umap)
iris.umap = umap(iris.data)
head(iris.umap$layout, 3)
plot.iris(iris.umap, iris.labels)


# real dataset in negative mode
a <- read.csv("MS/negative/metabolome_neg.csv", header = T, row.names = 1, 
              check.names = F)
a.data <- a[,-1]
a.labels <- as.factor(a[,1])
a.umap <- umap(a.data)
a.axis <- as.data.frame(a.umap$layout)
a.axis <- cbind(a.axis, a.labels)
names(a.axis) <- c("umap1", "umap2", "group")
pdf("MS/negative/umap.pdf", width = 8, height = 6)
plot.iris(a.umap, a.labels)
dev.off()
# by ggplot2
library(ggplot2)
p <- ggplot(a.axis, aes(x=umap1, y=umap2, group=group)) +
  geom_point(aes(shape=group, color=group, size=group))+
  scale_shape_manual(values=c(16, 17, 18))+
  scale_color_manual(values=c("#ff7f00", "#e377c2", "#17becf"))+
  scale_size_manual(values=c(2,3,4))+
  theme(legend.position="top")

p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(color="black",size=14, face="bold"),
               axis.text.y=element_text(size = 16, colour = "black", face="bold"))
p <- p + guides(fill=guide_legend(title=NULL))
p <- p + theme(axis.line.x=element_line(colour = "black", size = 1),axis.line.y=element_line(colour = "black", size = 1),
               legend.key.height=unit(1,"cm"),
               legend.key.width=unit(1,"cm"),
               legend.text=element_text(lineheight=0.8,face="bold",size=16),
               legend.title=element_text(size=20))
ggsave("MS/negative/umap_n.pdf", width = 6, height = 8)


# real dataset in positive mode
b <- read.csv("MS/positive/metabolome_pos.csv", header = T, row.names = 1,
              check.names = F)
b.data <- b[,-1]
b.labels <- as.factor(b[,1])
b.umap <- umap(b.data)
b.axis <- as.data.frame(b.umap$layout)
names(b.axis) <- c("umap1", "umap2")
b.axis <- cbind(b.axis, a.labels)
names(b.axis)[3] <- "group"
pdf("MS/positive/umap_P.pdf", width = 8, height = 6)
plot.iris(b.umap, b.labels)
dev.off()


library(ggplot2)
p <- ggplot(b.axis, aes(x=umap1, y=umap2, group=group)) +
  geom_point(aes(shape=group, color=group, size=group))+
  scale_shape_manual(values=c(16, 17, 18))+
  scale_color_manual(values=c("#ff7f00", "#e377c2", "#17becf"))+
  scale_size_manual(values=c(2,3,4))+
  theme(legend.position="top")

p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(color="black",size=14, face="bold"),
               axis.text.y=element_text(size = 16, colour = "black", face="bold"))
p <- p + guides(fill=guide_legend(title=NULL))
p <- p + theme(axis.line.x=element_line(colour = "black", size = 1),axis.line.y=element_line(colour = "black", size = 1),
               legend.key.height=unit(1,"cm"),
               legend.key.width=unit(1,"cm"),
               legend.text=element_text(lineheight=0.8,face="bold",size=16),
               legend.title=element_text(size=20))
ggsave("MS/positive/umap_p2.pdf", width = 6, height = 8)



# function of plot.iris
plot.iris <- function(x, labels,
                      main="A UMAP visualization of the Metabolomic dataset",
                      colors=c("#ff7f00", "#e377c2", "#17becf"),
                      pad=0.1, cex=0.8, pch=19, add=FALSE, legend.suffix="",
                      cex.main=1.2, cex.legend=1.2) {
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  } 
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  labels.u = unique(labels)
  legend.pos = "topleft"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomleft"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}


