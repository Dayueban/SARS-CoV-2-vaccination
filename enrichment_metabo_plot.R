# dotplot or barplot of enrichment analysis for metabolomics
rm(list = ls())
library(argparser)
library(ggplot2)
library(ggpubr)
a <- read.csv("msea_enrichment/msea_VB2SV_down.csv", header = T, check.names = F)
a$Enrichment_ratio <- a$hits / a$expected

a <- a[order(a$Enrichment_ratio, decreasing = T), ]

# showing top 20 metabolic pathway
a <- a[1:20, ]

#Drawing enrichment scatter plot
p <- ggplot(a, aes(x = reorder(pathway, Enrichment_ratio), y = Enrichment_ratio, colour = FDR, size = Enrichment_ratio))
p <- p + geom_point() +
  coord_flip()
p <- p + scale_colour_gradientn(colours = rainbow(4), guide = "colourbar") + expand_limits(color = seq(0,1,by = 0.25))
p <- p + ggtitle("Enriched Metabolites Sets in SV") + xlab("Pathway") +ylab("Enrichment ratio")
p <- p + theme_bw() + theme(axis.text = element_text(color = "black", size = 10),
                            axis.title = element_text(color = "black", size = 14))
p <- p + theme(panel.border = element_rect(colour = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.key = element_blank())
ggsave("msea_enrichment/msea_VB2SV_down.png", plot = p, width = 8, height = 6, type = 'cairo-png')


# second
b <- read.csv("msea_enrichment/msea_VB2SV_up.csv", header = T, check.names = F)
b$Enrichment_ratio <- b$hits / b$expected

b <- b[order(b$Enrichment_ratio, decreasing = T), ]

# showing top 20 metabolic pathway
b <- b[1:20, ]

#Drawing enrichment scatter plot
p1 <- ggplot(b, aes(x = reorder(pathway, Enrichment_ratio), y = Enrichment_ratio, colour = FDR, size = Enrichment_ratio))
p1 <- p1 + geom_point() +
  coord_flip()
p1 <- p1 + scale_colour_gradientn(colours = rainbow(4), guide = "colourbar") + expand_limits(color = seq(0,1,by = 0.25))
p1 <- p1 + ggtitle("Enriched Metabolites Sets in VB") + xlab("Pathway") +ylab("Enrichment ratio")
p1 <- p1 + theme_bw() + theme(axis.text = element_text(color = "black", size = 10),
                            axis.title = element_text(color = "black", size = 14))
p1 <- p1 + theme(panel.border = element_rect(colour = "black"))
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5), legend.key = element_blank())
ggsave("msea_enrichment/msea_VB2SV_up.png", plot = p1, width = 8, height = 6, type = 'cairo-png')

p2 <- ggarrange(p, p1, ncol = 2, labels = c("A", "B"))
ggsave("msea_enrichment/msea.png", plot = p2, width = 16, height = 6)
