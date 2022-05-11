################
#
#    hora
#    2021.11
#
################

#数据下载----
library(GEOquery)
GSE45006_1 <- getGEO(GEO = "GSE45006", destdir = ".", getGPL = F)
GSE45006_2 <- exprs(GSE45006_1[[1]])
GSE45006_3 <- pData(GSE45006_1[[1]])
data <- GSE45006_2[,c(5:8,17:20)]

#分组信息----
library(stringr)
GSE45006_group <- as.data.frame(str_split(GSE45006_3$source_name_ch1,',',simplify = T)[,2])
colnames(GSE45006_group) <- 'group'
rownames(GSE45006_group) <- rownames(GSE45006_3)
group <- as.data.frame(GSE45006_group[c(5:8,17:20),])
colnames(group) <- 'group'
rownames(group) <- rownames(GSE45006_group)[c(5:8,17:20)]

#id trans----
library(rat2302.db)
library(limma)
ids <- toTable(rat2302SYMBOL)

length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))

data <- as.data.frame(data)
GSE45006 <- data[rownames(data) %in% ids$probe_id,]
ids2 <- ids[match(rownames(GSE45006),ids$probe_id),]
GSE45006 <- as.matrix(GSE45006)
rownames(GSE45006) <- ids2$symbol
dimnames <- list(rownames(GSE45006),colnames(GSE45006))
d <- matrix(as.numeric(as.matrix(GSE45006)),nrow=nrow(GSE45006),dimnames=dimnames)
dat <- as.data.frame(avereps(d))


#pca----
library(ggplot2)

pc <- prcomp(t(dat))
pc1 <- as.data.frame(pc$x)
pc1$group <- c(rep('Sham',4),rep('SCI',4))
percentage <- round(pc$sdev / sum(pc$sdev) * 100, 2)
percentage <- paste(colnames(pc1), "(", paste( as.character(percentage), "%", ")", sep="") )
ggplot(pc1,aes(x=PC1, y=PC2,color=group))+
  xlab(percentage[1])  +
  ylab(percentage[2])  +
  stat_ellipse(aes(shape = NULL))+
  geom_point(size=2)

#boxplot
library(ggpubr)
rt <- t(dat)
type <- c(rep('sham',4),rep('SCI',4))
boxdata <- data.frame()
for(i in colnames(rt)){
  boxdata <- rbind(boxdata,cbind(expression=rt[,i],gene=i,type))
}
boxdata$sample <- rep(colnames(dat),14735)

ggboxplot(boxdata, x="sample", y="expression", color = "type", 
            ylab="Gene expression",
            xlab="",
            palette = c("green","red") )+rotate_x_text(60)
ggplot(data = boxdata,aes(x="sample", y="expression"))+
  geom_boxplot()
boxplot(dat,las=2)
