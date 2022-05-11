################
#
#    hora
#    2021.11
#
################

library(GEOquery)

#miRNA
GSE19890_1 <- getGEO(GEO = 'GSE19890',destdir = '.',getGPL = F)
GSE19890_2 <- exprs(GSE19890_1[[1]])
GSE19890_3 <- pData(GSE19890_1[[1]])
GPL9908 <- read.table('GPL9908.txt',header = T,sep = '\t',na.strings = c('NA'),quote = '',fill = T)
GSE19890_2 <- GSE19890_2[,6:15]
GSE19890 <- GSE19890_2[rownames(GSE19890_2) %in% GPL9908$ID, ]
GPL9908_1 <- GPL9908[match(rownames(GSE19890_2),GPL9908$ID),]
rownames(GSE19890) <- GPL9908_1$miRNA_ID_LIST

dimnames <- list(rownames(GSE19890),colnames(GSE19890))
d <- matrix(as.numeric(as.matrix(GSE19890)),nrow=nrow(GSE19890),dimnames=dimnames)
GSE19890 <- as.data.frame(avereps(d))
  #DIFF
library(limma) 
group <- as.data.frame(c(rep('SCI',5),rep('Sham',5)))
rownames(group) <- colnames(GSE19890)
colnames(group) <- 'Group'
design <- model.matrix(~0+factor(group$Group))
colnames(design) <- levels(factor(group$Group))
rownames(design) <- rownames(group)
cm <- makeContrasts('SCI-Sham',levels = design)
fit <- lmFit(GSE19890,design)
fit <- contrasts.fit(fit,cm)
fit <- eBayes(fit)
allDEG <- topTable(fit,coef=1,n=Inf)
allDEG <- na.omit(allDEG)

diff <- GSE19890[match(rownames(allDEG),rownames(GSE19890)),]
GSE19890_selDEG1 <- allDEG[allDEG$P.Value<0.05,]
GSE45006 <- cbind(diff,allDEG)
write.csv(GSE45006,'D:/研究生/25.铁死亡、自噬、焦亡课题/result/miRNA-diff.csv',quote = F)



#circRNA
GSE156999_1 <- getGEO(GEO = 'GSE156999',destdir = '.',getGPL = F)
GSE156999_3 <- pData(GSE156999_1[[1]])
GSE156999_2 <- read.table('GSE156999.txt',header = T,sep = '\t')
GSE156999_anno <- read.table('circRNA.txt',header = T,sep = '\t')
GSE156999 <- GSE156999_2[GSE156999_2$ID %in% GSE156999_anno$ID, ]
GSE156999_anno2 <- GSE156999_anno[match(GSE156999$ID,GSE156999_anno$ID),]
GSE156999 <- as.matrix(GSE156999)
rownames(GSE156999) <- GSE156999_anno2$circRNA
GSE156999 <- GSE156999[,-1]
dimnames <- list(rownames(GSE156999),colnames(GSE156999))
d <- matrix(as.numeric(as.matrix(GSE156999)),nrow=nrow(GSE156999),dimnames=dimnames)
GSE156999 <- as.data.frame(avereps(d))

  #diff
library(edgeR)
data <- GSE156999
group <- factor(c(rep(2,3), rep(1,3)))

y <- DGEList(counts = data,group = group)
head(y)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=F]
y <- calcNormFactors(y)
head(y$samples)
y <- estimateDisp(y)
et <- exactTest(y)

allDEG <- as.data.frame(topTags(et,n=nrow(data)))

GSE156999_selDEG <- allDEG[allDEG$PValue<0.05,]
diff <- data[match(rownames(allDEG),rownames(data)),]
GSE133093 <- cbind(diff,allDEG)
write.csv(GSE133093,'D:/研究生/25.铁死亡、自噬、焦亡课题/circRNA-diff.csv',quote = F)

#pcr
pcr <- read.table('pcr.txt',header = T,sep = '\t')
stat3 <- pcr[c(1,2),]
tlr4 <- pcr[c(3,4),]
hmox1 <- pcr[c(5,6),]
hif1a <- pcr[c(7,8),]
cybb <- pcr[c(9,10),]

library(ggplot2)
ggplot(data = stat3, aes(x=group, y=STAT3,fill=group))+
  geom_bar()+
  geom_jitter()

ggplot(cybb,aes(x=group, y=mean,fill=group)) +
  geom_bar(stat = "identity", position=position_dodge(), width=0.5) + 
  geom_errorbar(aes(ymin=mean-SEM, ymax=mean+SEM), width=0.2, 
                color="black", 
                position=position_dodge(0.6))+
  geom_point()
