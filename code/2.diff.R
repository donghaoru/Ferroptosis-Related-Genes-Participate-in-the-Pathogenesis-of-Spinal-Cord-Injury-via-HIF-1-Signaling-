################
#
#    hora
#    2021.11
#
################

group$group <- c(rep('Sham',4),rep('SCI',4))

library(limma)
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- rownames(group)
cm <- makeContrasts('SCI-sham',levels = design)
fit <- lmFit(dat,design)
fit <- contrasts.fit(fit,cm)
fit <- eBayes(fit)
allDEG <- topTable(fit,coef=1,n=Inf)
allDEG <- na.omit(allDEG)

diff <- dat[match(rownames(allDEG),rownames(dat)),]
selDEG <- topTable(fit,p.value = 0.05,lfc = 1,n=Inf)
GSE45006 <- cbind(diff,allDEG)
write.csv(GSE45006,'D:/研究生/25.铁死亡、自噬、焦亡课题/result/GSE45006.csv',quote = F)

#heatmap
expro <- dat[gene,]
expro <- scale(expro)
anno <- group
library(pheatmap)
pheatmap(expro,annotation_col = anno)

expro <- GSE156999[rownames(GSE156999_selDEG),]
colnames(expro) <- GSE156999_3$geo_accession
expro <- scale(expro)
anno <- as.data.frame(c(rep('SCI',3),rep('Sham',3)))
rownames(anno) <- colnames(expro)
colnames(anno) <- 'Group'
library(pheatmap)
pheatmap(expro,annotation_col = anno)

expro <- GSE19890[rownames(GSE19890_selDEG1),]
rownames(expro) <- rownames(GSE19890_selDEG)
expro <- scale(expro)
anno <- group

library(pheatmap)
pheatmap(expro,annotation_col = anno,cluster_cols = F)


#volcano
library(ggpubr)
library(ggthemes)

dat <- allDEG
dat$Group = as.factor(ifelse(dat$P.Value < 0.05 & dat$logFC <= -log2(1), 'Down',
                                 ifelse(dat$P.Value< 0.05 & dat$logFC >= log2(1), 'Up', 'Not')))
dat$v=-log2(dat$P.Value)

ggscatter(dat, legend = "right",#ylim = c(0, 75), xlim=c(-5,5),
          x = "logFC", y = "v",
          xlab = 'LogFC', ylab="-Log2(Pvalue)",
          color = "Group", size = 2, 
          palette = c("#00AFBB", "#999999", "#FC4E07")
        )+theme_base()+
  geom_hline(yintercept = -log2(0.05),linetype='dashed')+
  geom_vline(xintercept = c(-log2(1.1),log2(1.1)),linetype='dashed')


dat <- allDEG
dat$Group = as.factor(ifelse(dat$adj.P.Val < 0.05 & dat$logFC <= -log2(2), 'Down',
                             ifelse(dat$adj.P.Val< 0.05 & dat$logFC >= log2(2), 'Up', 'Not')))
dat$v=-log2(dat$adj.P.Val)
dat$X <- rownames(dat)

ggscatter(dat, legend = "right",#ylim = c(0, 75), xlim=c(-5,5),
          x = "logFC", y = "v",
          xlab = 'LogFC', ylab="-Log2(adj.Pvalue)",
          color = "Group", size = 1, 
          palette = c("#00AFBB", "#999999", "#FC4E07"),
          label = 'X',
          repel = T,
          label.select =c('Stat3','Tlr4','Hif1a','Hmox1','Cybb')
           )+theme_base()+
  geom_hline(yintercept = -log2(0.05),linetype='dashed')+
  geom_vline(xintercept = c(-log2(2),log2(2)),linetype='dashed')
