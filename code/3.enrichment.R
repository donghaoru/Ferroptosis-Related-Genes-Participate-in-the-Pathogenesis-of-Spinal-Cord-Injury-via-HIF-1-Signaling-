################
#
#    hora
#    2021.11
#
################

library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)
library(GOplot)
library(circlize)
library(enrichplot)
#DEG ----
DEG.id <- mapIds(x=org.Rn.eg.db,  
                 keys=rownames(selDEG),
                 keytype ='SYMBOL',
                 column='ENTREZID')
DEG.id <- na.omit(DEG.id)
KEGG <- enrichKEGG(gene = DEG.id,
                   organism = "rno",
                   pvalueCutoff =1,
                   qvalueCutoff = 1)
dotplot(KEGG, showCategory = 15)

GO <- enrichGO(gene = DEG.id,
               OrgDb = org.Rn.eg.db,
               pvalueCutoff =1,
               ont="all",
               readable =T,
               qvalueCutoff = 1)
dotplot(GO,showCategory = 5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')

k <- as.data.frame(KEGG)
kk <- data.frame(Category = "All",ID = k$ID,Term = k$Description, Genes = gsub("/", ", ", k$geneID), adj_pval = k$p.adjust)

a <- allDEG
eg <- bitr(geneID = rownames(a),
           fromType = 'SYMBOL',
           toType = 'ENTREZID',
           OrgDb = org.Rn.eg.db)
a$ENTREZID <- eg[match(rownames(a),eg$SYMBOL),2]
a <- na.omit(a[order(a$logFC,decreasing = T),])
genelist <- a$logFC
names(genelist) <- a$ENTREZID
GSEA <- gseKEGG(genelist,organism = 'rno',pvalueCutoff = 1)
dotplot(GSEA,showCategory = 15,split=".sign")+facet_wrap(~.sign,scales = "free")
path <- 'rno04066'
gseaplot2(GSEA,path,pvalue_table=T)

barplot(GO,color='pvalue', drop = TRUE, showCategory =5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
barplot(KEGG,color='pvalue', showCategory = 15)
#visualize----
ego=as.data.frame(top_KEGG)
ego <- ego[1:5,]
go=data.frame(Category = "All",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)

genelist <- data.frame(ID = toupper(rownames(topDEG)), logFC = topDEG$logFC)
genelist <- data.frame(ID = DEG.id, logFC = select$logFC)
genelist <- data.frame(ID = rownames(topDEG), logFC = topDEG$logFC)
row.names(genelist)=genelist[,1]
a <- go[c(1,5,10,15),]
circ <- circle_dat(go,genelist)
termNum = 4                                     
geneNum = nrow(genelist)                        

chord <- chord_dat(circ,genelist[1:geneNum,], go$Term[1:termNum])

GOChord(chord, 
        space = 0.001,           #基因之间的间距
        gene.order = 'logFC',    #按照logFC值对基因排序
        gene.space = 0.25,       #基因名跟圆圈的相对距离
        gene.size = 5,           #基因名字体大小 
        border.size = 0.1,       #线条粗细
        process.label = 7.5)     #term字体大小
#select-----
select <- selDEG[gene,]
DEG.id <- mapIds(x=org.Rn.eg.db,  
                    keys=rownames(select),
                    keytype ='SYMBOL',
                    column='ENTREZID')
DEG.id <- na.omit(DEG.id)
sel_KEGG <- enrichKEGG(gene = DEG.id,
                   organism = "rno",
                   pvalueCutoff =0.05,
                   qvalueCutoff = 1)
dotplot(sel_KEGG, showCategory = 15)

sel_GO <- enrichGO(gene = DEG.id,
               OrgDb = org.Rn.eg.db,
               pvalueCutoff =0.05,
               ont="all",
               readable =T,
               qvalueCutoff = 1)
dotplot(sel_GO,showCategory = 5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')

a <- select
eg <- bitr(geneID = rownames(a),
           fromType = 'SYMBOL',
           toType = 'ENTREZID',
           OrgDb = org.Rn.eg.db)
a$ENTREZID <- eg[match(rownames(a),eg$SYMBOL),2]
a <- na.omit(a[order(a$logFC,decreasing = T),])
genelist <- a$logFC
names(genelist) <- a$ENTREZID
sel_GSEA <- gseKEGG(genelist,organism = 'rno',pvalueCutoff = 1)
dotplot(sel_GSEA,showCategory = 15,color='pvalue',split=".sign")+facet_wrap(~.sign,scales = "free")


#top10
DEG.id <- mapIds(x=org.Rn.eg.db,  
                 keys=g,
                 keytype ='SYMBOL',
                 column='ENTREZID')
DEG.id <- na.omit(DEG.id)
top_KEGG <- enrichKEGG(gene = DEG.id,
                       organism = "rno",
                       TERM2GENE = DEG.id,
                       TERM2NAME = g,
                       pvalueCutoff =1,
                       qvalueCutoff = 1)

de_ekp <- enricher(DEG.id,
                   TERM2GENE = DEG.id,
                   TERM2NAME = g,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1)

dotplot(top_KEGG, showCategory = 15)

top_GO <- enrichGO(gene = DEG.id,
                   OrgDb = org.Rn.eg.db,
                   pvalueCutoff =1,
                   ont="all",
                   readable =T,
                   qvalueCutoff = 1)
dotplot(top_GO,showCategory = 5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')

top_DEG <- allDEG[g,]
a <- top_DEG
eg <- bitr(geneID = rownames(a),
           fromType = 'SYMBOL',
           toType = 'ENTREZID',
           OrgDb = org.Rn.eg.db)
a$ENTREZID <- eg[match(rownames(a),eg$SYMBOL),2]
a <- na.omit(a[order(a$logFC,decreasing = T),])
genelist <- a$logFC
names(genelist) <- a$ENTREZID
top_GSEA <- gseKEGG(genelist,organism = 'rno',pvalueCutoff = 1)
dotplot(sel_GSEA,showCategory = 15,color='pvalue',split=".sign")+facet_wrap(~.sign,scales = "free")

fc <- topDEG$logFC
names(fc) <- rownames(topDEG)
cnetplot(top_GO, circular = TRUE, colorEdge = TRUE,foldChange = fc,showCategory = 4)

top_KEGG <- setReadable(top_KEGG, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
cnetplot(top_KEGG, circular = TRUE, colorEdge = TRUE,foldChange = fc,showCategory = 4)


save(allDEG,dat,GO,group,GSEA,KEGG,selDEG,sel_GO,sel_KEGG,sel_GSEA,top_GSEA,top_KEGG,top_GO,gene,g,file = 'all.Rdata')
save(allDEG,GSE156999,GSE156999_selDEG,file = 'circRNA.Rdata')
save(allDEG,GSE19890,group,GSE19890_selDEG,file = 'miRNA.Rdata')
#Stat3 Tlr4 Hmox1 Hif1a Cybb