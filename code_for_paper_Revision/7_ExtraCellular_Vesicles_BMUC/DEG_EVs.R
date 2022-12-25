library(DESeq2)
library(ggplot2)

group=read.csv('/data/Manuscript for NM/20220826 NC revision/Exosome BM vs UC/group.csv',row.names = 'sample')
RawData_Exosome=read.csv('/data/Manuscript for NM/20220826 NC revision/Exosome BM vs UC/RawData-BM_UCExosome.csv',row.names = 'SYMBOL')


row.names(group)=colnames(RawData_Exosome)
dds_EXO <- DESeqDataSetFromMatrix(countData = as.matrix(RawData_Exosome),
                                  colData = as.matrix(group),
                                  design= ~ UC_vs_Other)

dds_EXO_UC <- DESeq(dds_EXO)
res_EXO_UC <- as.matrix(results(dds_EXO_UC))


data <- data.frame(logFC=res_EXO_UC$log2FoldChange,padj=res_EXO_UC$padj)
row.names(data)=row.names(res_EXO_UC)
data$sig[(data$padj > 0.05|is.na(data$padj))|(data$logFC < 1)& data$logFC > -1] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 1] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -1] <- "down"
write.csv(data,'/data/Manuscript for NM/20220826 NC revision/Exosome BM vs UC/Volcanoplot.csv')
data=read.csv('/data/Manuscript for NM/20220826 NC revision/Exosome BM vs UC/Volcanoplot.csv',row.names='Gene')

library(ggrepel)
ggplot(data,aes(logFC,-1*log10(padj),color = sig))+geom_point(size=4)+labs(x="log2(FoldChange)",y="-log10(P.adj)")+ scale_color_manual(values =c("#BC3C28","grey","#0072B5"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+geom_text_repel(aes(x=logFC,y=-1*log10(padj),label=label),max.overlaps = 100)

