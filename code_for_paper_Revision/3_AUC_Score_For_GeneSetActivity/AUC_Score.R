library(Seurat)
library(pheatmap)
library(AUCelll)
library(GSEABase)


MSC=readRDS('/data/gaoyuchen/MSC_no_C8_C567re-Order.rds')

cells_rankings <- AUCell_buildRankings(MSC@assays$RNA@data)

##load gene set, e.g. same geneset used in GSVA analysis
geneset <- getGmt("/data/gaoyuchen/NC_Revision_AUCScore/AUC_GeneSets.gmt")  

##Calculate AUC
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings,aucMaxRank = ceiling(0.1*nrow(cells_rankings)))

##Calculate Average Matrix
AUC=as.data.frame(cells_AUC@assays@data@listData$AUC)

meta=MSC@meta.data
C1=row.names(subset(meta,Cluster=='C1'))
C2=row.names(subset(meta,Cluster=='C2'))
C3=row.names(subset(meta,Cluster=='C3'))
C4=row.names(subset(meta,Cluster=='C4'))
C5=row.names(subset(meta,Cluster=='C5'))
C6=row.names(subset(meta,Cluster=='C6'))
C7=row.names(subset(meta,Cluster=='C7'))

AUC=cbind(rowMeans(AUC[,C1]),rowMeans(AUC[,C2]),rowMeans(AUC[,C3]),rowMeans(AUC[,C4]),rowMeans(AUC[,C5])
          ,rowMeans(AUC[,C6]),rowMeans(AUC[,C7]))

##Heatmap for Fig 1c
pheatmap(AUC,cluster_rows = F,cluster_cols = F,scale = 'row',colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),border=FALSE,cellwidth = 10, cellheight = 10,breaks=unique(c(seq(-2,2, length=100))))

write.csv(AUC,'/data/gaoyuchen/NC_Revision_AUCScore/AUC_Score_result.csv')
