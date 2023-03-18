library(Seurat)
library(pheatmap)
library(AUCell)
library(GSEABase)


MSC=readRDS('/data/MSC_no_C8_C567re-Order.rds') 

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


####PseudoBulk Differential analysis of AUC####
AUC=as.data.frame(cells_AUC@assays@data@listData$AUC)
library(DESeq2)

AUC_pseudobulk=cbind(rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='AD1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='AD2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='AD3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='BM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='BM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='BM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='PM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='PM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='PM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='UC1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='UC2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C1'&group=='UC3'))]),
                     
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='AD1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='AD2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='AD3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='BM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='BM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='BM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='PM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='PM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='PM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='UC1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='UC2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C2'&group=='UC3'))]),
                     
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='AD1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='AD2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='AD3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='BM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='BM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='BM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='PM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='PM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='PM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='UC1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='UC2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C3'&group=='UC3'))]),
                     
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='AD1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='AD2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='AD3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='BM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='BM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='BM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='PM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='PM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='PM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='UC1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='UC2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C4'&group=='UC3'))]),
                     
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='AD1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='AD2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='AD3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='BM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='BM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='BM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='PM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='PM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='PM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='UC1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='UC2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C5'&group=='UC3'))]),
                     
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='AD1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='AD2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='AD3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='BM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='BM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='BM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='PM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='PM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='PM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='UC1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='UC2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C6'&group=='UC3'))]),
                     
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='AD1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='AD2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='AD3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='BM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='BM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='BM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='PM1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='PM2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='PM3'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='UC1'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='UC2'))]),
                     rowSums(AUC[,row.names(subset(meta,Cluster=='C7'&group=='UC3'))]))

colnames(AUC_pseudobulk)=c('C1_AD1','C1_AD2','C1_AD3','C1_BM1','C1_BM2','C1_BM3','C1_PM1','C1_PM2','C1_PM3','C1_UC1','C1_UC2','C1_UC3',
                           'C2_AD1','C2_AD2','C2_AD3','C2_BM1','C2_BM2','C2_BM3','C2_PM1','C2_PM2','C2_PM3','C2_UC1','C2_UC2','C2_UC3',
                           'C3_AD1','C3_AD2','C3_AD3','C3_BM1','C3_BM2','C3_BM3','C3_PM1','C3_PM2','C3_PM3','C3_UC1','C3_UC2','C3_UC3',
                           'C4_AD1','C4_AD2','C4_AD3','C4_BM1','C4_BM2','C4_BM3','C4_PM1','C4_PM2','C4_PM3','C4_UC1','C4_UC2','C4_UC3',
                           'C5_AD1','C5_AD2','C5_AD3','C5_BM1','C5_BM2','C5_BM3','C5_PM1','C5_PM2','C5_PM3','C5_UC1','C5_UC2','C5_UC3',
                           'C6_AD1','C6_AD2','C6_AD3','C6_BM1','C6_BM2','C6_BM3','C6_PM1','C6_PM2','C6_PM3','C6_UC1','C6_UC2','C6_UC3',
                           'C7_AD1','C7_AD2','C7_AD3','C7_BM1','C7_BM2','C7_BM3','C7_PM1','C7_PM2','C7_PM3','C7_UC1','C7_UC2','C7_UC3')

design=read.csv('/data/NC_Revision_AUCScore/Design_cluster_AUC.csv',row.names = 'symbol')  ##design file for DEseq2; containing batch information

dds_Clusters=DESeqDataSetFromMatrix(countData = as.matrix (round(AUC_pseudobulk)),
                                    colData = as.matrix(design),
                                    design= ~ batch + C1)
##C1##
dds_C1=DESeq(dds_Clusters)
res_C1=as.matrix(results(dds_C1))

##C2##
design(dds_Clusters) <- formula(~batch + C2)
dds_C2 <- DESeq(dds_Clusters)
res_C2 <- as.matrix(results(dds_C2))

##C3##
design(dds_Clusters) <- formula(~batch + C3)
dds_C3 <- DESeq(dds_Clusters)
res_C3 <- as.matrix(results(dds_C3))

##C4##
design(dds_Clusters) <- formula(~batch + C4)
dds_C4 <- DESeq(dds_Clusters)
res_C4 <- as.matrix(results(dds_C4))

##C5##
design(dds_Clusters) <- formula(~batch + C5)
dds_C5 <- DESeq(dds_Clusters)
res_C5 <- as.matrix(results(dds_C5))

##C6##
design(dds_Clusters) <- formula(~batch + C6)
dds_C6 <- DESeq(dds_Clusters)
res_C6 <- as.matrix(results(dds_C6))

##C7##
design(dds_Clusters) <- formula(~batch + C7)
dds_C7 <- DESeq(dds_Clusters)
res_C7 <- as.matrix(results(dds_C7))
