library(Seurat)
library(ggplot2)
library(DESeq2)

####'Pseudobulk' strategy for differential expression testing (considering batch effect)####

MSC=readRDS('/data/MSC_no_C8_C567re-Order.rds') ## RDS file: 12 MSC sample from different tissue as shown in Fig. 1

pseudobulk_MSC=cbind(rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='AD1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='AD2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='AD3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='BM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='BM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='BM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='PM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='PM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='PM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='UC1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='UC2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C1'&group=='UC3')@assays$RNA@counts)),
                     
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='AD1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='AD2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='AD3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='BM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='BM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='BM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='PM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='PM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='PM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='UC1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='UC2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C2'&group=='UC3')@assays$RNA@counts)),
                     
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='AD1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='AD2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='AD3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='BM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='BM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='BM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='PM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='PM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='PM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='UC1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='UC2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C3'&group=='UC3')@assays$RNA@counts)),
                     
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='AD1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='AD2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='AD3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='BM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='BM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='BM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='PM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='PM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='PM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='UC1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='UC2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C4'&group=='UC3')@assays$RNA@counts)),
                     
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='AD1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='AD2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='AD3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='BM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='BM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='BM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='PM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='PM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='PM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='UC1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='UC2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C5'&group=='UC3')@assays$RNA@counts)),
                     
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='AD1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='AD2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='AD3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='BM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='BM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='BM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='PM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='PM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='PM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='UC1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='UC2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C6'&group=='UC3')@assays$RNA@counts)),
                     
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='AD1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='AD2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='AD3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='BM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='BM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='BM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='PM1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='PM2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='PM3')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='UC1')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='UC2')@assays$RNA@counts)),
                     rowSums(as.matrix(subset(MSC,Cluster=='C7'&group=='UC3')@assays$RNA@counts)))

colnames(pseudobulk_MSC)=c('C1_AD1','C1_AD2','C1_AD3','C1_BM1','C1_BM2','C1_BM3','C1_PM1','C1_PM2','C1_PM3','C1_UC1','C1_UC2','C1_UC3',
                           'C2_AD1','C2_AD2','C2_AD3','C2_BM1','C2_BM2','C2_BM3','C2_PM1','C2_PM2','C2_PM3','C2_UC1','C2_UC2','C2_UC3',
                           'C3_AD1','C3_AD2','C3_AD3','C3_BM1','C3_BM2','C3_BM3','C3_PM1','C3_PM2','C3_PM3','C3_UC1','C3_UC2','C3_UC3',
                           'C4_AD1','C4_AD2','C4_AD3','C4_BM1','C4_BM2','C4_BM3','C4_PM1','C4_PM2','C4_PM3','C4_UC1','C4_UC2','C4_UC3',
                           'C5_AD1','C5_AD2','C5_AD3','C5_BM1','C5_BM2','C5_BM3','C5_PM1','C5_PM2','C5_PM3','C5_UC1','C5_UC2','C5_UC3',
                           'C6_AD1','C6_AD2','C6_AD3','C6_BM1','C6_BM2','C6_BM3','C6_PM1','C6_PM2','C6_PM3','C6_UC1','C6_UC2','C6_UC3',
                           'C7_AD1','C7_AD2','C7_AD3','C7_BM1','C7_BM2','C7_BM3','C7_PM1','C7_PM2','C7_PM3','C7_UC1','C7_UC2','C7_UC3')

design=read.csv('/data/1_DEG analysis/Design_cluster_revision2.csv',
                row.names = 'symbol')  ##design file for DEseq2 analysis (contain batch information)

dds_Clusters=DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_MSC),
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
