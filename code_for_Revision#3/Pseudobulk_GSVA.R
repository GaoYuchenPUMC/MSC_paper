library(Seurat)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(DESeq2)

#### 'Pseudobulk' expression for each cluster from every sample ######

MSC=readRDS('/data/MSC_no_C8_C567re-Order.rds') ## RDS file 12 MSC sample from different tissue as shown in Fig. 1

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

#### limma remove batch effect ####
design=read.csv('/data/NC_Revision3_GSVA/Design_cluster_GSVA.csv',row.names = 'symbol')
pseudobulk_MSC=removeBatchEffect(pseudobulk_MSC,batch = design$batch)

#### GSVA analysis ####
geneset <- getGmt("/data/NC_Revision3_GSVA/GSVA_GeneSets.gmt")

gsva_mat_batchcorrect <- gsva(expr=pseudobulk_MSC, 
                 gset.idx.list=geneset, 
                 kcdf="Poisson" ,    #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T)

#### Pseudobulk differential analysis using limma ####
library(dplyr)
list <- c(rep("C1", 12), rep("C2", 12),rep("C3", 12),rep("C4", 12),
          rep("C5", 12),rep("C6", 12),rep("C7", 12)) %>% factor(., levels = c('C1','C2','C3','C4','C5','C6','C7'), ordered = F)
list <- model.matrix(~factor(list)+0)
colnames(list) <- c('C1','C2','C3','C4','C5','C6','C7')
df.fit <- lmFit(gsva_mat_batchcorrect, list)

##C1
df.matrix <- makeContrasts(contrasts = "C1", levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
C1_Output <- topTable(fit,n = Inf, adjust = "fdr")
##C2
df.matrix <- makeContrasts(contrasts = "C2", levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
C2_Output <- topTable(fit,n = Inf, adjust = "fdr")
##C3
df.matrix <- makeContrasts(contrasts = "C3", levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
C3_Output <- topTable(fit,n = Inf, adjust = "fdr")
##C4
df.matrix <- makeContrasts(contrasts = "C4", levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
C4_Output <- topTable(fit,n = Inf, adjust = "fdr")
##C5
df.matrix <- makeContrasts(contrasts = "C5", levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
C5_Output <- topTable(fit,n = Inf, adjust = "fdr")
##C6
df.matrix <- makeContrasts(contrasts = "C6", levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
C6_Output <- topTable(fit,n = Inf, adjust = "fdr")
##C7
df.matrix <- makeContrasts(contrasts = "C7", levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
C7_Output <- topTable(fit,n = Inf, adjust = "fdr")


#### Heatmap showing differential pathway (calculated by GSVA) ####
gsva_plot=cbind(rowMeans(gsva_mat_batchcorrect[,1:12]),rowMeans(gsva_mat_batchcorrect[,13:24]),rowMeans(gsva_mat_batchcorrect[,25:36]),
                rowMeans(gsva_mat_batchcorrect[,37:48]),rowMeans(gsva_mat_batchcorrect[,49:60]),rowMeans(gsva_mat_batchcorrect[,61:72]),
                rowMeans(gsva_mat_batchcorrect[,73:84]))
pheatmap(gsva_plot,cluster_rows = F,cluster_cols = F,scale = 'row',
         colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),
         border=FALSE,cellwidth = 10, cellheight = 10,breaks=unique(c(seq(-1.8,1.7, length=100))))
