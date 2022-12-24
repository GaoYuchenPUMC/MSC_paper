library(Seurat)
library(ggplot2)

##Load the reference dataset, 12 MSC sample as show in Fig. 1
MSC=readRDS('/data/12sample_MSC/MSC_no_C8_C567re-Order.rds')
DefaultAssay(MSC)='integrated'

####Seurat v4 Reference Mapping: FindTransferAnchor; Website:https://satijalab.org/seurat/articles/multimodal_reference_mapping.html####
B469=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/RDS files/B469.rds')
B537=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/RDS files/B537.rds')
B593=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/RDS files/B593.rds')
B801=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/RDS files/B801.rds')
B842=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/RDS files/B842.rds')
B844=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/RDS files/B844.rds')
B857=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/RDS files/B857.rds')
B927=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/RDS files/B927.rds')

AllBM=list(B469,B537,B593,B801,B842,B844,B857,B927)

MSC <- RunPCA(MSC, npcs = 30, verbose = FALSE,features = VariableFeatures(object = MSC))
MSC <- RunUMAP(MSC, reduction = "pca", dims = 1:30,return.model=TRUE,n.neighbors = 12L)
MSC <- FindNeighbors(object = MSC, reduction = "pca",dims = 1:30, 
                     k.param = 30, cache.index = TRUE, return.neighbor = TRUE, l2.norm = TRUE)

anchors <- list()
for (i in 1:length(AllBM)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = MSC,
    query = AllBM[[i]],
    k.filter = NA,
    reference.reduction = "pca", 
    dims = 1:30
  )
}

for (i in 1:length(AllBM)) {
  AllBM[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = AllBM[[i]],
    reference = MSC, 
    refdata = list(
      Cluster='Cluster'),
    reference.reduction = "pca",
    reduction.model = "umap"
  )
}

# Merge the batches 
AllBM <- merge(AllBM[[1]], AllBM[2:length(AllBM)], merge.dr = "ref.umap")
AllBM =NormalizeData(AllBM, normalization.method = "LogNormalize", scale.factor = 10000)

# Cluster with a total number of less than 100 were removed
# Then Name the new Cluster in BM-MSC dataset
meta=AllBM@meta.data

meta$predicted.Cluster[meta$predicted.Cluster=='C1']='BM1'
meta$predicted.Cluster[meta$predicted.Cluster=='C2']='BM2'
meta$predicted.Cluster[meta$predicted.Cluster=='C3']='BM3'
meta$predicted.Cluster[meta$predicted.Cluster=='C4']='BM4'
meta$predicted.Cluster[meta$predicted.Cluster=='C5']='BM5'
meta$predicted.Cluster[meta$predicted.Cluster=='C7']='BM6'

AllBM@meta.data=meta

####'Pseudobulk' strategy for differential expression testing####
pseudobulk_AllBM=cbind(rowSums(as.matrix(subset(AllBM,Cluster=='BM1'&group=='B469')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM1'&group=='B537')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM1'&group=='B593')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM1'&group=='B844')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM1'&group=='B801')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM1'&group=='B842')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM1'&group=='B857')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM1'&group=='B927')@assays$RNA@counts)),
                       
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM2'&group=='B469')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM2'&group=='B537')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM2'&group=='B593')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM2'&group=='B844')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM2'&group=='B801')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM2'&group=='B842')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM2'&group=='B927')@assays$RNA@counts)),
                       
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM3'&group=='B469')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM3'&group=='B537')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM3'&group=='B593')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM3'&group=='B844')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM3'&group=='B801')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM3'&group=='B842')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM3'&group=='B857')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM3'&group=='B927')@assays$RNA@counts)),
                       
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM4'&group=='B469')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM4'&group=='B537')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM4'&group=='B593')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM4'&group=='B844')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM4'&group=='B801')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM4'&group=='B842')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM4'&group=='B857')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM4'&group=='B927')@assays$RNA@counts)),
                       
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM5'&group=='B469')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM5'&group=='B537')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM5'&group=='B593')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM5'&group=='B844')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM5'&group=='B801')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM5'&group=='B842')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM5'&group=='B857')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM5'&group=='B927')@assays$RNA@counts)),
                       
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM6'&group=='B469')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM6'&group=='B537')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM6'&group=='B593')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM6'&group=='B844')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM6'&group=='B801')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM6'&group=='B842')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM6'&group=='B857')@assays$RNA@counts)),
                       rowSums(as.matrix(subset(AllBM,Cluster=='BM6'&group=='B927')@assays$RNA@counts)))


##B857-Aged3 dont have BM2
colnames(pseudobulk_AllBM)=c('BM1_Y1','BM1_Y2','BM1_Y3','BM1_Y4','BM1_A1','BM1_A2','BM1_A3','BM1_A4',
                             'BM2_Y1','BM2_Y2','BM2_Y3','BM2_Y4','BM2_A1','BM2_A2','BM2_A4',
                             'BM3_Y1','BM3_Y2','BM3_Y3','BM3_Y4','BM3_A1','BM3_A2','BM3_A3','BM3_A4',
                             'BM4_Y1','BM4_Y2','BM4_Y3','BM4_Y4','BM4_A1','BM4_A2','BM4_A3','BM4_A4',
                             'BM5_Y1','BM5_Y2','BM5_Y3','BM5_Y4','BM5_A1','BM5_A2','BM5_A3','BM5_A4',
                             'BM6_Y1','BM6_Y2','BM6_Y3','BM6_Y4','BM6_A1','BM6_A2','BM6_A3','BM6_A4')
write.csv(pseudobulk_AllBM,'C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/PseudoBulk DEG/Psudobulk expression YoungAgeBM.csv')

design=read.csv('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/PseudoBulk DEG/Design_cluster.csv')

library(DESeq2)
dds_Clusters=DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_AllBM),
                                    colData = as.matrix(design),
                                    design= ~ BM1)
##BM1##
dds_BM1=DESeq(dds_Clusters)
res_BM1=as.matrix(results(dds_BM1,alpha = 0.05))

##BM2##
dds_Clusters=DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_AllBM),
                                    colData = as.matrix(design),
                                    design= ~ BM2)
dds_BM2 <- DESeq(dds_Clusters)
res_BM2 <- as.matrix(results(dds_BM2))

##BM3##
dds_Clusters=DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_AllBM),
                                    colData = as.matrix(design),
                                    design= ~ BM3)
dds_BM3 <- DESeq(dds_Clusters)
res_BM3 <- as.matrix(results(dds_BM3))

##BM4##
dds_Clusters=DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_AllBM),
                                    colData = as.matrix(design),
                                    design= ~ BM4)
dds_BM4 <- DESeq(dds_Clusters)
res_BM4 <- as.matrix(results(dds_BM4))

##BM5##
dds_Clusters=DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_AllBM),
                                    colData = as.matrix(design),
                                    design= ~ BM5)
dds_BM5 <- DESeq(dds_Clusters)
res_BM5 <- as.matrix(results(dds_BM5))

##BM6##
dds_Clusters=DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_AllBM),
                                    colData = as.matrix(design),
                                    design= ~ BM6)
dds_BM6 <- DESeq(dds_Clusters)
res_BM6 <- as.matrix(results(dds_BM6))

####GO dotplot####
GO=read.csv('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/GO analysis/Dotplot GO logFC025.csv')
GO$Padjust=p.adjust(GO$Pvalue,method = 'BH')
GO$GO.term=factor(GO$GO.term,levels = c(rev(GO$GO.term)))
ggplot(GO,aes(Cluster,GO.term))+ geom_point(size=5,aes(color=Padjust))+ scale_colour_gradient(low="red",high="#F3A319")

write.csv(GO,'C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/GO analysis/Dotplot GO logFC025.csv')

####Gene dotplot####
AllBM@meta.data$Cluster=factor(AllBM@meta.data$Cluster,levels = c('BM1','BM2','BM3','BM4','BM5','BM6'))
DotPlot(AllBM,features = c('DNMT1','SUZ12','E2F1',
                           'EZH2','TOP2A','HMGB2',
                           'LMNB1','HDAC2','UBE2S',
                           'HSPA5','PDIA5','IL6','GLB1',
                           'IGFBP2','CDKN2A','COL6A1','IGFBP3'),cols = c("lightgrey", "red"),col.max = 1,dot.scale = 6,group.by = 'Cluster',scale.by = 'size' )+ RotatedAxis()

