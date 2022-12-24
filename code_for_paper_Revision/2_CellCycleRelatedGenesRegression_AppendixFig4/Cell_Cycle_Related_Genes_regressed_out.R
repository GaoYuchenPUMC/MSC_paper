#### Proliferation or Senescence: 1-Remove cellcycle related genes and redo dimensionality reduction analysis ####

AD1=readRDS('/data/gaoyuchen/10X_MSC/AD1.rds')
AD2=readRDS('/data/gaoyuchen/10X_MSC/AD2.rds')
AD3=readRDS('/data/gaoyuchen/10X_MSC/AD3.rds')

BM1=readRDS('/data/gaoyuchen/10X_MSC/BM1.rds')
BM2=readRDS('/data/gaoyuchen/10X_MSC/BM2.rds')
BM3=readRDS('/data/gaoyuchen/10X_MSC/BM3.rds')

PM1=readRDS('/data/gaoyuchen/10X_MSC/PM1.rds')
PM2=readRDS('/data/gaoyuchen/10X_MSC/PM2.rds')
PM3=readRDS('/data/gaoyuchen/10X_MSC/PM3.rds')

UC1=readRDS('/data/gaoyuchen/10X_MSC/UC1.rds')
UC2=readRDS('/data/gaoyuchen/10X_MSC/UC2.rds')
UC3=readRDS('/data/gaoyuchen/10X_MSC/UC3.rds')

Allsample=list(AD1,AD2,AD3,BM1,BM2,BM3,PM1,PM2,PM3,UC1,UC2,UC3)
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Allsample)
##Blacklist remove from features
s.genes <- (cc.genes$s.genes)
g2m.genes <- (cc.genes$g2m.genes)
ProlifGene=intersect(c(s.genes,g2m.genes),features)

features=setdiff(features,c(ProlifGene))

##FindAchor with new feature geneset that remove Cell cycle related genes
MSC.anchors <- FindIntegrationAnchors(object.list = Allsample, anchor.features = features)
MSC <- IntegrateData(anchorset = MSC.anchors)

##Dimensional reduction with new features
DefaultAssay(MSC) <- "integrated"
MSC <- ScaleData(MSC, vars.to.regress =c( "percent.mt",'nCount_RNA'))
MSC <- RunPCA(MSC, npcs = 30, verbose = FALSE,features = VariableFeatures(object = MSC))
ElbowPlot(MSC,ndims = 30)

MSC <- RunUMAP(MSC, reduction = "pca", dims = 1:30)
MSC<- FindNeighbors(MSC, reduction = "pca", dims = 1:30)
MSC<- FindClusters(MSC, resolution =0.4)
DimPlot(MSC, reduction = "umap", label = TRUE, pt.size = 0.5)

##UMAP project new and old##
DimPlot(MSC,label = T,pt.size = 0.6)+scale_color_manual(values = c('#009251',"#7CAE00","#72BE9C",'#5088BF',"#FF5BC6",'orange2','#EC494E'))

meta=MSC@meta.data
meta_old=MSC_old@meta.data
meta_old=meta_old[rownames(meta),]
meta$Cluster_old=meta_old$Cluster
MSC@meta.data=meta

DimPlot(MSC,label = T,pt.size = 0.6,group.by = 'Cluster_old')+scale_color_manual(values = c('#009251',"#7CAE00","#72BE9C",'#5088BF',"#FF5BC6",'orange2','#EC494E'))


##Fractions##
100*length(row.names(subset(meta,group=='AD1'&seurat_clusters==3)))/length(row.names(subset(meta,group=='AD1')))
100*length(row.names(subset(meta,group=='AD2'&seurat_clusters==3)))/length(row.names(subset(meta,group=='AD2')))
100*length(row.names(subset(meta,group=='AD3'&seurat_clusters==3)))/length(row.names(subset(meta,group=='AD3')))
100*length(row.names(subset(meta,group=='BM1'&seurat_clusters==3)))/length(row.names(subset(meta,group=='BM1')))
100*length(row.names(subset(meta,group=='BM2'&seurat_clusters==3)))/length(row.names(subset(meta,group=='BM2')))
100*length(row.names(subset(meta,group=='BM3'&seurat_clusters==3)))/length(row.names(subset(meta,group=='BM3')))

100*length(row.names(subset(meta,group=='PM1'&seurat_clusters==3)))/length(row.names(subset(meta,group=='PM1')))
100*length(row.names(subset(meta,group=='PM2'&seurat_clusters==3)))/length(row.names(subset(meta,group=='PM2')))
100*length(row.names(subset(meta,group=='PM3'&seurat_clusters==3)))/length(row.names(subset(meta,group=='PM3')))
100*length(row.names(subset(meta,group=='UC1'&seurat_clusters==3)))/length(row.names(subset(meta,group=='UC1')))
100*length(row.names(subset(meta,group=='UC2'&seurat_clusters==3)))/length(row.names(subset(meta,group=='UC2')))
100*length(row.names(subset(meta,group=='UC3'&seurat_clusters==3)))/length(row.names(subset(meta,group=='UC3')))

##Compare the barcode in each clusters of MSC/MSC remove cell cycle geenes from PCAanalysis##
MSC_old=readRDS('/data/gaoyuchen/MSC_no_C8_C567re-Order.rds')

C1_same=intersect(colnames(subset(MSC,idents='3')),colnames(subset(MSC_old,idents='C1')))
100*length(C1_same)/length(colnames(subset(MSC_old,idents='C1')))

C2_same=intersect(colnames(subset(MSC,idents='6')),colnames(subset(MSC_old,idents='C2')))
100*length(C2_same)/length(colnames(subset(MSC_old,idents='C2')))

C3_same=intersect(colnames(subset(MSC,idents='1')),colnames(subset(MSC_old,idents='C3')))
100*length(C3_same)/length(colnames(subset(MSC_old,idents='C3')))

C4_same=intersect(colnames(subset(MSC,idents='2')),colnames(subset(MSC_old,idents='C4')))
100*length(C4_same)/length(colnames(subset(MSC_old,idents='C4')))

C5_same=intersect(colnames(subset(MSC,idents='4')),colnames(subset(MSC_old,idents='C5')))
100*length(C5_same)/length(colnames(subset(MSC_old,idents='C5')))

C6_same=intersect(colnames(subset(MSC,idents='5')),colnames(subset(MSC_old,idents='C6')))
100*length(C6_same)/length(colnames(subset(MSC_old,idents='C6')))
 
C7_same=intersect(colnames(subset(MSC,idents='0')),colnames(subset(MSC_old,idents='C7')))
100*length(C7_same)/length(colnames(subset(MSC_old,idents='C7')))

##Pearson correlation analysis###

cov_matrix=cbind(AverageExpression(subset(MSC,idents='C1'))$RNA,AverageExpression(subset(MSC,idents='C2'))$RNA,AverageExpression(subset(MSC,idents='C3'))$RNA,
                 AverageExpression(subset(MSC,idents='C4'))$RNA,AverageExpression(subset(MSC,idents='C5'))$RNA,AverageExpression(subset(MSC,idents='C6'))$RNA,
                 AverageExpression(subset(MSC,idents='C7'))$RNA,AverageExpression(subset(MSC_old,idents='C1'))$RNA,AverageExpression(subset(MSC_old,idents='C2'))$RNA,
                 AverageExpression(subset(MSC_old,idents='C3'))$RNA,AverageExpression(subset(MSC_old,idents='C4'))$RNA,AverageExpression(subset(MSC_old,idents='C5'))$RNA,
                 AverageExpression(subset(MSC_old,idents='C6'))$RNA,AverageExpression(subset(MSC_old,idents='C7'))$RNA)

colnames(cov_matrix)=c('C1_new','C2_new','C3_new','C4_new','C5_new','C6_new','C7_new',
                       'C1','C2','C3','C4','C5','C6','C7')
cov_pearson <- cor(cov_matrix, method = 'pearson')
write.csv(cov_pearson,'/data/gaoyuchen/NC_Revision/Reviewer2_CellCycle_Parts/PCA_remove/correlation_pearson.csv')

pheatmap(cov_pearson[c(8:14),c(1:7)],cluster_rows = F,cluster_cols = F,colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),border=FALSE,cellwidth = 10, cellheight = 10,breaks=unique(c(seq(0.95,1, length=100))))

##Pseudobulk' strategy for DEG in PCA remove---use code above###

##Dotplot##
DotPlot(MSC,features = c('EZH2','DNMT1',
                         'LMNB1','MKI67',
                         'UBE2S','PTTG1',
                         'HSPB1','HSPA9','GLB1','PDIA6',
                         'TP53','B4GALT1','CTNNB1',
                         'CDKN1A','COL5A2','IGFBP5','COL12A1'),cols = c("lightgrey", "red"),col.max = 1.5,scale.min = )+ RotatedAxis()+coord_flip()

####2-Regress Out Cell cycle####
MSC=readRDS('/data/gaoyuchen/MSC_no_C8_C567re-Order.rds')
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(MSC)='RNA'
MSC <- CellCycleScoring(MSC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(MSC)='integrated'
MSC <- ScaleData(MSC, vars.to.regress = c("S.Score", "G2M.Score","percent.mt",'nCount_RNA','group'))

MSC <- RunPCA(MSC, npcs = 30, verbose = FALSE,features = VariableFeatures(object = MSC))

MSC <- RunUMAP(MSC, reduction = "pca", dims = 1:25)
MSC<- FindNeighbors(MSC, reduction = "pca", dims = 1:25)
MSC<- FindClusters(MSC, resolution =0.35)
DimPlot(MSC, reduction = "umap", label = TRUE, pt.size = 0.5)

##new meta and old meta comparing
MSC_old=readRDS('/data/gaoyuchen/MSC_no_C8_C567re-Order.rds')

meta=read.csv('/data/gaoyuchen/NC_Revision/Reviewer2_CellCycle_Parts/Regress_cellcycle/META.csv',row.names = 'barcode')
meta_old=MSC_old@meta.data
meta$Cluster_old=meta_old$Cluster
MSC@meta.data=meta

Idents(MSC)='Cluster'
DimPlot(MSC,label = T,pt.size = 0.6)+scale_color_manual(values = c('#5088BF','orange2',"#FF5BC6","#EC494E",'#009251',"#7CAE00"))

Idents(MSC)='Cluster_old'
DimPlot(MSC,label = T,pt.size = 0.6)+scale_color_manual(values = c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E"))

##Fractions##
100*length(row.names(subset(meta,group=='AD1'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='AD1')))
100*length(row.names(subset(meta,group=='AD2'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='AD2')))
100*length(row.names(subset(meta,group=='AD3'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='AD3')))
100*length(row.names(subset(meta,group=='BM1'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='BM1')))
100*length(row.names(subset(meta,group=='BM2'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='BM2')))
100*length(row.names(subset(meta,group=='BM3'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='BM3')))

100*length(row.names(subset(meta,group=='PM1'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='PM1')))
100*length(row.names(subset(meta,group=='PM2'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='PM2')))
100*length(row.names(subset(meta,group=='PM3'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='PM3')))
100*length(row.names(subset(meta,group=='UC1'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='UC1')))
100*length(row.names(subset(meta,group=='UC2'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='UC2')))
100*length(row.names(subset(meta,group=='UC3'&seurat_clusters=='0')))/length(row.names(subset(meta,group=='UC3')))

##Pearson correlation analysis###

MSC_old=readRDS('/data/gaoyuchen/MSC_no_C8_C567re-Order.rds')
cov_matrix=cbind(AverageExpression(subset(MSC,idents='2'))$RNA,AverageExpression(subset(MSC,idents='3'))$RNA,AverageExpression(subset(MSC,idents='1'))$RNA,
                 AverageExpression(subset(MSC,idents='5'))$RNA,AverageExpression(subset(MSC,idents='4'))$RNA,AverageExpression(subset(MSC,idents='0'))$RNA,
                 AverageExpression(subset(MSC_old,idents='C1'))$RNA,AverageExpression(subset(MSC_old,idents='C2'))$RNA,
                 AverageExpression(subset(MSC_old,idents='C3'))$RNA,AverageExpression(subset(MSC_old,idents='C4'))$RNA,AverageExpression(subset(MSC_old,idents='C5'))$RNA,
                 AverageExpression(subset(MSC_old,idents='C6'))$RNA,AverageExpression(subset(MSC_old,idents='C7'))$RNA)

colnames(cov_matrix)=c('non-Senescence 1','non-Senescence 2','C4_new','C5_new','C6_new','C7_new',
                       'C1','C2','C3','C4','C5','C6','C7')
cov_pearson <- cor(cov_matrix, method = 'pearson')
write.csv(cov_pearson,'/data/gaoyuchen/NC_Revision/Reviewer2_CellCycle_Parts/Regress_cellcycle/correlation_pearson.csv')

pheatmap(cov_pearson[c(7:13),c(1:6)],cluster_rows = F,cluster_cols = F,colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),border=FALSE,cellwidth = 10, cellheight = 10,breaks=unique(c(seq(0.95,1, length=100))))

