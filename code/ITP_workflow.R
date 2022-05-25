library('Seurat')
library('ggplot2')
##############Creat Seurat objects, performs QC filtering of the cells, and Data Norlmalization###############
######HD1--SampleID: NC2#####
NC2.data <- Read10X(data.dir = "/data/ITP_Tcells/N2PB/filtered_feature_bc_matrix/")
NC2 <- CreateSeuratObject(counts = NC2.data, project = "10X_NC2", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
NC2[["percent.mt"]] <- PercentageFeatureSet(NC2, pattern = "^MT-")
ggplot(NC2@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(NC2@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
NC2 <- subset(NC2, subset = nFeature_RNA > 250 & nFeature_RNA < 2900 & percent.mt < 11& nCount_RNA<12500)
##Normalize
NC2 <- NormalizeData(NC2, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
NC2 <- FindVariableFeatures(NC2, selection.method = "vst", nfeatures = 1300)

NC2$group<- "NC2"
NC2=RenameCells(NC2,add.cell.id = "NC2")

######HD2--SampleID: NC3#####
NC3.data <- Read10X(data.dir = "/data/ITP_Tcells/N3PB/filtered_feature_bc_matrix/")
NC3 <- CreateSeuratObject(counts = NC3.data, project = "10X_NC3", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
NC3[["percent.mt"]] <- PercentageFeatureSet(NC3, pattern = "^MT-")
ggplot(NC3@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(NC3@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
NC3 <- subset(NC3, subset = nFeature_RNA > 250 & nFeature_RNA < 2900 & percent.mt < 7 & nCount_RNA<12500)
##Normalize
NC3 <- NormalizeData(NC3, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
NC3 <- FindVariableFeatures(NC3, selection.method = "vst", nfeatures = 1300)

NC3$group<- "NC3"
NC3=RenameCells(NC3,add.cell.id = "NC3")

######HD3--SampleID: NC4#####
NC4.data <- Read10X(data.dir = "/data/ITP_Tcells/N4PB/filtered_feature_bc_matrix/")
NC4 <- CreateSeuratObject(counts = NC4.data, project = "10X_NC4", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
NC4[["percent.mt"]] <- PercentageFeatureSet(NC4, pattern = "^MT-")
ggplot(NC4@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(NC4@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
NC4 <- subset(NC4, subset = nFeature_RNA > 250 & nFeature_RNA < 3250 & percent.mt < 11 & nCount_RNA<12500)
##Normalize
NC4 <- NormalizeData(NC4, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
NC4 <- FindVariableFeatures(NC4, selection.method = "vst", nfeatures = 1300)

NC4$group<- "NC4"
NC4=RenameCells(NC4,add.cell.id = "NC4")

######P1--SampleID: P1A#####
P1A.data <- Read10X(data.dir = "/data/ITP_Tcells/P1A/filtered_feature_bc_matrix/")
P1A <- CreateSeuratObject(counts = P1A.data, project = "10X_P1A", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
P1A[["percent.mt"]] <- PercentageFeatureSet(P1A, pattern = "^MT-")
ggplot(P1A@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(P1A@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
P1A <- subset(P1A, subset = nFeature_RNA > 250 & nFeature_RNA < 2750 & percent.mt < 8.5 & nCount_RNA<10000)
##Normalize
P1A <- NormalizeData(P1A, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
P1A <- FindVariableFeatures(P1A, selection.method = "vst", nfeatures = 1300)

P1A$group<- "P1A"
P1A=RenameCells(P1A,add.cell.id = "P1A")

######P2--SampleID: P5A#####
P5A.data <- Read10X(data.dir = "/data/ITP_Tcells/P5A/filtered_feature_bc_matrix/")
P5A <- CreateSeuratObject(counts = P5A.data, project = "10X_P5A", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
P5A[["percent.mt"]] <- PercentageFeatureSet(P5A, pattern = "^MT-")
ggplot(P5A@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(P5A@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
P5A <- subset(P5A, subset = nFeature_RNA > 250 & nFeature_RNA < 2750 & percent.mt < 10 & nCount_RNA<10000)
##Normalize
P5A <- NormalizeData(P5A, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
P5A <- FindVariableFeatures(P5A, selection.method = "vst", nfeatures = 1300)

P5A$group<- "P5A"
P5A=RenameCells(P5A,add.cell.id = "P5A")

######P3--SampleID: P6A######
P6A.data <- Read10X(data.dir = "/data/ITP_Tcells/P6A/filtered_feature_bc_matrix/")
P6A <- CreateSeuratObject(counts = P6A.data, project = "10X_P6A", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
P6A[["percent.mt"]] <- PercentageFeatureSet(P6A, pattern = "^MT-")
ggplot(P6A@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(P6A@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
P6A <- subset(P6A, subset = nFeature_RNA > 250 & nFeature_RNA < 2750 & percent.mt < 7 & nCount_RNA<12500)
##Normalize
P6A <- NormalizeData(P6A, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
P6A <- FindVariableFeatures(P6A, selection.method = "vst", nfeatures = 1300)

P6A$group<- "P6A"
P6A=RenameCells(P6A,add.cell.id = "P6A")

######FindAnchor and Data Integration#######
Allsample=list(NC2,NC3,NC4,
               P1A,P5A,P6A)
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Allsample)
##Blacklist
RpGene=features[(grep(features, pattern = "^RP([0-9]+-|[LS])"))]
MtGene=features[(grep(features, pattern = "^MT"))]

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ProlifGene=intersect(c(s.genes,g2m.genes),features)

TCR=read.csv('/data/ITP_Tcells/TCRecptor.csv',header = T,sep=',')
TCRGene=intersect(c(as.vector(TCR$Approved.symbol)),features)

features=setdiff(features,c(RpGene,MtGene,ProlifGene,TCRGene))

PBMC.anchors <- FindIntegrationAnchors(object.list = Allsample, anchor.features = features)
PBMC <- IntegrateData(anchorset = PBMC.anchors)

##Integrate analysis
DefaultAssay(PBMC)='integrated'
PBMC <- ScaleData(PBMC, vars.to.regress = c("percent.mt", "nCount_RNA"))
PBMC <- RunPCA(PBMC, npcs = 30, verbose = FALSE)
ElbowPlot(PBMC,ndims=30)
PBMC <- RunUMAP(PBMC, reduction = "pca", dims = 1:25)
PBMC <- FindNeighbors(PBMC, reduction = "pca", dims = 1:25)
PBMC <- FindClusters(PBMC, resolution = 0.6)
DimPlot(PBMC, reduction = "umap", label = TRUE,pt.size = 0.5)

DefaultAssay(PBMC)='RNA'
markers=FindAllMarkers(PBMC,logfc.threshold = 0.3,only.pos = T)

saveRDS(PBMC,'/data/ITP_Tcells/ALL_PBMC_20220412.rds')

##remove doublets and Re-run the UMAP## 
PBMC=subset(PBMC,idents =c(2,11,5,18,28,16,13,19,10,3,15,7,0,29,1,6,24,21,34,4,17,26,30,8,9,12,20) )
DefaultAssay(PBMC)='integrated'
PBMC <- ScaleData(PBMC, vars.to.regress = c("percent.mt", "nCount_RNA"))
PBMC <- RunPCA(PBMC, npcs = 30, verbose = FALSE)
ElbowPlot(PBMC,ndims=30)
PBMC <- RunUMAP(PBMC, reduction = "pca", dims = 1:25)
PBMC <- FindNeighbors(PBMC, reduction = "pca", dims = 1:25)
PBMC <- FindClusters(PBMC, resolution = 0.6)
DimPlot(PBMC, reduction = "umap", label = TRUE,pt.size = 0.5)

####Cell Proportion: B cells####
meta=PBMC@meta.data
##niave B-- inside B population
Bcells=subset(meta,seurat_clusters==6|seurat_clusters==12|seurat_clusters==18)
100*length(row.names(subset(meta,seurat_clusters==6&group=='NC2')))/length(row.names(subset(Bcells,group=='NC2')))
100*length(row.names(subset(meta,seurat_clusters==6&group=='NC3')))/length(row.names(subset(Bcells,group=='NC3')))
100*length(row.names(subset(meta,seurat_clusters==6&group=='NC4')))/length(row.names(subset(Bcells,group=='NC4')))

100*length(row.names(subset(meta,seurat_clusters==6&group=='P1A')))/length(row.names(subset(Bcells,group=='P1A')))
100*length(row.names(subset(meta,seurat_clusters==6&group=='P5A')))/length(row.names(subset(Bcells,group=='P5A')))
100*length(row.names(subset(meta,seurat_clusters==6&group=='P6A')))/length(row.names(subset(Bcells,group=='P6A')))

##Memory B-- inside B population
100*length(row.names(subset(meta,seurat_clusters==12&group=='NC2')))/length(row.names(subset(Bcells,group=='NC2')))
100*length(row.names(subset(meta,seurat_clusters==12&group=='NC3')))/length(row.names(subset(Bcells,group=='NC3')))
100*length(row.names(subset(meta,seurat_clusters==12&group=='NC4')))/length(row.names(subset(Bcells,group=='NC4')))

100*length(row.names(subset(meta,seurat_clusters==12&group=='P1A')))/length(row.names(subset(Bcells,group=='P1A')))
100*length(row.names(subset(meta,seurat_clusters==12&group=='P5A')))/length(row.names(subset(Bcells,group=='P5A')))
100*length(row.names(subset(meta,seurat_clusters==12&group=='P6A')))/length(row.names(subset(Bcells,group=='P6A')))

##Plasma B-- inside B population
100*length(row.names(subset(meta,seurat_clusters==18&group=='NC2')))/length(row.names(subset(Bcells,group=='NC2')))
100*length(row.names(subset(meta,seurat_clusters==18&group=='NC3')))/length(row.names(subset(Bcells,group=='NC3')))
100*length(row.names(subset(meta,seurat_clusters==18&group=='NC4')))/length(row.names(subset(Bcells,group=='NC4')))

100*length(row.names(subset(meta,seurat_clusters==18&group=='P1A')))/length(row.names(subset(Bcells,group=='P1A')))
100*length(row.names(subset(meta,seurat_clusters==18&group=='P5A')))/length(row.names(subset(Bcells,group=='P5A')))
100*length(row.names(subset(meta,seurat_clusters==18&group=='P6A')))/length(row.names(subset(Bcells,group=='P6A')))

####Cell Proportion  T cells####
meta=PBMC@meta.data
##niave CD4 T-- inside T population
Tcells=subset(meta,seurat_clusters==3|seurat_clusters==10|seurat_clusters==14|seurat_clusters==1|seurat_clusters==9|seurat_clusters==5|
                seurat_clusters==19|seurat_clusters==8|seurat_clusters==0|seurat_clusters==11|seurat_clusters==20|seurat_clusters==4|
                seurat_clusters==7|seurat_clusters==16)
100*length(row.names(subset(meta,seurat_clusters==3&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #12.10575
100*length(row.names(subset(meta,seurat_clusters==3&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) # 19.3579
100*length(row.names(subset(meta,seurat_clusters==3&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #13.48563

100*length(row.names(subset(meta,seurat_clusters==3&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #11.26513
100*length(row.names(subset(meta,seurat_clusters==3&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #7.718201
100*length(row.names(subset(meta,seurat_clusters==3&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #6.629173

##niave CD8 T-- inside T population
100*length(row.names(subset(meta,seurat_clusters==10&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #4.499072
100*length(row.names(subset(meta,seurat_clusters==10&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) # 7.237842
100*length(row.names(subset(meta,seurat_clusters==10&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #12.28199

100*length(row.names(subset(meta,seurat_clusters==10&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #2.461984
100*length(row.names(subset(meta,seurat_clusters==10&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #5.241645
100*length(row.names(subset(meta,seurat_clusters==10&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #5.46847

##Memory CD4 T-- inside T population
100*length(row.names(subset(meta,seurat_clusters==1&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #19.57328
100*length(row.names(subset(meta,seurat_clusters==1&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) # 17.87614
100*length(row.names(subset(meta,seurat_clusters==1&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #20.8057

100*length(row.names(subset(meta,seurat_clusters==1&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #17.76146
100*length(row.names(subset(meta,seurat_clusters==1&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #13.44073
100*length(row.names(subset(meta,seurat_clusters==1&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #13.38997

##T reg-- inside T population
100*length(row.names(subset(meta,seurat_clusters==14&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #1.785714
100*length(row.names(subset(meta,seurat_clusters==14&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) # 2.070669
100*length(row.names(subset(meta,seurat_clusters==14&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #1.621223

100*length(row.names(subset(meta,seurat_clusters==14&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #2.482673
100*length(row.names(subset(meta,seurat_clusters==14&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #1.41861
100*length(row.names(subset(meta,seurat_clusters==14&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #1.59148

##Effect CD8-- inside T population
100*length(row.names(subset(meta,seurat_clusters==5&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #11.9782
100*length(row.names(subset(meta,seurat_clusters==5&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) # 5.775076
100*length(row.names(subset(meta,seurat_clusters==5&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #4.937362

100*length(row.names(subset(meta,seurat_clusters==5&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #11.68925
100*length(row.names(subset(meta,seurat_clusters==5&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #9.1849
100*length(row.names(subset(meta,seurat_clusters==5&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #8.639464

##NKT-- inside T population
100*length(row.names(subset(meta,seurat_clusters==9&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #3.049629
100*length(row.names(subset(meta,seurat_clusters==9&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #9.080547
100*length(row.names(subset(meta,seurat_clusters==9&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #10.19406

100*length(row.names(subset(meta,seurat_clusters==9&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #1.520637
100*length(row.names(subset(meta,seurat_clusters==9&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #6.01106
100*length(row.names(subset(meta,seurat_clusters==9&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #8.268517

##Cytoto CD4-- inside T population
100*length(row.names(subset(meta,seurat_clusters==8&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #4.012059
100*length(row.names(subset(meta,seurat_clusters==8&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #5.832067
100*length(row.names(subset(meta,seurat_clusters==8&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #3.193319

100*length(row.names(subset(meta,seurat_clusters==8&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #11.21341
100*length(row.names(subset(meta,seurat_clusters==8&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #5.289733
100*length(row.names(subset(meta,seurat_clusters==8&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #3.936819

##memoryCD8-- inside T population
100*length(row.names(subset(meta,seurat_clusters==11&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #4.34833
100*length(row.names(subset(meta,seurat_clusters==11&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #3.305471
100*length(row.names(subset(meta,seurat_clusters==11&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #2.579219

100*length(row.names(subset(meta,seurat_clusters==11&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #3.961932
100*length(row.names(subset(meta,seurat_clusters==11&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #2.21207
100*length(row.names(subset(meta,seurat_clusters==11&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #5.444538

##CytotoxicCD8-- inside T population
100*length(row.names(subset(meta,seurat_clusters==0&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #11.65353
100*length(row.names(subset(meta,seurat_clusters==0&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #6.591945
100*length(row.names(subset(meta,seurat_clusters==0&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #7.221813

100*length(row.names(subset(meta,seurat_clusters==0&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #28.08524
100*length(row.names(subset(meta,seurat_clusters==0&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #31.01707
100*length(row.names(subset(meta,seurat_clusters==0&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #28.47912

##KLRC2NK-- inside T population
100*length(row.names(subset(meta,seurat_clusters==4&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #15.96707
100*length(row.names(subset(meta,seurat_clusters==4&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #11.60714
100*length(row.names(subset(meta,seurat_clusters==4&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #10.31688

100*length(row.names(subset(meta,seurat_clusters==4&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) # 4.272266
100*length(row.names(subset(meta,seurat_clusters==4&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) # 14.23419
100*length(row.names(subset(meta,seurat_clusters==4&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #11.91815

##CD56-NK-- inside T population
100*length(row.names(subset(meta,seurat_clusters==7&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #9.450371
100*length(row.names(subset(meta,seurat_clusters==7&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #9.707447
100*length(row.names(subset(meta,seurat_clusters==7&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #11.86441

100*length(row.names(subset(meta,seurat_clusters==7&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) # 4.3964
100*length(row.names(subset(meta,seurat_clusters==7&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #3.34215
100*length(row.names(subset(meta,seurat_clusters==7&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #5.181285

##CD56+NK-- inside T population
100*length(row.names(subset(meta,seurat_clusters==16&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #1.124768
100*length(row.names(subset(meta,seurat_clusters==16&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #1.367781
100*length(row.names(subset(meta,seurat_clusters==16&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #1.154507

100*length(row.names(subset(meta,seurat_clusters==16&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) # 0.548257
100*length(row.names(subset(meta,seurat_clusters==16&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #0.4808848
100*length(row.names(subset(meta,seurat_clusters==16&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #0.6222329

##IFN-- inside T population
100*length(row.names(subset(meta,seurat_clusters==19&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #0.359462
100*length(row.names(subset(meta,seurat_clusters==19&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #0.1139818
100*length(row.names(subset(meta,seurat_clusters==19&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #0.2210759

100*length(row.names(subset(meta,seurat_clusters==19&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #0.2275784
100*length(row.names(subset(meta,seurat_clusters==19&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #0.2404424
100*length(row.names(subset(meta,seurat_clusters==19&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #0.3470145

##Divding-- inside T population
100*length(row.names(subset(meta,seurat_clusters==20&group=='NC2')))/length(row.names(subset(Tcells,group=='NC2'))) #0.09276438
100*length(row.names(subset(meta,seurat_clusters==20&group=='NC3')))/length(row.names(subset(Tcells,group=='NC3'))) #0.07598784
100*length(row.names(subset(meta,seurat_clusters==20&group=='NC4')))/length(row.names(subset(Tcells,group=='NC4'))) #0.1228199

100*length(row.names(subset(meta,seurat_clusters==20&group=='P1A')))/length(row.names(subset(Tcells,group=='P1A'))) #0.1137892
100*length(row.names(subset(meta,seurat_clusters==20&group=='P5A')))/length(row.names(subset(Tcells,group=='P5A'))) #0.1683097
100*length(row.names(subset(meta,seurat_clusters==20&group=='P6A')))/length(row.names(subset(Tcells,group=='P6A'))) #0.08376212

