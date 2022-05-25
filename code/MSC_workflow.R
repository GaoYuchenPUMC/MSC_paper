library('Seurat')
library('ggplot2')
##############Creat Seurat objects, performs QC filtering of the cells, and Data Norlmalization###############
######AD1#####
AD1.data <- Read10X(data.dir = "/data/AD1/filtered_feature_bc_matrix/")
AD1 <- CreateSeuratObject(counts = AD1.data, project = "10X_AD1", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
AD1[["percent.mt"]] <- PercentageFeatureSet(AD1, pattern = "^MT-")
ggplot(AD1@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(AD1@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
AD1 <- subset(AD1, subset = nFeature_RNA > 3600 & nFeature_RNA < 5600 & percent.mt < 3.8& nCount_RNA<55000& nCount_RNA>17000)
##Normalize
AD1 <- NormalizeData(AD1, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
AD1 <- FindVariableFeatures(AD1, selection.method = "vst", nfeatures = 1500)

AD1$group<- "AD1"
AD1=RenameCells(AD1,add.cell.id = "AD1")

######AD2#####
AD2.data <- Read10X(data.dir = "/data/AD2/filtered_feature_bc_matrix/")
AD2 <- CreateSeuratObject(counts = AD2.data, project = "10X_AD2", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
AD2[["percent.mt"]] <- PercentageFeatureSet(AD2, pattern = "^MT-")
ggplot(AD2@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(AD2@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
AD2 <- subset(AD2, subset = nFeature_RNA > 3600 & nFeature_RNA < 5750 & percent.mt < 4& nCount_RNA<45000& nCount_RNA>16000)
##Normalize
AD2 <- NormalizeData(AD2, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
AD2 <- FindVariableFeatures(AD2, selection.method = "vst", nfeatures = 1500)

AD2$group<- "AD2"
AD2=RenameCells(AD2,add.cell.id = "AD2")

######AD3#####
AD3.data <- Read10X(data.dir = "/data/AD3/filtered_feature_bc_matrix/")
AD3 <- CreateSeuratObject(counts = AD3.data, project = "10X_AD3", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
AD3[["percent.mt"]] <- PercentageFeatureSet(AD3, pattern = "^MT-")
ggplot(AD3@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(AD3@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
AD3 <- subset(AD3, subset = nFeature_RNA > 3250 & nFeature_RNA < 5750 & percent.mt < 3.8& nCount_RNA<50000 & nCount_RNA>15000)
##Normalize
AD3 <- NormalizeData(AD3, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
AD3 <- FindVariableFeatures(AD3, selection.method = "vst", nfeatures = 1500)

AD3$group<- "AD3"
AD3=RenameCells(AD3,add.cell.id = "AD3")

######BM1#####
BM1.data <- Read10X(data.dir = "/data/BM1/filtered_feature_bc_matrix/")
BM1 <- CreateSeuratObject(counts = BM1.data, project = "10X_BM1", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
BM1[["percent.mt"]] <- PercentageFeatureSet(BM1, pattern = "^MT-")
ggplot(BM1@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(BM1@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
BM1 <- subset(BM1, subset = nFeature_RNA > 3100 & nFeature_RNA < 5000 & percent.mt < 2.5 & nCount_RNA<40000& nCount_RNA>12500)
##Normalize
BM1 <- NormalizeData(BM1, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
BM1 <- FindVariableFeatures(BM1, selection.method = "vst", nfeatures = 1500)

BM1$group<- "BM1"
BM1=RenameCells(BM1,add.cell.id = "BM1")

######BM2#####
BM2.data <- Read10X(data.dir = "/data/BM2/filtered_feature_bc_matrix/")
BM2 <- CreateSeuratObject(counts = BM2.data, project = "10X_BM2", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
BM2[["percent.mt"]] <- PercentageFeatureSet(BM2, pattern = "^MT-")
ggplot(BM2@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(BM2@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
BM2 <- subset(BM2, subset = nFeature_RNA > 3800 & nFeature_RNA < 5500 & percent.mt < 4& nCount_RNA<37500& nCount_RNA>17500)
##Normalize
BM2 <- NormalizeData(BM2, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
BM2 <- FindVariableFeatures(BM2, selection.method = "vst", nfeatures = 1400)

BM2$group<- "BM2"
BM2=RenameCells(BM2,add.cell.id = "BM2")

######BM3#####
BM3.data <- Read10X(data.dir = "/data/BM3/filtered_feature_bc_matrix/")
BM3 <- CreateSeuratObject(counts = BM3.data, project = "10X_BM3", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
BM3[["percent.mt"]] <- PercentageFeatureSet(BM3, pattern = "^MT-")
ggplot(BM3@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(BM3@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
BM3 <- subset(BM3, subset = nFeature_RNA > 2800 & nFeature_RNA < 4100 & percent.mt < 4.5& nCount_RNA<20000& nCount_RNA>10000)
##Normalize
BM3 <- NormalizeData(BM3, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
BM3 <- FindVariableFeatures(BM3, selection.method = "vst", nfeatures = 1400)

BM3$group<- "BM3"
BM3=RenameCells(BM3,add.cell.id = "BM3")

######PM1#####
PM1.data <- Read10X(data.dir = "/data/PM1/filtered_feature_bc_matrix/")
PM1 <- CreateSeuratObject(counts = PM1.data, project = "10X_PM1", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
PM1[["percent.mt"]] <- PercentageFeatureSet(PM1, pattern = "^MT-")
ggplot(PM1@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(PM1@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
PM1 <- subset(PM1, subset = nFeature_RNA > 4300 & nFeature_RNA < 6300 & percent.mt < 4.5& nCount_RNA<55000& nCount_RNA>22500)
##Normalize
PM1 <- NormalizeData(PM1, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
PM1 <- FindVariableFeatures(PM1, selection.method = "vst", nfeatures = 1500)

PM1$group<- "PM1"
PM1=RenameCells(PM1,add.cell.id = "PM1")

######PM2#####
PM2.data <- Read10X(data.dir = "/data/PM2/filtered_feature_bc_matrix/")
PM2 <- CreateSeuratObject(counts = PM2.data, project = "10X_PM2", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
PM2[["percent.mt"]] <- PercentageFeatureSet(PM2, pattern = "^MT-")
ggplot(PM2@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(PM2@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
PM2 <- subset(PM2, subset = nFeature_RNA > 3500 & nFeature_RNA < 5600 & percent.mt < 4& nCount_RNA<47500& nCount_RNA>17500)
##Normalize
PM2 <- NormalizeData(PM2, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
PM2 <- FindVariableFeatures(PM2, selection.method = "vst", nfeatures = 1500)

PM2$group<- "PM2"
PM2=RenameCells(PM2,add.cell.id = "PM2")

######PM3#####
PM3.data <- Read10X(data.dir = "/data/PM3/filtered_feature_bc_matrix/")
PM3 <- CreateSeuratObject(counts = PM3.data, project = "10X_PM3", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
PM3[["percent.mt"]] <- PercentageFeatureSet(PM3, pattern = "^MT-")
ggplot(PM3@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(PM3@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
PM3 <- subset(PM3, subset = nFeature_RNA > 3750 & nFeature_RNA < 5800 & percent.mt < 3.7& nCount_RNA<50000& nCount_RNA>17500)
##Normalize
PM3 <- NormalizeData(PM3, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
PM3 <- FindVariableFeatures(PM3, selection.method = "vst", nfeatures = 1500)

PM3$group<- "PM3"
PM3=RenameCells(PM3,add.cell.id = "PM3")

######UC1#####
UC1.data <- Read10X(data.dir = "/data/UC1/filtered_feature_bc_matrix/")
UC1 <- CreateSeuratObject(counts = UC1.data, project = "10X_UC1", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
UC1[["percent.mt"]] <- PercentageFeatureSet(UC1, pattern = "^MT-")
ggplot(UC1@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(UC1@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
UC1 <- subset(UC1, subset = nFeature_RNA > 3800 & nFeature_RNA < 5800 & percent.mt < 4.2& nCount_RNA<53000& nCount_RNA>15000)
##Normalize
UC1 <- NormalizeData(UC1, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
UC1 <- FindVariableFeatures(UC1, selection.method = "vst", nfeatures = 1500)

UC1$group<- "UC1"
UC1=RenameCells(UC1,add.cell.id = "UC1")

######UC2#####
UC2.data <- Read10X(data.dir = "/data/UC2/filtered_feature_bc_matrix/")
UC2 <- CreateSeuratObject(counts = UC2.data, project = "10X_UC2", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
UC2[["percent.mt"]] <- PercentageFeatureSet(UC2, pattern = "^MT-")
ggplot(UC2@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(UC2@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
UC2 <- subset(UC2, subset = nFeature_RNA > 3750 & nFeature_RNA < 5850 & percent.mt < 3.24& nCount_RNA<50000& nCount_RNA>15000)
##Normalize
UC2 <- NormalizeData(UC2, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
UC2 <- FindVariableFeatures(UC2, selection.method = "vst", nfeatures = 1500)

UC2$group<- "UC2"
UC2=RenameCells(UC2,add.cell.id = "UC2")

######UC3#####
UC3.data <- Read10X(data.dir = "/data/UC3/filtered_feature_bc_matrix/")
UC3 <- CreateSeuratObject(counts = UC3.data, project = "10X_UC3", min.cells = 20, min.features = 200)
#We calculate the percentage of mitochondrial genes here and store it in percent.mt.
UC3[["percent.mt"]] <- PercentageFeatureSet(UC3, pattern = "^MT-")
ggplot(UC3@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(UC3@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))
UC3 <- subset(UC3, subset = nFeature_RNA > 3500 & nFeature_RNA < 5850 & percent.mt < 3& nCount_RNA<55000& nCount_RNA>15000)
##Normalize
UC3 <- NormalizeData(UC3, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
UC3 <- FindVariableFeatures(UC3, selection.method = "vst", nfeatures = 1500)

UC3$group<- "UC3"
UC3=RenameCells(UC3,add.cell.id = "UC3")

######FindAnchor and Data Integration#######
Allsample=list(AD1,AD2,AD3,BM1,BM2,BM3,
               PM1,PM2,PM3,UC1,UC2,UC3)
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Allsample)
MSC.anchors <- FindIntegrationAnchors(object.list = Allsample, anchor.features = features)
MSC <- IntegrateData(anchorset = MSC.anchors)

##Integrate analysis
DefaultAssay(MSC)='integrated'
MSC <- ScaleData(MSC, vars.to.regress = c("percent.mt", "nCount_RNA"))
MSC <- RunPCA(MSC, npcs = 30, verbose = FALSE)
ElbowPlot(MSC,ndims=30)
MSC <- RunUMAP(MSC, reduction = "pca", dims = 1:30,n.neighbors = 12L)
MSC <- FindNeighbors(MSC, reduction = "pca", dims = 1:30)
MSC <- FindClusters(MSC, resolution = 0.4)
DimPlot(MSC, reduction = "umap", label = TRUE,pt.size = 0.5)

DefaultAssay(MSC)='RNA'
markers=FindAllMarkers(MSC,logfc.threshold = 0.2,only.pos = T)

