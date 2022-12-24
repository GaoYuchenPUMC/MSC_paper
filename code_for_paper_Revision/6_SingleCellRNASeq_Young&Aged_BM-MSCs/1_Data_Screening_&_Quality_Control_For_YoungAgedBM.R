library(Seurat)
library(ggplot2)

####YoungBM1 B469####
B469.data <- Read10X(data.dir = "/data/NC_Revision_YoungAgeBM/CellRanger_Outs/Young_B469/outs/filtered_feature_bc_matrix")
B469 <- CreateSeuratObject(counts = B469.data, project = "10X_MSC", min.cells = 20, min.features = 200)
#We calculate the percentage of
# mitochondrial genes here and store it in percent.mt.
B469[["percent.mt"]] <- PercentageFeatureSet(B469, pattern = "^MT-")

###
ggplot(B469@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(B469@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))

##!!remove contamination
B469 <- subset(B469, subset = nFeature_RNA > 3000 & nFeature_RNA < 7000 & percent.mt < 7.5& nCount_RNA<55000 & nCount_RNA>10000)
##Normalization
B469 <- NormalizeData(B469, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
B469 <- FindVariableFeatures(B469, selection.method = "vst", nfeatures = 1500)
VariableFeaturePlot(B469)
##ScaleData
B469 <- ScaleData(B469, vars.to.regress =c( "percent.mt",'nCount_RNA'))
##
B469=RenameCells(B469,add.cell.id = "Young_B469")
B469$Type='Young'
B469$group='B469'
B469$Sample='YoungBM1'

##Dimensional Reduction
B469 <- RunPCA(B469, npcs = 30, verbose = FALSE,features = VariableFeatures(object = B469))
B469 <- RunUMAP(B469, reduction = "pca", dims = 1:20)

B469<- FindNeighbors(B469, reduction = "pca", dims = 1:20)
B469<- FindClusters(B469, resolution = 0.4)

DimPlot(B469, reduction = "umap", label = TRUE, pt.size = 1.3)

####YoungBM2 B537####
B537.data <- Read10X(data.dir = "/data/NC_Revision_YoungAgeBM/CellRanger_Outs/Young_B537/outs/filtered_feature_bc_matrix/")
B537 <- CreateSeuratObject(counts = B537.data, project = "10X_MSC", min.cells = 20, min.features = 200)
#We calculate the percentage of
# mitochondrial genes here and store it in percent.mt.
B537[["percent.mt"]] <- PercentageFeatureSet(B537, pattern = "^MT-")

###
ggplot(B537@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(B537@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))

##!!remove contamination
B537 <- subset(B537, subset = nFeature_RNA > 2250 & nFeature_RNA < 5750 & percent.mt < 6 & nCount_RNA<35000 & nCount_RNA>5000)
##Normalization
B537 <- NormalizeData(B537, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
B537 <- FindVariableFeatures(B537, selection.method = "vst", nfeatures = 1500)
VariableFeaturePlot(B537)
##ScaleData
B537 <- ScaleData(B537, vars.to.regress =c( "percent.mt",'nCount_RNA'))
##
B537=RenameCells(B537,add.cell.id = "Young_B537")
B537$Type='Young'
B537$group='B537'
B537$Sample='YoungBM2'

##Dimensional Reduction
B537 <- RunPCA(B537, npcs = 30, verbose = FALSE,features = VariableFeatures(object = B537))
B537 <- RunUMAP(B537, reduction = "pca", dims = 1:20)

B537<- FindNeighbors(B537, reduction = "pca", dims = 1:20)
B537<- FindClusters(B537, resolution = 0.4)

DimPlot(B537, reduction = "umap", label = TRUE, pt.size = 1.3)
####YoungBM3 B593####
B593.data <- Read10X(data.dir = "/data/NC_Revision_YoungAgeBM/CellRanger_Outs/Young_B593/outs/filtered_feature_bc_matrix/")
B593 <- CreateSeuratObject(counts = B593.data, project = "10X_MSC", min.cells = 20, min.features = 200)
#We calculate the percentage of
# mitochondrial genes here and store it in percent.mt.
B593[["percent.mt"]] <- PercentageFeatureSet(B593, pattern = "^MT-")

###
ggplot(B593@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(B593@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))

##!!remove contamination
B593 <- subset(B593, subset = nFeature_RNA > 4250 & nFeature_RNA < 7500 & percent.mt < 4.5 & nCount_RNA<55000 & nCount_RNA>0)
##Normalization
B593 <- NormalizeData(B593, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
B593 <- FindVariableFeatures(B593, selection.method = "vst", nfeatures = 1500)
VariableFeaturePlot(B593)
##ScaleData
B593 <- ScaleData(B593, vars.to.regress =c( "percent.mt",'nCount_RNA'))
##
B593=RenameCells(B593,add.cell.id = "Young_B593")
B593$Type='Young'
B593$group='B593'
B593$Sample='YoungBM3'

##Dimensional Reduction
B593 <- RunPCA(B593, npcs = 30, verbose = FALSE,features = VariableFeatures(object = B593))
B593 <- RunUMAP(B593, reduction = "pca", dims = 1:20)

B593<- FindNeighbors(B593, reduction = "pca", dims = 1:20)
B593<- FindClusters(B593, resolution = 0.4)

DimPlot(B593, reduction = "umap", label = TRUE, pt.size = 1.3)

####YoungBM4 B844####
B844.data <- Read10X(data.dir = "/data/NC_Revision_YoungAgeBM/CellRanger_Outs/Young_B844/outs/filtered_feature_bc_matrix/")
B844 <- CreateSeuratObject(counts = B844.data, project = "10X_MSC", min.cells = 20, min.features = 200)
#We calculate the percentage of
# mitochondrial genes here and store it in percent.mt.
B844[["percent.mt"]] <- PercentageFeatureSet(B844, pattern = "^MT-")

###
ggplot(B844@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(B844@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))

##!!remove contamination
B844 <- subset(B844, subset = nFeature_RNA > 3250 & nFeature_RNA < 7000 & percent.mt < 7 & nCount_RNA<60000 & nCount_RNA>10000)
##Normalization
B844 <- NormalizeData(B844, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
B844 <- FindVariableFeatures(B844, selection.method = "vst", nfeatures = 1500)
VariableFeaturePlot(B844)
##ScaleData
B844 <- ScaleData(B844, vars.to.regress =c( "percent.mt",'nCount_RNA'))
##
B844=RenameCells(B844,add.cell.id = "Young_B844")
B844$Type='Young'
B844$group='B844'
B844$Sample='YoungBM4'

##Dimensional Reduction
B844 <- RunPCA(B844, npcs = 30, verbose = FALSE,features = VariableFeatures(object = B844))
B844 <- RunUMAP(B844, reduction = "pca", dims = 1:20)

B844<- FindNeighbors(B844, reduction = "pca", dims = 1:20)
B844<- FindClusters(B844, resolution = 0.4)

DimPlot(B844, reduction = "umap", label = TRUE, pt.size = 1.3)

####AgedBM1 B801####
B801.data <- Read10X(data.dir = "/data/NC_Revision_YoungAgeBM/CellRanger_Outs/Aged_B801/outs/filtered_feature_bc_matrix/")
B801 <- CreateSeuratObject(counts = B801.data, project = "10X_MSC", min.cells = 20, min.features = 200)
#We calculate the percentage of
# mitochondrial genes here and store it in percent.mt.
B801[["percent.mt"]] <- PercentageFeatureSet(B801, pattern = "^MT-")

###
ggplot(B801@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(B801@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))

##!!remove contamination
B801 <- subset(B801, subset = nFeature_RNA > 3000 & nFeature_RNA < 6000 & percent.mt < 4.5 & nCount_RNA<45000 & nCount_RNA>10000)
##Normalization
B801 <- NormalizeData(B801, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
B801 <- FindVariableFeatures(B801, selection.method = "vst", nfeatures = 1500)
VariableFeaturePlot(B801)
##标准化---!! 3.0 scale data 默认只标准化var。genes （默认feature为find variable feature结果）
B801 <- ScaleData(B801, vars.to.regress =c( "percent.mt",'nCount_RNA'))
##
B801=RenameCells(B801,add.cell.id = "Aged_B801")
B801$Type='Aged'
B801$group='B801'
B801$Sample='AgedBM1'

##Dimensional Reduction
B801 <- RunPCA(B801, npcs = 30, verbose = FALSE,features = VariableFeatures(object = B801))
B801 <- RunUMAP(B801, reduction = "pca", dims = 1:20)

B801<- FindNeighbors(B801, reduction = "pca", dims = 1:20)
B801<- FindClusters(B801, resolution = 0.4)

DimPlot(B801, reduction = "umap", label = TRUE, pt.size = 1.3)
####AgedBM2 B842####
B842.data <- Read10X(data.dir = "/data/NC_Revision_YoungAgeBM/CellRanger_Outs/Aged_B842/outs/filtered_feature_bc_matrix/")
B842 <- CreateSeuratObject(counts = B842.data, project = "10X_MSC", min.cells = 20, min.features = 200)
#We calculate the percentage of
# mitochondrial genes here and store it in percent.mt.
B842[["percent.mt"]] <- PercentageFeatureSet(B842, pattern = "^MT-")

###
ggplot(B842@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(B842@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))

##!!remove contamination
B842 <- subset(B842, subset = nFeature_RNA > 1750 & nFeature_RNA < 5750 & percent.mt < 5 & nCount_RNA<35000 & nCount_RNA>5000)
##Normalization
B842 <- NormalizeData(B842, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
B842 <- FindVariableFeatures(B842, selection.method = "vst", nfeatures = 1500)
VariableFeaturePlot(B842)
##标准化---!! 3.0 scale data 默认只标准化var。genes （默认feature为find variable feature结果）
B842 <- ScaleData(B842, vars.to.regress =c( "percent.mt",'nCount_RNA'))
##
B842=RenameCells(B842,add.cell.id = "Aged_B842")
B842$Type='Aged'
B842$group='B842'
B842$Sample='AgedBM2'

##Dimensional Reduction
B842 <- RunPCA(B842, npcs = 30, verbose = FALSE,features = VariableFeatures(object = B842))
B842 <- RunUMAP(B842, reduction = "pca", dims = 1:20)

B842<- FindNeighbors(B842, reduction = "pca", dims = 1:20)
B842<- FindClusters(B842, resolution = 0.4)

DimPlot(B842, reduction = "umap", label = TRUE, pt.size = 1.3)
####AgedBM3 B857####
B857.data <- Read10X(data.dir = "/data/NC_Revision_YoungAgeBM/CellRanger_Outs/Aged_B857/outs/filtered_feature_bc_matrix/")
B857 <- CreateSeuratObject(counts = B857.data, project = "10X_MSC", min.cells = 20, min.features = 200)
#We calculate the percentage of
# mitochondrial genes here and store it in percent.mt.
B857[["percent.mt"]] <- PercentageFeatureSet(B857, pattern = "^MT-")

###
ggplot(B857@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,270000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(B857@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,270000,5000))+scale_y_continuous(breaks=seq(0,19000,250))+theme(axis.text.x= element_text(angle=90))

##!!remove contamination
B857 <- subset(B857, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 4 & nCount_RNA>0 & nCount_RNA<270000)
##Normalization
B857 <- NormalizeData(B857, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
B857 <- FindVariableFeatures(B857, selection.method = "vst", nfeatures = 1500)
VariableFeaturePlot(B857)
##标准化---!! 3.0 scale data 默认只标准化var。genes （默认feature为find variable feature结果）
B857 <- ScaleData(B857, vars.to.regress =c( "percent.mt",'nCount_RNA'))
##
B857=RenameCells(B857,add.cell.id = "Aged_B857")
B857$Type='Aged'
B857$group='B857'
B857$Sample='AgedBM3'

##Dimensional Reduction
B857 <- RunPCA(B857, npcs = 30, verbose = FALSE,features = VariableFeatures(object = B857))
B857 <- RunUMAP(B857, reduction = "pca", dims = 1:20)

B857<- FindNeighbors(B857, reduction = "pca", dims = 1:20)
B857<- FindClusters(B857, resolution = 0.4)

DimPlot(B857, reduction = "umap", label = TRUE, pt.size = 1.3)
####AgedBM4 B927####
B927.data <- Read10X(data.dir = "/data/NC_Revision_YoungAgeBM/CellRanger_Outs/Aged_B927/outs/filtered_feature_bc_matrix/")
B927 <- CreateSeuratObject(counts = B927.data, project = "10X_MSC", min.cells = 20, min.features = 200)
#We calculate the percentage of
# mitochondrial genes here and store it in percent.mt.
B927[["percent.mt"]] <- PercentageFeatureSet(B927, pattern = "^MT-")

###
ggplot(B927@meta.data,aes(x=nCount_RNA,y=percent.mt))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,25,1))+theme(axis.text.x= element_text(angle=90))
ggplot(B927@meta.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point(size=1,color='navy',shape=20)+scale_x_continuous(breaks=seq(0,170000,5000))+scale_y_continuous(breaks=seq(0,9000,250))+theme(axis.text.x= element_text(angle=90))

##!!remove contamination
B927 <- subset(B927, subset = nFeature_RNA > 3750 & nFeature_RNA < 8250 & percent.mt < 4 & nCount_RNA<110000 & nCount_RNA>15000)
##Normalization
B927 <- NormalizeData(B927, normalization.method = "LogNormalize", scale.factor = 10000)
##Identification of highly variable features (feature selection)
B927 <- FindVariableFeatures(B927, selection.method = "vst", nfeatures = 1500)
VariableFeaturePlot(B927)
##标准化---!! 3.0 scale data 默认只标准化var。genes （默认feature为find variable feature结果）
B927 <- ScaleData(B927, vars.to.regress =c( "percent.mt",'nCount_RNA'))
##
B927=RenameCells(B927,add.cell.id = "Aged_B927")
B927$Type='Aged'
B927$group='B927'

##Dimensional Reduction
B927 <- RunPCA(B927, npcs = 30, verbose = FALSE,features = VariableFeatures(object = B927))
B927 <- RunUMAP(B927, reduction = "pca", dims = 1:20)

B927<- FindNeighbors(B927, reduction = "pca", dims = 1:20)
B927<- FindClusters(B927, resolution = 0.4)

DimPlot(B927, reduction = "umap", label = TRUE, pt.size = 1.3)

