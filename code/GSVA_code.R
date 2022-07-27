###GSVA codes###
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)

##1.  Import of the reference gene set and Seurat MSC rds file

genesets <- getGmt('/public/home/gaoyuchen/GSVA_test/GSVA_Allpathway.gmt')
MSC=readRDS('/public/home/gaoyuchen/GSVA_test/MSC_integrated.rds')

##2. Obtain normalized count matrix
MSC = NormalizeData(MSC, normalization.method = "LogNormalize", scale.factor = 10000)
Cell_Gene_Matrix=MSC@assays$RNA@data

##3. GSVA analysis
GSVA <- gsva(as.matrix(Cell_Gene_Matrix), genesets, min.sz=10, max.sz=1000, abs.ranking=FALSE, verbose=TRUE,parallel.sz=5L)

##4. Differential analysis
meta=read.csv('/data/gaoyuchen/new_metadata.csv',row.names = 'barcode')
row.names(meta)=colnames(GSVA)

GSVA_Seurat <- CreateSeuratObject(counts = GSVA, meta.data = meta, project = "GSVA_singleCell")
Idents(GSVA_Seurat)='Cluster'

GSVA_subcluster=FindAllMarkers(GSVA_Seurat, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

