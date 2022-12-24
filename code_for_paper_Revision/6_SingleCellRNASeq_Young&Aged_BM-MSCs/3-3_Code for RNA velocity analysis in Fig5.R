####Velocyto R
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

AllBM=readRDS('/data/Single Cell YoungAge BM/AllBM6Cluster_Transferanchor.rds')
ldat <- read.loom.matrices("/data/Single Cell YoungAge BM/Velocyto/merged_BM.loom")


## Aged_B801:ACGATGTAGTACTGTCx  Aged_B842:ATCGCCTTCCGAACGCx  Aged_B857:ATCAGGTGTACGTTCAx   Aged_B927:ATTCAGGTCGGTTAGTx
## Young_B469:CAACCTCAGGAGTACCx  Young_B537:ATTACTCTCGAGCCTGx  Young_B593:AACGAAAAGAGAACCCx  Young_B844:ACTCCCACACGCGTGTx

library(stringr)
Age801_origin=str_sub(row.names(subset(AllBM@meta.data,group=='B801')),start=11)
Age842_origin=str_sub(row.names(subset(AllBM@meta.data,group=='B842')),start=11)
Age857_origin=str_sub(row.names(subset(AllBM@meta.data,group=='B857')),start=11)
Age927_origin=str_sub(row.names(subset(AllBM@meta.data,group=='B927')),start=11)

Young469_origin=str_sub(row.names(subset(AllBM@meta.data,group=='B469')),start=12)
Young537_origin=str_sub(row.names(subset(AllBM@meta.data,group=='B537')),start=12)
Young593_origin=str_sub(row.names(subset(AllBM@meta.data,group=='B593')),start=12)
Young844_origin=str_sub(row.names(subset(AllBM@meta.data,group=='B844')),start=12)

##Final Barcode
Age801_barcode=intersect(paste('Aged_B801:',Age801_origin,'x',sep = ''),colnames(ldat$spliced))
Age842_barcode=intersect(paste('Aged_B842:',Age842_origin,'x',sep = ''),colnames(ldat$spliced))
Age857_barcode=intersect(paste('Aged_B857:',Age857_origin,'x',sep = ''),colnames(ldat$spliced))
Age927_barcode=intersect(paste('Aged_B927:',Age927_origin,'x',sep = ''),colnames(ldat$spliced))

Young469_barcode=intersect(paste('Young_B469:',Young469_origin,'x',sep = ''),colnames(ldat$spliced))
Young537_barcode=intersect(paste('Young_B537:',Young537_origin,'x',sep = ''),colnames(ldat$spliced))
Young593_barcode=intersect(paste('Young_B593:',Young593_origin,'x',sep = ''),colnames(ldat$spliced))
Young844_barcode=intersect(paste('Young_B844:',Young844_origin,'x',sep = ''),colnames(ldat$spliced))

####1. Aged-BM Group####
barcode_age=c(Age801_barcode,Age842_barcode,Age857_barcode,Age927_barcode)
AllBM_age=subset(AllBM,Type=='Aged')

emb=AllBM_age@reductions$ref.umap@cell.embeddings
row.names(emb)=barcode_age

cluster.label=Idents(AllBM_age)
names(x = cluster.label) <- barcode_age

##color
AllBM_age@meta.data$Cluster = factor(AllBM_age@meta.data$Cluster,levels=c('BM1','BM2','BM3','BM4','BM5','BM6'), ordered=TRUE)

ident.colors <- c("#72BE9C",'orange2','#5088BF',"#7CAE00",'#009251','#EC494E')
names(x = ident.colors) <- levels(x = AllBM_age)
cell.colors <- ident.colors[Idents(object = AllBM_age)]
names(x = cell.colors) <- barcode_age

##splice unsplice
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode_age]
nmat <- nmat[,barcode_age]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))
##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,fit.quantile=fit.quantile,n.cores = 1)
##Plot AgedBM Velocity
show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)


####2. Young-BM Group####
barcode_Young=c(Young469_barcode,Young537_barcode,Young593_barcode,Young844_barcode)
AllBM_Young=subset(AllBM,Type=='Young')

emb=AllBM_Young@reductions$ref.umap@cell.embeddings
row.names(emb)=barcode_Young

cluster.label=Idents(AllBM_Young)
names(x = cluster.label) <- barcode_Young

##color
AllBM_Young@meta.data$Cluster = factor(AllBM_Young@meta.data$Cluster,levels=c('BM1','BM2','BM3','BM4','BM5','BM6'), ordered=TRUE)

ident.colors <- c("#72BE9C",'orange2','#5088BF',"#7CAE00",'#009251','#EC494E')
names(x = ident.colors) <- levels(x = AllBM_Young)
cell.colors <- ident.colors[Idents(object = AllBM_Young)]
names(x = cell.colors) <- barcode_Young

##splice unsplice
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode_Young]
nmat <- nmat[,barcode_Young]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))
##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,fit.quantile=fit.quantile,n.cores = 1)
##Plot YoungBM Velocity
show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)




