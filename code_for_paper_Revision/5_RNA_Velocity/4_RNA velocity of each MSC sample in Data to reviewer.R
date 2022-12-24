####RNA velocity of AD1 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
AD1_origin=str_sub(row.names(subset(MSC@meta.data,group=='AD1')),start=5)

##change barcode
ADbarcode1=intersect(paste('MSC_AD:',AD1_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=ADbarcode1

AD1=subset(MSC,group=='AD1')
emb=AD1@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(AD1)
names(x = cluster.label) <- barcode

AD1@meta.data$Cluster = factor(AD1@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = AD1)
cell.colors <- ident.colors[Idents(object = AD1)]
names(x = cell.colors) <- barcode

pca=AD1@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of AD2 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
AD2_origin=str_sub(row.names(subset(MSC@meta.data,group=='AD2')),start=5)

##change barcode
ADbarcode2=intersect(paste('MSC_AD2reseq:',AD2_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=ADbarcode2

AD2=subset(MSC,group=='AD2')
emb=AD2@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(AD2)
names(x = cluster.label) <- barcode

AD2@meta.data$Cluster = factor(AD2@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = AD2)
cell.colors <- ident.colors[Idents(object = AD2)]
names(x = cell.colors) <- barcode

pca=AD2@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of AD3 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
AD3_origin=str_sub(row.names(subset(MSC@meta.data,group=='AD3')),start=5)

##change barcode
ADbarcode3=intersect(paste('MSC_AD3:',AD3_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=ADbarcode3

AD3=subset(MSC,group=='AD3')
emb=AD3@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(AD3)
names(x = cluster.label) <- barcode

AD3@meta.data$Cluster = factor(AD3@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = AD3)
cell.colors <- ident.colors[Idents(object = AD3)]
names(x = cell.colors) <- barcode

pca=AD3@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)

####RNA velocity of BM1 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
BM1_origin=str_sub(row.names(subset(MSC@meta.data,group=='BM1')),start=5)

##change barcode
BMbarcode1=intersect(paste('MSC_BM:',BM1_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=BMbarcode1

BM1=subset(MSC,group=='BM1')
emb=BM1@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(BM1)
names(x = cluster.label) <- barcode

BM1@meta.data$Cluster = factor(BM1@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = BM1)
cell.colors <- ident.colors[Idents(object = BM1)]
names(x = cell.colors) <- barcode

pca=BM1@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of BM2 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
BM2_origin=str_sub(row.names(subset(MSC@meta.data,group=='BM2')),start=5)

##change barcode
BMbarcode2=intersect(paste('MSC_BM2reseq:',BM2_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=BMbarcode2

BM2=subset(MSC,group=='BM2')
emb=BM2@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(BM2)
names(x = cluster.label) <- barcode

BM2@meta.data$Cluster = factor(BM2@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = BM2)
cell.colors <- ident.colors[Idents(object = BM2)]
names(x = cell.colors) <- barcode

pca=BM2@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of BM3 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
BM3_origin=str_sub(row.names(subset(MSC@meta.data,group=='BM3')),start=5)

##change barcode
BMbarcode3=intersect(paste('MSC_BM3reseq:',BM3_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=BMbarcode3

BM3=subset(MSC,group=='BM3')
emb=BM3@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(BM3)
names(x = cluster.label) <- barcode

BM3@meta.data$Cluster = factor(BM3@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = BM3)
cell.colors <- ident.colors[Idents(object = BM3)]
names(x = cell.colors) <- barcode

pca=BM3@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of UC1 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
UC1_origin=str_sub(row.names(subset(MSC@meta.data,group=='UC1')),start=5)

##change barcode
UCbarcode1=intersect(paste('MSC_UL:',UC1_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=UCbarcode1

UC1=subset(MSC,group=='UC1')
emb=UC1@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(UC1)
names(x = cluster.label) <- barcode

UC1@meta.data$Cluster = factor(UC1@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = UC1)
cell.colors <- ident.colors[Idents(object = UC1)]
names(x = cell.colors) <- barcode

pca=UC1@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of UC2 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
UC2_origin=str_sub(row.names(subset(MSC@meta.data,group=='UC2')),start=5)

##change barcode
UCbarcode2=intersect(paste('MSC_UL2:',UC2_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=UCbarcode2

UC2=subset(MSC,group=='UC2')
emb=UC2@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(UC2)
names(x = cluster.label) <- barcode

UC2@meta.data$Cluster = factor(UC2@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = UC2)
cell.colors <- ident.colors[Idents(object = UC2)]
names(x = cell.colors) <- barcode

pca=UC2@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of UC3 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
UC3_origin=str_sub(row.names(subset(MSC@meta.data,group=='UC3')),start=5)

##change barcode
UCbarcode3=intersect(paste('MSC_UL3:',UC3_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=UCbarcode3

UC3=subset(MSC,group=='UC3')
emb=UC3@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(UC3)
names(x = cluster.label) <- barcode

UC3@meta.data$Cluster = factor(UC3@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = UC3)
cell.colors <- ident.colors[Idents(object = UC3)]
names(x = cell.colors) <- barcode

pca=UC3@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of PM1 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
PM1_origin=str_sub(row.names(subset(MSC@meta.data,group=='PM1')),start=5)

##change barcode
PMbarcode1=intersect(paste('MSC_PM_1:',PM1_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=PMbarcode1

PM1=subset(MSC,group=='PM1')
emb=PM1@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(PM1)
names(x = cluster.label) <- barcode

PM1@meta.data$Cluster = factor(PM1@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = PM1)
cell.colors <- ident.colors[Idents(object = PM1)]
names(x = cell.colors) <- barcode

pca=PM1@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of PM2 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
PM2_origin=str_sub(row.names(subset(MSC@meta.data,group=='PM2')),start=5)

##change barcode
PMbarcode2=intersect(paste('MSC_PM2:',PM2_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=PMbarcode2

PM2=subset(MSC,group=='PM2')
emb=PM2@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(PM2)
names(x = cluster.label) <- barcode

PM2@meta.data$Cluster = factor(PM2@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = PM2)
cell.colors <- ident.colors[Idents(object = PM2)]
names(x = cell.colors) <- barcode

pca=PM2@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
####RNA velocity of PM3 MSC####

##load loom files: Results from velocyto: Python##
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

library(stringr)
PM3_origin=str_sub(row.names(subset(MSC@meta.data,group=='PM3')),start=5)

##change barcode
PMbarcode3=intersect(paste('MSC_PM3:',PM3_origin,'x',sep = ''),colnames(ldat$spliced))
barcode=PMbarcode3

PM3=subset(MSC,group=='PM3')
emb=PM3@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(PM3)
names(x = cluster.label) <- barcode

PM3@meta.data$Cluster = factor(PM3@meta.data$Cluster,levels=c('C1',"C2",'C3','C4','C5','C6','C7'), ordered=TRUE)
ident.colors <- c('#009251',"#7CAE00","#72BE9C",'#5088BF','orange2',"#FF5BC6","#EC494E")
names(x = ident.colors) <- levels(x = PM3)
cell.colors <- ident.colors[Idents(object = PM3)]
names(x = cell.colors) <- barcode

pca=PM3@reductions$pca@cell.embeddings
row.names(pca)=barcode
cell.dist=as.dist(1-armaCor(t(pca)))

##splice unsplice files
emat <- ldat$spliced
nmat <- ldat$unspliced
emat <- emat[,barcode]
nmat <- nmat[,barcode]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

##Calculate Velocity
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores = 20)

show.velocity.on.embedding.cor(emb,rvel.cd,n=250,scale='sqrt',cell.colors=ac(cell.colors),
                               cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=1,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0,n.cores = 20)
