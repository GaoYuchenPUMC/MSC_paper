
##load loom files: Results from velocyto: Python##
ldat <- read.loom.matrices("/Data_MSC/Velocyto/merged_Loom/merged_All_12_MSCsamples.loom")
MSC= readRDS('/Data_MSC/MSC_intergated.rds')

##AD: MSC_AD:AAACCTGCAAAGGAAGx  MSC_AD2reseq:AGTAGTCTCGGCATCGx  MSC_AD3:CTACGTCCACAGCGTCx
##BM: MSC_BM:GGGCACTAGAACTGTAx  MSC_BM2reseq:CATCAAGCAGGCTGAAx  MSC_BM3reseq:AAACCTGCAAGGACTGx
##PM: MSC_PM_1:GTTCGGGGTTTGACTGx  MSC_PM2:AGGCCACCAAGCGCTCx   MSC_PM3:AAAGTAGCACCCATGGx
##UC: MSC_UL:AAGGCAGTCAAACGGGx  MSC_UL2:ACTGCTCAGATCGATAx  MSC_UL3:CATGACACAATCTGCAx

library(stringr)
BM1_origin=str_sub(row.names(subset(MSC@meta.data,group=='BM1')),start=5)
BM2_origin=str_sub(row.names(subset(MSC@meta.data,group=='BM2')),start=5)
BM3_origin=str_sub(row.names(subset(MSC@meta.data,group=='BM3')),start=5)
AD1_origin=str_sub(row.names(subset(MSC@meta.data,group=='AD1')),start=5)
AD2_origin=str_sub(row.names(subset(MSC@meta.data,group=='AD2')),start=5)
AD3_origin=str_sub(row.names(subset(MSC@meta.data,group=='AD3')),start=5)
PM1_origin=str_sub(row.names(subset(MSC@meta.data,group=='PM1')),start=5)
PM2_origin=str_sub(row.names(subset(MSC@meta.data,group=='PM2')),start=5)
PM3_origin=str_sub(row.names(subset(MSC@meta.data,group=='PM3')),start=5)
UC1_origin=str_sub(row.names(subset(MSC@meta.data,group=='UC1')),start=5)
UC2_origin=str_sub(row.names(subset(MSC@meta.data,group=='UC2')),start=5)
UC3_origin=str_sub(row.names(subset(MSC@meta.data,group=='UC3')),start=5)

##change barcode
BMbarcode1=intersect(paste('MSC_BM:',BM1_origin,'x',sep = ''),colnames(ldat$spliced))
BMbarcode2=intersect(paste('MSC_BM2reseq:',BM2_origin,'x',sep = ''),colnames(ldat$spliced))
BMbarcode3=intersect(paste('MSC_BM3reseq:',BM3_origin,'x',sep = ''),colnames(ldat$spliced))

ADbarcode1=intersect(paste('MSC_AD:',AD1_origin,'x',sep = ''),colnames(ldat$spliced))
ADbarcode2=intersect(paste('MSC_AD2reseq:',AD2_origin,'x',sep = ''),colnames(ldat$spliced))
ADbarcode3=intersect(paste('MSC_AD3:',AD3_origin,'x',sep = ''),colnames(ldat$spliced))

PMbarcode1=intersect(paste('MSC_PM_1:',PM1_origin,'x',sep = ''),colnames(ldat$spliced))
PMbarcode2=intersect(paste('MSC_PM2:',PM2_origin,'x',sep = ''),colnames(ldat$spliced))
PMbarcode3=intersect(paste('MSC_PM3:',PM3_origin,'x',sep = ''),colnames(ldat$spliced))

UCbarcode1=intersect(paste('MSC_UL:',UC1_origin,'x',sep = ''),colnames(ldat$spliced))
UCbarcode2=intersect(paste('MSC_UL2:',UC2_origin,'x',sep = ''),colnames(ldat$spliced))
UCbarcode3=intersect(paste('MSC_UL3:',UC3_origin,'x',sep = ''),colnames(ldat$spliced))

barcode=c(ADbarcode1,ADbarcode2,ADbarcode3,BMbarcode1,BMbarcode2,BMbarcode3,PMbarcode1,PMbarcode2,PMbarcode3,UCbarcode1,UCbarcode2,UCbarcode3)

emb=MSC@reductions$umap@cell.embeddings
row.names(emb)=barcode

cluster.label=Idents(MSC)
names(x = cluster.label) <- barcode

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = MSC)))
names(x = ident.colors) <- levels(x = MSC)
cell.colors <- ident.colors[Idents(object = MSC)]
names(x = cell.colors) <- barcode

pca=MSC@reductions$pca@cell.embeddings
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



