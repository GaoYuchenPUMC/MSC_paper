library(monocle)
library(Seurat)
library(ggplot2)

##1.Create the msc-monocle file
MSC=readRDS('/data/MSC_integrated.rds')
exprs=as.matrix(MSC@assays$RNA@data)
phData=MSC@meta.data
feaData=as.data.frame(row.names(exprs))
row.names(feaData)=feaData[,1]
colnames(feaData)[1]='gene_short_name'
pd <- new("AnnotatedDataFrame", data = phData)
fd <- new("AnnotatedDataFrame", data = feaData)
MSC_monocle <- newCellDataSet(exprs,phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())

##2.Estimate size factors and dispersions 
MSC_monocle <- estimateSizeFactors(MSC_monocle)
MSC_monocle <- estimateDispersions(MSC_monocle)

##4.Filtering cells and genes
MSC_monocle <- detectGenes(MSC_monocle, min_expr = 0.1)
print(head(fData(MSC_monocle)))
expressed_genes <- row.names(subset(fData(MSC_monocle),num_cells_expressed >= 10))

##6.Trajectory step 1: choose genes that define a cell's progress 
diff_test_res <- differentialGeneTest(MSC_monocle[expressed_genes,],fullModelFormulaStr = "~Cluster",cores = 20)
ordering_genes <-row.names(subset(diff_test_res,qval<1e-167)) ###The top 1250 DiffGenes

MSC_monocle <- setOrderingFilter(MSC_monocle, ordering_genes)
plot_ordering_genes(MSC_monocle)
##7.Trajectory step 2: reduce data dimensionality 
MSC_monocle <- reduceDimension(MSC_monocle, max_components = 2,method = 'DDRTree',cores=20)
##8.Trajectory step 3: order cells along the trajectory 
MSC_monocle <- orderCells(MSC_monocle) 
plot_cell_trajectory(MSC_monocle, color_by = "Cluster",cell_size = 1)+scale_color_manual(values = c('#009251',"#7CAE00","#72BE9C",'#5088BF',"#FF5BC6",'orange2',"#EC494E"))

##9.Set The Root Stage
MSC_monocle <-readRDS('/data/gaoyuchen/NC_Revision/Reviewer2_Monocle/monocle_order.rds')
plot_cell_trajectory(MSC_monocle, color_by = "State",cell_size = 1)
plot_cell_trajectory(MSC_monocle, color_by = "Pseudotime",cell_size = 1)

 ## We set the root based on where certain marker genes are expressed 
 ## Based on the biological knowledge of the system, 
 ## a highly proliferative non-senescent MSC population(C1 C2 C3)
 ## express high levels of proliferation markers. These cells are mainly 
 ## distributed in State 5 State 6 and State 7 after pseudotiem analysis,and we set the root state to State 7
MSC_monocle <- orderCells(MSC_monocle, root_state =7)
plot_cell_trajectory(MSC_monocle, color_by = "Pseudotime")

genes_curve=c('DNMT1','MKI67','CKS2','TOP2A','LMNB1','CDKN1A','TP53','COL3A1','IGFBP5')
plot_genes_in_pseudotime(MSC_monocle[genes_curve,],color_by = 'Cluster',ncol = 2,panel_order = c('DNMT1','MKI67','CKS2','TOP2A','LMNB1','CDKN1A','TP53','COL3A1','IGFBP5'))+scale_color_manual(values = c('#009251',"#7CAE00","#72BE9C",'#5088BF',"#FF5BC6",'orange2',"#EC494E"))

##10.Extract the information of Pseudotime point of each cell, draw the pseudotime heatmap of SCENIC
pseudotime_cell=MSC_monocle@phenoData@data

library(Seurat)
library(stringi)
library(stringr)
library(pheatmap)

regulonAUC=readRDS('/data/AliCloud_Data/Data_MSC/SCENIC_06_07_Cellranger12Aggr/int/3.4_regulonAUC.Rds')
AUC_eachCell=regulonAUC@assays@data@listData$AUC

AD1barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='AD1')),start=5),'-1',sep=''),colnames(AUC_eachCell))
AD2barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='AD2')),start=5),'-2',sep=''),colnames(AUC_eachCell))
AD3barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='AD3')),start=5),'-3',sep=''),colnames(AUC_eachCell))
BM1barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='BM1')),start=5),'-4',sep=''),colnames(AUC_eachCell))
BM2barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='BM2')),start=5),'-5',sep=''),colnames(AUC_eachCell))
BM3barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='BM3')),start=5),'-6',sep=''),colnames(AUC_eachCell))
PM1barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='PM1')),start=5),'-7',sep=''),colnames(AUC_eachCell))
PM2barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='PM2')),start=5),'-8',sep=''),colnames(AUC_eachCell))
PM3barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='PM3')),start=5),'-9',sep=''),colnames(AUC_eachCell))
UC1barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='UC1')),start=5),'-10',sep=''),colnames(AUC_eachCell))
UC2barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='UC2')),start=5),'-11',sep=''),colnames(AUC_eachCell))
UC3barcode=intersect(paste(str_sub(row.names(subset(pseudotime_cell,group=='UC3')),start=5),'-12',sep=''),colnames(AUC_eachCell))

row.names(pseudotime_cell)=c(AD1barcode,AD2barcode,AD3barcode,BM1barcode,BM2barcode,BM3barcode,
                             PM1barcode,PM2barcode,PM3barcode,UC1barcode,UC2barcode,UC3barcode)
pseudotime_cell=pseudotime_cell[order(pseudotime_cell$Pseudotime),]
barcode=row.names(pseudotime_cell)

Route_Gene=c('CTCF (16g)','HMG20B (10g)','EZH2 (30g)','ZEB1_extended (16g)','FOSL1 (13g)','E2F7 (246g)',
             'MYBL1 (453g)','E2F8 (151g)','E2F1 (1670g)','SOX11 (91g)','MYBL2 (1294g)','RAD21_extended (3818g)',
             'HDAC2_extended (4300g)','TFDP1 (3464g)','SMARCA4_extended (4683g)','TBX2_extended (23g)','HMGA2 (102g)',
             'ZNF121 (16g)','XBP1 (318g)','FOXP1_extended (18g)','NFE2L2 (15g)',
             'RORA (11g)','MLXIP (11g)','KLF10 (18g)','FOS (15g)','MEF2A (68g)','KLF4 (20g)','EGR3 (22g)',
             'FOXS1_extended (32g)','FOXC1 (15g)','TEAD1 (48g)','EMX2_extended (32g)','KLF6_extended (29g)',
             'TCF4 (27g)','STAT1 (632g)','FOXK1 (25g)','MAFG_extended (20g)','FOSB (88g)','JUNB (131g)',
             'MAF (29g)','EGR1 (79g)','CEBPB (59g)','SOX4 (34g)','RUNX2 (15g)')

AUC_route=AUC_eachCell[Route_Gene,barcode]

annotation_col =cbind(as.data.frame(pseudotime_cell[barcode,]),c(1:45955))
annotation_col =annotation_col[,c(11,14)]
colnames(annotation_col)[2]='Pseudotime'

ann_colors = list(Cluster = c(C1='#009251',C2="#7CAE00",C3="#72BE9C",C4='#5088BF',C7="#EC494E",C5='orange2',C6="#FF5BC6"),
                  Pseudotime=viridis::viridis(20))

pheatmap((AUC_route),show_colnames = F,cluster_rows = F,cluster_cols = F,scale = 'row',
         color = colorRampPalette(c("#1677B2", "white", "#CB0826"))(100),
         border=FALSE,breaks=unique(c(seq(-2.5,2.5, length=100))),annotation_col = annotation_col,annotation_colors = ann_colors)
