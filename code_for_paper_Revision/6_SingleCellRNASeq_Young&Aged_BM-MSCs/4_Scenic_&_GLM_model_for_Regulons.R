####Running SCENIC for Young Aged BM####

#1 Directories
setwd("/data/NC_Revision_YoungAgeBM/Scenic")  
dir.create("int")

#2 Input Expression matrix (on windows Seurat4)
library(Seurat)
#use raw count data
AllBM=readRDS('C:/Users/admin/Desktop/Manuscript for NM/20220826 NC revision/Single Cell YoungAge BM/AllBM6Cluster_Transferanchor.rds')
exprMat=as.matrix(AllBM@assays$RNA@counts)

cellInfo <- data.frame(seuratCluster=Idents(AllBM))
head(cellInfo)
cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "seuratCluster"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"

##subsample 1000 cells from each sample
meta=AllBM@meta.data

Young1_subsample=sample(row.names(meta[meta$group=='B469',]),1000)
Young2_subsample=sample(row.names(meta[meta$group=='B537',]),1000)
Young3_subsample=sample(row.names(meta[meta$group=='B593',]),1000)
Young4_subsample=sample(row.names(meta[meta$group=='B844',]),1000)
Aged1_subsample=sample(row.names(meta[meta$group=='B801',]),1000)
Aged2_subsample=sample(row.names(meta[meta$group=='B842',]),1000)
Aged3_subsample=sample(row.names(meta[meta$group=='B857',]),1000)
Aged4_subsample=sample(row.names(meta[meta$group=='B927',]),1000)

barcode_subsample=c(Young1_subsample,Young2_subsample,Young3_subsample,Young4_subsample,Aged1_subsample,Aged2_subsample,Aged3_subsample,Aged4_subsample)
exprMat_subsample=exprMat[,barcode_subsample]

cellInfo_subsample=cbind(barcode_subsample,as.vector(cellInfo[barcode_subsample,]))
row.names(cellInfo_subsample)=barcode_subsample
cellInfo_subsample=as.data.frame(cellInfo_subsample[,-1])
colnames(cellInfo_subsample)[1]="CellType"

colVars <- list(CellType=c("2"="forestgreen",
                           "1"="darkorange",
                           "4"="magenta4",
                           "5"="hotpink",
                           "0"="red3",
                           "6"="skyblue"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]

#4 Initialize SCENIC settings
library(SCENIC)
org="hgnc" # or hgnc--human, or dmel--fly  mgi--mouse
dbDir="/data/AliCloud_Data/Data_MSC/SCENIC_MouseBrain/Human_RcisTarget_Database" # RcisTarget databases location
myDatasetTitle="SCENIC of YoungAgedBM" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#5 Gene filter/selection
exprMat_subsample=read.csv('int/exprMat_subsample.csv')
row.names(exprMat_subsample)=exprMat_subsample$X
exprMat_subsample=exprMat_subsample[,-1]

exprMat_subsample=as.matrix(exprMat_subsample)
genesKept <- geneFiltering(exprMat_subsample, scenicOptions=scenicOptions,
                           minCountsPerGene=200,
                           minSamples=20)
##Before proceeding to the network inference, check whether any known relevant genes are filtered-out (if any relevant gene is missing, doublecheck
##whether the filters are appropiate):
interestingGenes <- c('GATA2','FOS','HMGA2','BRCA2','HMGA1','THAP11','DLX5')
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_subsample_filtered <- exprMat_subsample[genesKept, ]
dim(exprMat_subsample_filtered)

#6 Correlation
runCorrelation(exprMat_subsample_filtered, scenicOptions)

#7 GENIE3 
exprMat_subsample_filtered <- log2(exprMat_subsample_filtered+1)
# Run GENIE3 
runGenie3(exprMat_subsample_filtered, scenicOptions)

#8 Build and score the GRN (runSCENIC_…) The activity of the regulatory network, trained on this subset of 8000 cells, 
#  was evaluated on all the cells in the dataset with AUCell
exprMat=read.csv('int/exprMat.csv')
row.names(exprMat)=exprMat$X
exprMat=exprMat[,-1]

logMat <- log2(exprMat+1)
dim(logMat)
rm(exprMat)
logMat=as.matrix(logMat)

library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, logMat,skipHeatmap = TRUE,skipTsne = TRUE)

##10.Export to loom/SCope

# DGEM (Digital gene expression matrix)
# (non-normalized counts)
# exprMat <- get_dgem(open_loom(loomPath))
# dgem <- exprMat
# head(colnames(dgem))  #should contain the Cell ID/name

# Export:
library(SCopeLoomR)
scenicOptions@fileNames$output["loomFile",] <- "output/YoungAged_BM_SCENIC.loom"

fileName <- getOutName(scenicOptions, "loomFile")
scenicOptions@inputDatasetInfo$cellInfo = cellInfo
export2scope(scenicOptions, exprMat)


####Perform GLM analysis####
meta=read.csv('/data/NC_Revision_YoungAgeBM/Scenic/metaAllBM.csv',row.names = 'barcode')

##Cluster 123 vs Other##
meta_BM123vsOther=read.csv('/data/NC_Revision_YoungAgeBM/Scenic/meta_for_GLM_BM123vsOthers.csv',row.names = 'barcode')
AUC_BM123vsOther=AUC[,as.vector(row.names(meta_BM123vsOther))]
AUC_BM123vsOther=cbind(t(AUC_BM123vsOther),meta_BM123vsOther$BM123vsOthers)
colnames(AUC_BM123vsOther)[366]='BM123_vs_other'

colnames(AUC_BM123vsOther)=gsub("[ (.*)]", "", colnames(AUC_BM123vsOther))
AUC_BM123vsOther=as.data.frame(AUC_BM123vsOther)

compare_t_list=list()
for(i in colnames(AUC_BM123vsOther[,c(1:365)])){
  f=as.formula(paste(i,'BM123_vs_other',sep='~'))
  model=glm(f,data = AUC_BM123vsOther)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:365)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}

colnames(compare_t_mat)=row.names(AUC)
write.csv(t(compare_t_mat),'/data/NC_Revision_YoungAgeBM/Scenic/GLM/BM123_vs_other.csv')
rm(AUC_BM123vsOther,meta_BM123vsOther,compare_t_list,compare_t_mat)

##Cluster 4 vs 123
meta_BM4vsBM123=read.csv('/data/NC_Revision_YoungAgeBM/Scenic/meta_for_GLM_BM4vsBM123.csv',row.names = 'barcode')
AUC_BM4vsBM123=AUC[,as.vector(row.names(meta_BM4vsBM123))]
AUC_BM4vsBM123=cbind(t(AUC_BM4vsBM123),meta_BM4vsBM123$BM4vsBM123)
colnames(AUC_BM4vsBM123)[366]='BM4_vs_BM123'

colnames(AUC_BM4vsBM123)=gsub("[ (.*)]", "", colnames(AUC_BM4vsBM123))
AUC_BM4vsBM123=as.data.frame(AUC_BM4vsBM123)

compare_t_list=list()
for(i in colnames(AUC_BM4vsBM123[,c(1:365)])){
  f=as.formula(paste(i,'BM4_vs_BM123',sep='~'))
  model=glm(f,data = AUC_BM4vsBM123)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:365)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}
colnames(compare_t_mat)=row.names(AUC)
write.csv(t(compare_t_mat),'/data/NC_Revision_YoungAgeBM/Scenic/GLM/BM4_vs_BM123.csv')
rm(AUC_BM4vsBM123,meta_BM4vsBM123,compare_t_list,compare_t_mat)

##Cluster 5 vs 4
meta_BM5vsBM4=read.csv('/data/NC_Revision_YoungAgeBM/Scenic/meta_for_GLM_BM5vsBM4.csv',row.names = 'barcode')
AUC_BM5vsBM4=AUC[,as.vector(row.names(meta_BM5vsBM4))]
AUC_BM5vsBM4=cbind(t(AUC_BM5vsBM4),meta_BM5vsBM4$BM5vsBM4)
colnames(AUC_BM5vsBM4)[366]='BM5_vs_BM4'

colnames(AUC_BM5vsBM4)=gsub("[ (.*)]", "", colnames(AUC_BM5vsBM4))
AUC_BM5vsBM4=as.data.frame(AUC_BM5vsBM4)

compare_t_list=list()
for(i in colnames(AUC_BM5vsBM4[,c(1:365)])){
  f=as.formula(paste(i,'BM5_vs_BM4',sep='~'))
  model=glm(f,data = AUC_BM5vsBM4)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:365)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}
colnames(compare_t_mat)=row.names(AUC)
write.csv(t(compare_t_mat),'/data/NC_Revision_YoungAgeBM/Scenic/GLM/BM5_vs_BM4.csv')
rm(AUC_BM5vsBM4,meta_BM5vsBM4,compare_t_list,compare_t_mat)


##Cluster 6 vs 4
meta_BM6vsBM4=read.csv('/data/NC_Revision_YoungAgeBM/Scenic/meta_for_GLM_BM6vsBM4.csv',row.names = 'barcode')
AUC_BM6vsBM4=AUC[,as.vector(row.names(meta_BM6vsBM4))]
AUC_BM6vsBM4=cbind(t(AUC_BM6vsBM4),meta_BM6vsBM4$BM6vsBM4)
colnames(AUC_BM6vsBM4)[366]='BM6_vs_BM4'

colnames(AUC_BM6vsBM4)=gsub("[ (.*)]", "", colnames(AUC_BM6vsBM4))
AUC_BM6vsBM4=as.data.frame(AUC_BM6vsBM4)

compare_t_list=list()
for(i in colnames(AUC_BM6vsBM4[,c(1:365)])){
  f=as.formula(paste(i,'BM6_vs_BM4',sep='~'))
  model=glm(f,data = AUC_BM6vsBM4)
  compare_t_list[[i]]=summary(model)$coefficients[2,]
}


compare_t_mat=as.data.frame(compare_t_list[1])
for(i in c(2:365)){
  data=compare_t_list[i]
  compare_t_mat=cbind(compare_t_mat,data)
}

colnames(compare_t_mat)=row.names(AUC)
write.csv(t(compare_t_mat),'/data/NC_Revision_YoungAgeBM/Scenic/GLM/BM6_vs_BM4.csv')
rm(AUC_BM6vsBM4,meta_BM6vsBM4,compare_t_list,compare_t_mat)

