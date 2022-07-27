#1 Directories
setwd("/Data_MSC/SCENIC_5_27_MSC/")
dir.create("int")

#2 Input Expression matrix
library(Seurat)

MSC=readRDS('/Data_MSC/MSC_integrated.rds')
exprMat=as.matrix(MSC@assays$RNA@counts)

cellInfo <- data.frame(seuratCluster=Idents(MSC))
head(cellInfo)
cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "seuratCluster"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"

##3 subset sample: 1000 cells per sample

meta=MSC@meta.data

AD1_subsample=sample(row.names(meta[meta$group=='AD1',]),1000)
AD2_subsample=sample(row.names(meta[meta$group=='AD2',]),1000)
AD3_subsample=sample(row.names(meta[meta$group=='AD3',]),1000)
BM1_subsample=sample(row.names(meta[meta$group=='BM1',]),1000)
BM2_subsample=sample(row.names(meta[meta$group=='BM2',]),1000)
BM3_subsample=sample(row.names(meta[meta$group=='BM3',]),1000)
PM1_subsample=sample(row.names(meta[meta$group=='PM1',]),1000)
PM2_subsample=sample(row.names(meta[meta$group=='PM2',]),1000)
PM3_subsample=sample(row.names(meta[meta$group=='PM3',]),1000)
UC1_subsample=sample(row.names(meta[meta$group=='UC1',]),1000)
UC2_subsample=sample(row.names(meta[meta$group=='UC2',]),1000)
UC3_subsample=sample(row.names(meta[meta$group=='UC3',]),1000)

barcode_subsample=c(AD1_subsample,AD2_subsample,AD3_subsample,BM1_subsample,BM2_subsample,BM3_subsample,PM1_subsample,PM2_subsample,PM3_subsample,UC1_subsample,UC2_subsample,UC3_subsample)
exprMat_subsample=exprMat[,barcode_subsample]

cellInfo_subsample=cbind(barcode_subsample,as.vector(cellInfo[barcode_subsample,]))
row.names(cellInfo_subsample)=barcode_subsample
cellInfo_subsample=as.data.frame(cellInfo_subsample[,-1])
colnames(cellInfo_subsample)[1]="CellType"

saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(cellInfo_subsample, file="int/cellInfo_subsample.Rds")

colVars <- list(CellType=c("2"="forestgreen",
                           "1"="darkorange",
                           "4"="magenta4",
                           "5"="hotpink",
                           "0"="red3",
                           "6"="skyblue",
                           "3"="darkblue",'7'='black'))

colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")

##4 Initialize SCENIC settings
library(SCENIC)
org="hgnc" # or hgnc--human, or dmel--fly  mgi--mouse
dbDir="/Data_MSC/SCENIC_MouseBrain/Human_RcisTarget_Database" # RcisTarget databases location
myDatasetTitle="SCENIC example on MSC" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=30)
# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

##5 Gene filter/selection
# (Adjust minimum values according to your dataset) 

exprMat_subsample=as.matrix(exprMat_subsample)
genesKept <- geneFiltering(exprMat_subsample, scenicOptions=scenicOptions,
                           minCountsPerGene=200,
                           minSamples=20)
##Before proceeding to the network inference, check whether any known relevant genes are filtered-out (if any relevant gene is missing, doublecheck
##whether the filters are appropiate):
interestingGenes <- c("Sox9", "Sox10", "Dlx5")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_subsample_filtered <- exprMat_subsample[genesKept, ]
dim(exprMat_subsample_filtered)

##6 Correlation
runCorrelation(exprMat_subsample_filtered, scenicOptions)

##7 GENIE (To run GRNBoost (in Python) instead of GENIE3. See ?exportsForGRNBoost for details)

## If launched in a new session, you will need to reload...
# setwd("...")
# loomPath <- "..."
# loom <- open_loom(loomPath, mode="r")
# exprMat <- get_dgem(loom)
# close_loom(loom)
# genesKept <- loadInt(scenicOptions, "genesKept")
# exprMat_filtered <- exprMat[genesKept,]
# library(SCENIC)
# scenicOptions <- readRDS("int/scenicOptions.Rds")
# Optional: add log (if it is not logged/normalized already)
exprMat_subsample_filtered <- log2(exprMat_subsample_filtered+1)
# Run GENIE3 
runGenie3(exprMat_subsample_filtered, scenicOptions)

##8 Build and score the GRN (runSCENIC_…) 
# Optional: log expression (for TF expression plot, it does not affect any other calculation)
MSC=readRDS('/Data_MSC/MSC_12sampleCCA_515.rds')
logMat=as.matrix(MSC@assays$RNA@counts)
dim(logMat)

library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, logMat)

##9 Clustering / dimensionality reduction on the regulon activity

  # For toy dataset  nPcs <- c(5)
 nPcs <- c(5,15,30,50)

scenicOptions@settings$seed <- 123 # same seed for all of them
scenicOptions@settings$nCores <- 10
  # Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,30,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,30,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC") 


  ## Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(4,4))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.4)

  ## Using only "high-confidence" regulons (normally similar)
par(mfrow=c(length(nPcs), 4))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.4)

  #The chosen t-SNE can then be saved as default to use for plots (can also be “binary”, see below):

scenicOptions@settings$defaultTsne$aucType <- "AUC"  ##or "oHC_AUC"
scenicOptions@settings$defaultTsne$dims <- 15
scenicOptions@settings$defaultTsne$perpl <- 50
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

##10 Optional steps: Binarize the network activity (regulon on/off)
MSC=readRDS('/Data_MSC/MSC_12sampleCCA_515.rds')
logMat=as.matrix(MSC@assays$RNA@counts)
dim(logMat)

aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

  # Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
  ##Once you have optimized the thresholds, run runSCENIC_4_aucell_binarize to binarize the AUC, and generate some extra figures and clusterings:
  # scenicOptions@settings$devType="png"
runSCENIC_4_aucell_binarize(scenicOptions)




##10.Export to loom/SCope

 # DGEM (Digital gene expression matrix)
 # (non-normalized counts)
 # exprMat <- get_dgem(open_loom(loomPath))
 # dgem <- exprMat
 # head(colnames(dgem))  #should contain the Cell ID/name

 # Export:
library(SCopeLoomR)
scenicOptions@fileNames$output["loomFile",] <- "output/MSC_12sample_SCENIC.loom"

fileName <- getOutName(scenicOptions, "loomFile")
scenicOptions@inputDatasetInfo$cellInfo = cellInfo
export2scope(scenicOptions, exprMat)

