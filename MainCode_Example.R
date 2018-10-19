# Main code example

# load the libraries
library(stringr)
library(CoRegNet)
library(doParallel)
library(igraph)
library(lsei)
library(lars)
library(matrixStats)
library(RCurl)
library(limma)
library(impute)
# library(scoop) # does not work yet
               # should be installed by another way
# source("/Users/Magali/Desktop/recherche/LIONS/codejulien/scoop/R/scoop.R")

# Define your path and target directory
path <- "/Users/mchampion/Desktop/LIONS/"
pathEM <- paste0(path,"algoEM/")
TargetDirectory <- paste0(path,"results/")
source(paste0(path,"LIONS/MainCode.R"))
source(paste0(path,"LIONS/TCGA_data.R"))

# load the data
# first type of data (CIT)
MA_normal <- read.table(paste0(path,"Data/expression/NormalExpressionData_AllGenes.txt"))
MA_cancer <- read.table(paste0(path,"Data/expression/TumorExpressionData_AllGenes.txt"))
rownames(MA_cancer) <- str_replace_all(rownames(MA_cancer),"[A-Z]","") # anonymous data

# second type of data (TCGA)
DataDirectory <- paste0(path,"Data/TCGA")
load(paste0(path,"Data/TCGA/BatchData.rda"))
DataSetDirectories <- Download_CancerSite(CancerSite = "BLCA", TargetDirectory = DataDirectory,downloadData = FALSE)
ProcessedData <- Preprocess_CancerSite(CancerSite = "BLCA",DataSetDirectories = DataSetDirectories)
load(paste0(DataDirectory,"/ProcessedData_BLCA.Rdata"))
MA <- t(ProcessedData$MA_TCGA)

# load the different subtypes
Subtypes <- read.table(paste0(path,"Data/Subtypes/Classif_TCGA_mars_2018.csv"),sep=",",header=TRUE)
samples <- as.character(Subtypes$ID)
samples <- substr(samples,1,nchar(samples)-1)
Subtypes <- Subtypes[,4]
Subtypes <- as.character(Subtypes)
names(Subtypes) <- samples
Subtypes <- Subtypes[rownames(MA)]

# CNV data 
load(paste0(DataDirectory,"/ProcessedData_BLCA.Rdata"))
CNV <- ProcessedData$CNV_TCGA_thresholded
colnames(CNV) <- paste0(colnames(CNV),"-01")
CNV <- t(CNV)

# list of transcription factors
TFs <- read.table(paste0(path,"Data/TCGA/AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1 

# set parameters
VarMax <- 0.75 # filter for genes, keep only the top VarMax % variant genes

# Run the main code (using hLicorn for reconstucting the network)
subtypes <- names(which(Subtypes==names(table(Subtypes))[1]))
Results <- LIONS_Main_Code(MA_cancer = MA[subtypes,],
                           MA_normal = MA[!(rownames(MA)%in%subtypes),],
                           TFs = TFs,Method="hLICORN",
                           CNV_matrix = CNV,CNV_correction = TRUE,
                           LicornThresholds=list(minCoregSupport=0.5,searchThresh=0.5),
                           TargetDirectory = TargetDirectory, pathEM = pathEM)