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
path <- "/Users/Magali/Desktop/recherche/LIONS/"
pathEM <- paste0(path,"algoEM/")
TargetDirectory <- paste0(path,"results/")
source(paste0(path,"LIONS_project/MainCode.R"))
source(paste0(path,"LIONS_project/TCGA_data.R"))

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
MA_cancer = MA_normal <- t(ProcessedData$MA_TCGA)
CNV_matrix <- t(ProcessedData$CNV_TCGA)

# list of transcription factors
TFs <- read.table(paste0(path,"Data/expression/AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1 

# set parameters
VarMax <- 0.75 # filter for genes, keep only the top VarMax % variant genes

# Run the main code (using the Lasso for reconstucting the network)
Results <- LIONS_Main_Code(MA_cancer = MA_cancer,MA_normal = MA_normal,TFs = TFs,
                            VarMax = VarMax,Method="Lasso",
                           LassoThresholds=list(subsamples = 100),
                            TargetDirectory = TargetDirectory,pathEM = pathEM)

# Run the main code (using hLicorn for reconstucting the network)
Results <- LIONS_Main_Code(MA_cancer = t(ProcessedData$MA_TCGA), MA_normal = t(ProcessedData$MA_TCGA), TFs = TFs,
                           VarMax = VarMax, Method="hLICORN",CNV_matrix = CNV_matrix,CNV_correction = TRUE,
                           LicornThresholds=list(minCoregSupport=0.25,searchThresh=0.25),
                           TargetDirectory = TargetDirectory, pathEM = pathEM)