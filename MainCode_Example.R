# Main code example

# load the libraries
library(stringr)
library(CoRegNet)
library(doParallel)
library(igraph)
library(lsei)

# Define your path and target directory
path <- "/Users/Magali/Desktop/recherche/LIONS/"
pathEM <- paste0(path,"algoEM/")
TargetDirectory <- "/Users/Magali/Desktop/recherche/LIONS/results/"
FileResults <- "Results_1sttest"
source(paste0(path,"LIONS_project/MainCode.R"))

# load the data
MA_normal <- read.table(paste0(path,"Data/expression/NormalExpressionData_AllGenes.txt"))
MA_cancer <- read.table(paste0(path,"Data/expression/TumorExpressionData_AllGenes.txt"))
rownames(MA_cancer) <- str_replace_all(rownames(MA_cancer),"[A-Z]","") # anonymous data
TFs <- read.table(paste0(path,"Data/expression/AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1 

# set parameters
VarMax <- 0.1 # keep the top 75% variant genes
LicornThresholds <- list(searchThresh=0.75,minCoregSupport = 0.5) # Licorn parameters, to check the sparsity

# Run the main code
Results <- LIONS_Main_Code(MA_cancer,MA_normal,TFs,
                            VarMax,LicornThresholds=list(searchThresh,minCoregSupport),
                            TargetDirectory,pathEM,FileResults)