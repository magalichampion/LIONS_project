# Main code example

# load the libraries
library(stringr)
library(CoRegNet)
library(doParallel)
library(igraph)
library(lsei)
library(lars)
library(matrixStats)
# library(scoop) # does not work yet
               # should be installed by another way
# source("/Users/Magali/Desktop/recherche/LIONS/codejulien/scoop/R/scoop.R")

# Define your path and target directory
path <- "/Users/Magali/Desktop/recherche/LIONS/"
pathEM <- paste0(path,"algoEM/")
TargetDirectory <- paste0(path,"results/")
source(paste0(path,"LIONS_project/MainCode.R"))

# load the data
# first type of data (CIT)
MA_normal <- read.table(paste0(path,"Data/expression/NormalExpressionData_AllGenes.txt"))
MA_cancer <- read.table(paste0(path,"Data/expression/TumorExpressionData_AllGenes.txt"))
rownames(MA_cancer) <- str_replace_all(rownames(MA_cancer),"[A-Z]","") # anonymous data
TFs <- read.table(paste0(path,"Data/expression/AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1 

# second type of data (TCGA)

# set parameters
VarMax <- 0.75 # keep the top VarMax% variant genes (to make it going faster)

# Run the main code (using the Lasso for reconstucting the network)
Results <- LIONS_Main_Code(MA_cancer = MA_cancer,MA_normal = MA_normal,TFs = TFs,
                            VarMax = VarMax,Method="Lasso",
                           LassoThresholds=list(subsamples = 100),
                            TargetDirectory = TargetDirectory,pathEM = pathEM)

# Run the main code (using hLicorn for reconstucting the network)
Results <- LIONS_Main_Code(MA_cancer = MA_cancer, MA_normal = MA_normal, TFs = TFs,
                           VarMax = VarMax, Method="hLICORN",
                           TargetDirectory = TargetDirectory, pathEM = pathEM)