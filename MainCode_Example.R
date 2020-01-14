############## Main code for reproducing the results of paper #########################
### Identification of deregulated TFs involved in specific bladder cancer subtypes ####

# load the libraries
library(stringr)
library(CoRegNet)
library(doParallel)
library(igraph)
library(matrixStats)

# Define your path and target directory
TargetDirectory <- "results/"
pathEM <- "algoEM/"
DataDirectory <- "data_TCGA/"
source("MainCode.R") 

# load the TCGA data
load(paste0(DataDirectory,"ProcessedData_BLCA.Rdata"))
MA <- t(ProcessedData$MA_TCGA)
CNV <- t(ProcessedData$CNV_TCGA)

# load the different subtypes
Subtypes <- read.table(paste0(DataDirectory,"Classif_TCGA_mars_2018.csv"),sep=",",header=TRUE)
samples <- as.character(Subtypes$ID)
samples <- substr(samples,1,nchar(samples)-1)
Subtypes <- Subtypes[,4]
Subtypes <- as.character(Subtypes)
names(Subtypes) <- samples
Subtypes <- Subtypes[rownames(MA)]

# list of transcription factors
TFs <- read.table(paste0(DataDirectory,"AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1 

# set parameters
VarMax <- 0.75 # filter for genes, keep only the top VarMax % variant genes

# Run the main code 
subtypes <- names(which(Subtypes==unique(Subtypes)[1])) # we work on subtype 1
Results <- DeregGenes(MA_cancer = MA[subtypes,],
                           MA_normal = MA[!(rownames(MA)%in%subtypes),],
                           TFs = TFs,
                           CNV_matrix = CNV,CNV_correction = TRUE,
                           LicornThresholds=list(minCoregSupport=0.5,searchThresh=0.5),
                           TargetDirectory = TargetDirectory, pathEM = pathEM)