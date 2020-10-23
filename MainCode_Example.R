############# Main code for reproducing the results of paper ##############
## Identification of deregulation mechanisms specific to cancer subtypes ##

# load the libraries
library(CoRegNet)
library(doParallel)
library(stringr)
library(glmnet)

# Define your path and target directory
TargetDirectory <- "results/" # before running the algorithm, a "results" folder must be created
pathEM <- "algoEM/"
DataDirectory <- "data_TCGA/"
source("MainCode.R") 

# load the TCGA data
load(paste0(DataDirectory,"ProcessedData_BLCA.Rdata"))
MA <- t(ProcessedData$MA_TCGA)
CNV <- t(ProcessedData$CNV_TCGA)

# load the different subtypes
load(paste0(DataDirectory,"/ConsensusSubtypes.Rdata"))
Subtypes <- subtypes
names(Subtypes) <- substr(names(Subtypes), 0, 15)
Subtypes <- Subtypes[rownames(MA)]

# list of transcription factors
TFs <- read.table(paste0(DataDirectory,"AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1 

# set parameters
VarMax <- 0.75 # filter for genes, keep only the top VarMax % variant genes

# Run the main code 
subtypes <- names(which(Subtypes==unique(Subtypes)[1])) # for the 1st subtype
Results <- DeregGenes(MA_cancer = MA[subtypes,],
                      MA_normal = MA[!(rownames(MA)%in%subtypes),],
                      TFs = TFs,
                      CNV_matrix = CNV,CNV_correction = TRUE,
                      VarMax=VarMax,
                      LicornThresholds=list(minCoregSupport=0.5,searchThresh=0.5),
                      TargetDirectory = TargetDirectory, pathEM = pathEM)
