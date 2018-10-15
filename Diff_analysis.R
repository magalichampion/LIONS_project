# code for differential analysis (test)
path <- "/Users/Magali/Desktop/recherche/LIONS/"
pathEM <- paste0(path,"algoEM/")
TargetDirectory <- paste0(path,"results/")
library(limma)

# load the data
load(paste0(DataDirectory,"/ProcessedData_BLCA.Rdata"))
MA_cancer <- ProcessedData$MA_TCGA

TFs <- read.table(paste0(path,"Data/expression/AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1 

Subgroups_full <- read.table(paste0(path,"Data/Subtypes/BLCA-TP.mergedcluster.csv"),sep=",",header=TRUE)
rownames(Subgroups_full) <- Subgroups_full$SampleName  
Subgroups <- Subgroups_full$mRNAseq_cHierarchical
names(Subgroups) <- rownames(Subgroups_full)
Subgroups <- Subgroups[colnames(MA_cancer)]
Subgroups <- as.character(Subgroups)
names(Subgroups) <- colnames(MA_cancer)

# run the main code
# create the design matrix
design <- cbind(ST1=1*(Subgroups=="1"),ST2=1*(Subgroups=="2"),ST3=1*(Subgroups=="3"))

# fit the model
fit <- lmFit(MA_cancer,design)

# create the contrast matrix
Contrasts <- makeContrasts(ST1-ST2,ST1-ST3,ST2-ST3,levels=design)
fit2 <- contrasts.fit(fit,Contrasts)

# run differential analysis
Fit <- eBayes(fit2)

# see the results
topTable(Fit,adjust="BH",p.value = 0.05,number=nrow(MA_cancer))
Fit$p.value

# ST1 vs ST2
gene <- rownames(topTable(Fit,adjust="BH",p.value = 0.05,number=nrow(MA_cancer),coef=1))
length(intersect(rownames(Beta),gene))
