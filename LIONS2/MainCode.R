############## here is the main code for LIONS ##############

LIONS_Main_Code <- function(MA_cancer,MA_normal,TFs,
                            CNV_matrix, CNV_correction = FALSE,
                            VarMax=1,
                            LicornThresholds=list(minCoregSupport=0.5,searchThresh=0.5),
                            TargetDirectory,pathEM){
##### List of inputs ####
# Data:
#   - MA_cancer: matrix of gene expression from tumor samples (genes in columns, samples in rows)
#   - MA_normal: matrix of gene expression from normal tissues (genes in columns, samples in rows)
#   - TFs: list of transcription factors
#   - CNV_matrix: matrix containing CNV data (genes in coumns, samples in rows)
#
# CNV_correction: TRUE/FALSE depending on whether you want to correct gene expression data from copy number effects (if TRUE, the program needs CNV data)
#
# Parameters:
#   - VarMax: keep the top VarMax genes (by default = 1)
#   - minCoregSupport: (internal parameter of licorn - step 1) threshold for defining groups of co-regulated gene (by default = 0.5)
#   - searchThresh: (internal parameter of licorn - step 2) threshold for defining a signifcant regulation between co-regulators and their targets (by default = 0.75)
#
# Others:
#   - Target directory: where to store the outputs and data
#   - pathEM: where the EM algorithm is located

if (is.null(LicornThresholds$minCoregSupport)){
  minCoregSupport <- 0.5
} else {
  minCoregSupport <- LicornThresholds$minCoregSupport 
}
if (is.null(LicornThresholds$searchThresh)){
  searchThresh <- 0.5
} else {
  searchThresh <- LicornThresholds$searchThresh 
}

#######################################
##### 1st step: Cleaning the data #####
#######################################
  
Genes <- intersect(colnames(MA_cancer),colnames(MA_normal))
if (length(Genes)==0){
  stop("\n\t There is no genes in commun between MA_cancer and MA_normal. Please check your data.")
}
if (length(intersect(Genes,TFs))==0){
  stop("\n\t There is no TFs in the gene expression data matrix. Please check your data.")
}

cat("\n\t Restricting the number of genes to the",VarMax*100,"% with highest variance.")
# from the normal data
vars <- apply(MA_normal, 2, var)
MA_normal <- MA_normal[,order(vars, decreasing=TRUE)[1:floor(VarMax*ncol(MA_normal))]]
# from the cancer data
vars <- apply(MA_cancer, 2, var)
MA_cancer <- MA_cancer[,order(vars, decreasing=TRUE)[1:floor(VarMax*ncol(MA_cancer))]]

Genes <- intersect(colnames(MA_cancer),colnames(MA_normal))
if (length(Genes)==0){
  stop("\n\t There is no gene in commun between MA_cancer and MA_normal after thresholding the data. Please re-run the algorithm with a less contrained VarMax parameter.")
}
TFs <- intersect(Genes,TFs)
Targets <- Genes[is.na(match(Genes,TFs))]
if (length(TFs)==0){
  stop("\n\t There is no TFs with high variance. Please re-run the algorithm with a less constrained VarMax parameter.")
}
MA_cancer <- MA_cancer[,Genes]
MA_normal <- MA_normal[,Genes]
cat("\n\t Summary:",length(Genes)," genes with expression data from ",nrow(MA_cancer)," tumor samples and ",nrow(MA_normal)," normal tissues.")

if (CNV_correction==TRUE){
  Genes <- intersect(colnames(CNV_matrix),Genes)
  Samples_cancer <- intersect(rownames(CNV_matrix),rownames(MA_cancer))
  Samples_normal <- intersect(rownames(CNV_matrix),rownames(MA_normal))
  
  cat(paste0("There is ",length(Genes)," genes, ",length(Samples_cancer)," tumor samples and ",length(Samples_normal)," normal samples with both gene expression and CNV data."))
  MA_cancer <- MA_cancer[Samples_cancer,Genes]
  MA_normal <- MA_normal[Samples_normal,Genes]
  CNV_normal <- CNV_matrix[Samples_normal,Genes]
  CNV_cancer <- CNV_matrix[Samples_cancer,Genes]
  
  TFs <- intersect(Genes,TFs)
  
  # correction of the normal expression data
  MA_normal <- DataCorrection(MA = MA_normal,CNV_matrix = CNV_normal)
}

# standardize data
MA_normal = MA_normal - matrix(1,nrow(MA_normal),1) %*% colMeans(MA_normal)
MA_normal = MA_normal / sqrt( (1/nrow(MA_normal)) * matrix(1,nrow(MA_normal),1) %*% colSums(MA_normal^2))

##########################################################################
##### 2nd step: Inferring the regulatory network from normal tissues #####
##########################################################################
cat("\n\t Inferring the gene regulatory network running LICORN on normal tissues.")

# connect TFs to their targets
dummyNet <- hLICORN(numericalExpression = t(MA_normal), TFlist = TFs,
                      minCoregSupport = minCoregSupport, searchThresh = searchThresh) # this is a coregnet class object
GRNnetwork <- Transform_Network(dummyNet)

# connect TFs to other TFs
GRNnetwork2 <- c()
for (i in (1:length(TFs))){
  # Check for an error (this means that there is no interaction)
  possibleError <- tryCatch({
    dummyNet <- hLICORN(numericalExpression = t(MA_normal[,TFs]),TFlist = TFs[-i],
                      minCoregSupport = minCoregSupport, searchThresh = searchThresh)
    GRNnetwork2 <- rbind(GRNnetwork2,Transform_Network(dummyNet))
  }
  ,
  error=function(e) {
    e
    GRNnetwork2 <- GRNnetwork2
  })
}

# collapse networks
GRNnetwork2[,1] <- str_replace_all(GRNnetwork2[,1],as.character(GRNnetwork2[,1]),paste0(GRNnetwork2[,1],"_TF"))
GRNnetwork <- rbind(GRNnetwork,GRNnetwork2)

# create the adjacency matrix
GRN <- rbind(as.matrix(GRNnetwork[, c(1,2)]),
             as.matrix(GRNnetwork[, c(1,3)]))
GRN <- GRN[GRN[, 2] != "", ]
rep <- sapply(strsplit(GRN[, 2], " "), length)
GRN <- cbind(rep(GRN[, 1], rep), unlist(strsplit(GRN[, 2], " ")))
Adj_matrix <- as_adjacency_matrix(graph_from_edgelist(GRN, directed = TRUE),sparse=FALSE)
Adj_matrix <- Adj_matrix[!(rowSums(Adj_matrix) == 0), ]
Adj_matrix <- Adj_matrix[,!(colSums(Adj_matrix) == 0)]
Adj_matrix <- Adj_matrix[order(rownames(Adj_matrix)), ]
cat("\n\t The inferred network is made of",nrow(GRN),"edges, which connect",length(unique(GRN[,2])),"TFs to",length(unique(GRN[,1])),"target genes.")  
cat("\n\t Each target gene is associated with an averaged number of",mean(rowSums(Adj_matrix)),"TF. Conversely, a TF is associated with an averaged number of",mean(colSums(Adj_matrix)),"target genes.")

# update cancer data
MA_cancer_TFs <- matrix(0,ncol=(nrow(Adj_matrix)+ncol(Adj_matrix)),nrow=nrow(MA_cancer))
rownames(MA_cancer_TFs) = rownames(MA_cancer)
colnames(MA_cancer_TFs) <- c(colnames(Adj_matrix),rownames(Adj_matrix))
InterTFs <- colnames(MA_cancer_TFs)[grep("_TF",colnames(MA_cancer_TFs))] # TFs that are also targets
InterTargets <- intersect(colnames(MA_cancer),colnames(MA_cancer_TFs)) # all other genes
MA_cancer_TFs[,InterTFs] <- as.matrix(MA_cancer[,str_replace_all(InterTFs,"_TF","")],nrow(MA_cancer),length(InterTFs))
MA_cancer_TFs[,InterTargets] <- as.matrix(MA_cancer[,InterTargets],nrow(MA_cancer),length(InterTargets)) 

if (CNV_correction == TRUE){
  # we only correct the expression of target genes
  TargetGenes <- rownames(Adj_matrix)
  TargetGenes_fullname <- str_replace_all(TargetGenes,"_TF","")
  MA_Targets <- DataCorrection(MA = MA_cancer_TFs[,TargetGenes],CNV_matrix = CNV_cancer[,TargetGenes_fullname])
  MA_cancer_TFs[,TargetGenes] <- as.matrix(MA_Targets,nrow=nrow(MA_Targets),ncol=ncol(MA_Targets))
}

# standardize data
MA_cancer_TFs = MA_cancer_TFs - matrix(1,nrow(MA_cancer_TFs),1) %*% colMeans(MA_cancer_TFs)
MA_cancer_TFs = MA_cancer_TFs / sqrt( (1/nrow(MA_cancer_TFs)) * matrix(1,nrow(MA_cancer_TFs),1) %*% colSums(MA_cancer_TFs^2))

# save data
write.table(MA_cancer_TFs,file=paste0(TargetDirectory,"Expression.txt"),sep=" ",quote=FALSE)
write.table(GRNnetwork,file=paste0(TargetDirectory,"Network.txt"),sep=" ; ",quote=FALSE,col.names=colnames(GRNnetwork),row.names=FALSE)

#############################################################################
### 3rd step: computing a deregulation score for all genes in each sample ###
#############################################################################
cat("\n\t Computing a deregulation score for all genes in each sample through an EM algorithm.")
system(paste0("java -jar ",pathEM,"ddt.jar -network ",TargetDirectory,"Network.txt -expression ",TargetDirectory,"Expression.txt -scores ",TargetDirectory,"Scores.txt"))
Score <- read.table(paste0(TargetDirectory,"Scores.txt"),sep=',')
colnames(Score) <- rownames(MA_cancer_TFs)

#############################################
### 4th step: Finding the deregulated TFs ###
#############################################
cat("\n\t Identifying the deregulated TFs through a penalized linear model.")

Score <- Score[!(rowSums(Score) == 0), ]
Score <- Score[order(rownames(Score)), ]
Score <- as.matrix(Score)
#Scores_log <- log(as.matrix(Score))
#Scores_log <- log(as.matrix(Score/(1-Score)))
#Adj_matrix_adj <- t(t(Adj_matrix)/colSums(Adj_matrix))

OptimizationLoop <- function(i){
  fit <- glmnet(Adj_matrix,Score[,i],lambda=0,upper.limits=1,lower.limits=0)
  #fitcv <- cv.glmnet(Adj_matrix,Score[,i],upper.limits=1,lower.limits=0)
  #fit <- glmnet(Adj_matrix,Score[,i],lambda=fitcv$lambda.1se,upper.limits=1,lower.limits=0)
  beta <- as.vector(fit$beta)
  names(beta) <- colnames(Adj_matrix)
  beta
}
Beta <- matrix(unlist(mclapply(X=1:ncol(Score), FUN=OptimizationLoop)), ncol = ncol(Score), byrow = FALSE)
rownames(Beta) <- colnames(Adj_matrix)
colnames(Beta) <- colnames(Score)

# save final results
Results <- list(Beta=Beta,Score=Scores_log,Adj_matrix=Adj_matrix)
return(Results)
}

# this code transforms a coregnet object into a matrix that can be used as an input of the EM algorithm 
# (matrix with three columns: gene, activator, inhibitor)
Transform_Network <- function(dummyNet){
  if (length(dummyNet@adjacencyList$bygene)>1){
    Network <- foreach(j=1:length(dummyNet@adjacencyList$bygene),.combine='rbind') %do% {
      gene = dummyNet@adjacencyList$bygene[[j]]
      Acts = paste(unlist(gene$act),collapse=' ')  
      Inhs = paste(unlist(gene$rep),collapse=' ') 
      c(names(dummyNet@adjacencyList$bygene)[j],Acts,Inhs)
    }
  } else {
    gene = dummyNet@adjacencyList$bygene[[1]]
    Acts = paste(unlist(gene$act),collapse=' ')  
    Inhs = paste(unlist(gene$rep),collapse=' ') 
    Network <- matrix(c(names(dummyNet@adjacencyList$bygene)[1],Acts,Inhs),ncol=3,nrow=1)
  }
  colnames(Network) <- c("gene","activator","inhibitor")
  Network <- data.frame(Network)
  rownames(Network) <- NULL
  return(Network)
}

## this function is useful to correct gene expression for CNV
DataCorrection <- function(MA,CNV_matrix){
  cat("Correcting gene expression using CNV.\n")
  CorrectedExpression <- foreach(j=1:ncol(CNV_matrix),.combine='cbind') %do% {
    Predictions <-resid(lm(MA[,j]~as.numeric(CNV_matrix[,j])))
  }
  rownames(CorrectedExpression) <- rownames(MA)
  colnames(CorrectedExpression) <- colnames(MA)
  return(CorrectedExpression)
}
