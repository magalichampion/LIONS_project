############## here is the main code for LIONS ##############

LIONS_Main_Code <- function(MA_cancer,MA_normal,TFs,
                            VarMax,LicornThresholds=list(searchThresh,minCoregSupport),
                            TargetDirectory,pathEM,FileResults){
##### List of inputs ####
# Data:
#   - MA_cancer: matrix of gene expression from tumor samples (genes in columns, samples in rows)
#   - MA_normal: matrix of gene expression from normal tissues (genes in columns, samples in rows)
#   - TFs: list of transcription factors
#
# Parameters:
#   - VarMax: keep the top VarMax genes (by default = 1)
#   - LicornThresholds: list of Licorn internal parameters searchThresh and minCoregSupport
#
# Others:
#   - Target directory: where to store the outputs and data
#   - pathEM: where the EM algorithm is located
#   - FileResults: give a name to the results

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

cat("\n\t Restrict the number of genes to the",VarMax*100,"% with highest variance.")
# from the normal data
vars <- apply(MA_normal, 2, var)
MA_normal <- MA_normal[,order(vars, decreasing=TRUE)[1:floor(VarMax*ncol(MA_normal))]]
# from the cancer data
vars <- apply(MA_cancer, 2, var)
MA_cancer <- MA_cancer[,order(vars, decreasing=TRUE)[1:floor(VarMax*ncol(MA_cancer))]]

Genes <- intersect(colnames(MA_cancer),colnames(MA_normal))
if (length(Genes)==0){
  stop("\n\t There is no genes in commun between MA_cancer and MA_normal. Please re-run the algorithm with a less contrained VarMax parameter.")
}
TFs <- intersect(Genes,TFs)
if (length(TFs)==0){
  stop("\n\t There is no TFs with high variance. Please re-run the algorithm with a less constrained VarMax parameter.")
}
MA_cancer <- MA_cancer[,Genes]
MA_normal <- MA_normal[,Genes]
cat("\n\t Summary:",length(Genes),"genes with expression data from",nrow(MA_cancer),"tumor samples and",nrow(MA_normal),"normal tissues.")

# standardize data
MA_normal = MA_normal - matrix(1,nrow(MA_normal),1) %*% colMeans(MA_normal)
MA_normal = MA_normal / sqrt( (1/nrow(MA_normal)) * matrix(1,nrow(MA_normal),1) %*% colSums(MA_normal^2))

##########################################################################
##### 2nd step: Inferring the regulatory network from normal tissues #####
##########################################################################
cat("\n\t Infer the gene regulatory network running hLICORN on normal tissues.")

# connect TFs to their targets
dummyNet <- hLICORN(numericalExpression = t(MA_normal),TFlist = TFs,minCoregSupport = LicornThresholds$minCoregSupport,
                    searchThresh = LicornThresholds$searchThresh) # this is a coregnet class object
GRNnetwork <- Transform_Network(dummyNet)

# connect TFs to other TFs
GRNnetwork2 <- c()
for (i in (1:length(TFs))){
  # Check for an error (this means that there is no interaction)
  possibleError <- tryCatch({
    dummyNet <- hLICORN(numericalExpression = t(MA_normal[,TFs]),TFlist = TFs[-i],minCoregSupport = LicornThresholds$minCoregSupport,
                          searchThresh = LicornThresholds$searchThresh)
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
 
# update cancer data
MA_cancer_TFs <- matrix(0,ncol=(nrow(Adj_matrix)+ncol(Adj_matrix)),nrow=nrow(MA_cancer))
rownames(MA_cancer_TFs) = rownames(MA_cancer)
colnames(MA_cancer_TFs) <- c(colnames(Adj_matrix),rownames(Adj_matrix))
InterTargets <- intersect(colnames(MA_cancer),colnames(MA_cancer_TFs)) 
InterTFs <- colnames(MA_cancer_TFs)[grep("_TF",colnames(MA_cancer_TFs))]
MA_cancer_TFs[,InterTargets] <- as.matrix(MA_cancer[,InterTargets],nrow(MA_cancer),length(InterTargets)) 
MA_cancer_TFs[,InterTFs] <- as.matrix(MA_cancer[,str_replace_all(InterTFs,"_TF","")],nrow(MA_cancer),length(InterTFs))

# standardize data
MA_cancer_TFs = MA_cancer_TFs - matrix(1,nrow(MA_cancer_TFs),1) %*% colMeans(MA_cancer_TFs)
MA_cancer_TFs = MA_cancer_TFs / sqrt( (1/nrow(MA_cancer_TFs)) * matrix(1,nrow(MA_cancer_TFs),1) %*% colSums(MA_cancer_TFs^2))

# save data
write.table(MA_cancer_TFs,file=paste0(TargetDirectory,"Expression.txt"),sep=" ",quote=FALSE)
write.table(GRNnetwork,file=paste0(TargetDirectory,"Network.txt"),sep=" ; ",quote=FALSE,col.names=colnames(GRNnetwork),row.names=FALSE)

#############################################################################
### 3rd step: computing a deregulation score for all genes in each sample ###
#############################################################################
cat("\n\t Compute a deregulation score for all genes in each sample trhough an EM algorithm.")
system(paste0("java -jar ",pathEM,"ddt.jar -network ",TargetDirectory,"Network.txt -expression ",TargetDirectory,"Expression.txt -scores ",TargetDirectory,"Scores.txt"))
Score <- read.table(paste0(TargetDirectory,"Scores.txt"),sep=',')
colnames(Score) <- substring(colnames(Score),2)

#############################################
### 4th step: Finding the deregulated TFs ###
#############################################
cat("\n\t Identify the deregulated TFs through a penalized linear model.")
Score <- Score[!(rowSums(Score) == 0), ]
Score <- Score[order(rownames(Score)), ]
Scores_log <- scale(log10(as.matrix(Score)))

OptimizationLoop <- function(i){
  lsei(a=Adj_matrix,b=Scores_log[,i],
       c=rep(1,ncol(Adj_matrix)),d=1,
       e=diag(ncol(Adj_matrix)),f=rep(0,ncol(Adj_matrix)))
}
Beta <- matrix(unlist(mclapply(X=1:ncol(Scores_log), FUN=OptimizationLoop)), ncol = ncol(Scores_log), byrow = FALSE)
rownames(Beta) <- colnames(Adj_matrix)
colnames(Beta) <- colnames(Scores_log)

# save final results
Results <- list(Beta=Beta,Score=Score,Adj_matrix=Adj_matrix)
save(Results,file=c(paste0(TargetDirectory,FileResults,".Rdata")))
}

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
