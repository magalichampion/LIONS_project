############## here is the main code for LIONS ##############

LIONS_Main_Code <- function(MA_cancer,MA_normal,TFs,
                            CNV_matrix, CNV_correction = FALSE,
                            Method="hLICORN",
                            VarMax=1,
                            LicornThresholds=list(minCoregSupport=0.5,searchThresh=0.5),
                            LassoThresholds=list(subsamples=ncol(MA_normal)/2,weakness=0.5,maxScore=0.75,n.lambda=50,min.ratio=1e-4),
                            TreeCoopLassoThresholds=list(subsamples=100,weakness=0.5,n.lambda=50,min.ratio=1e-4),
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
# Methods: to infer a GRN, you can use one of the following method
#   - hLICORN: implemented in the R package Coregnet (by default)
#   - Lasso + stability selection
#   - Cooperative Lasso: implemented in the R package scoop (not yet)
#
# Parameters:
#   - VarMax: keep the top VarMax genes (by default = 1)
#
#   For hLICORN: LicornThresholds (internal parameters of Licorn)
#   - minCoregSupport: step 1 of Licorn, threshold for defining groups of co-regulated gene (by default = 0.5)
#   - searchThresh: step 2 of Licorn, threshold for defining a signifcant regulation between co-regulators and their targets (by default = 0.75)
#   
#   For Lasso: Lasso thresholds (optimized version by J. Chiquet)
#   - Subsamples: number of subsamples for the sability selection step - should be large but may take some time to run (by default = n/2)
#   - weakness: (by default = 0.5)
#   - maxScore: selection of edges with a probability larger than maxScore (by default = 0.75)
#   - n.lambda: number of penalities (by default = 50)
#   - min.ratio: minimal penalty as a ratio of the largest (by default = 1e-4)
#   - mc.cores: the number of cores (by default = 1)
#  
# Others:
#   - Target directory: where to store the outputs and data
#   - pathEM: where the EM algorithm is located

if (Method=="hLICORN"){
   if (is.null(LicornThresholds$minCoregSupport)){
    minCoregSupport <- 0.5
  } else {
    minCoregSupport <- LicornThresholds$minCoregSupport 
  }
  if (is.null(LicornThresholds$searchThresh)){
    searchThresh <- 0.75
  } else {
    searchThresh <- LicornThresholds$searchThresh 
  }
}

if (Method == "Lasso"){
  if (is.null(LassoThresholds$subsamples)){
    subsamples <- ncol(MA_normal)/2
  } else {
    subsamples <- LassoThresholds$subsamples
  }
  if (is.null(LassoThresholds$weakness)){
    weakness <- 0.5
  } else {
    weakness <- LassoThresholds$weakness
  }
  if (is.null(LassoThresholds$maxScore)){
    maxScore <- 0.75
  } else {
    maxScore <- LassoThresholds$maxScore
  }
  if (is.null(LassoThresholds$n.lambda)){
    n.lambda <- 50
  } else {
    n.lambda <- LassoThresholds$n.lambda
  }
  if (is.null(LassoThresholds$min.ratio)){
    min.ratio <- 1e-4
  } else {
    min.ratio <- LassoThresholds$min.ratio
  }
  if (is.null(LassoThresholds$mc.cores)){
    mc.cores <- 1
  } else {
    mc.cores <- LassoThresholds$mc.cores
  }
}

if (Method == "Tree-CoopLasso"){
  if (is.null(TreeCoopLasso$subsamples)){
    subsamples <- 100
  } else {
    subsamples <- TreeCoopLasso$subsamples
  }
  if (is.null(TreeCoopLasso$weakness)){
    weakness <- 0.5
  } else {
    weakness <- TreeCoopLasso$weakness
  }
  if (is.null(TreeCoopLasso$n.lambda)){
    n.lambda <- 50
  } else {
    n.lambda <- TreeCoopLasso$n.lambda
  }
  if (is.null(TreeCoopLasso$min.ratio)){
    min.ratio <- 1e-4
  } else {
    min.ratio <- TreeCoopLasso$min.ratio
  }
  if (is.null(TreeCoopLasso$mc.cores)){
    mc.cores <- 1
  } else {
    mc.cores <- TreeCoopLasso$mc.cores
  }
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
MA_cancer <- MA_cancer[,Genes]
MA_normal <- MA_normal[,Genes]

if (CNV_correction==TRUE){
  Genes <- intersect(colnames(CNV_matrix),Genes)
  cat(paste0("There is ",length(Genes)," genes with both gene expression and CNV data."))
  MA_cancer <- MA_cancer[,Genes]
  MA_normal <- MA_normal[,Genes]
  CNV_matrix <- CNV_matrix[,Genes]
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
  stop("\n\t There is no genes in commun between MA_cancer and MA_normal. Please re-run the algorithm with a less contrained VarMax parameter.")
}
TFs <- intersect(Genes,TFs)
Targets <- Genes[is.na(match(Genes,TFs))]
if (length(TFs)==0){
  stop("\n\t There is no TFs with high variance. Please re-run the algorithm with a less constrained VarMax parameter.")
}
MA_cancer <- MA_cancer[,Genes]
MA_normal <- MA_normal[,Genes]
cat("\n\t Summary:",length(Genes),"genes with expression data from",nrow(MA_cancer),"tumor samples and",nrow(MA_normal),"normal tissues.")

if (CNV_correction == TRUE){ # correction for the target genes only
  MA_normal_corrected <- DataCorrection(MA = MA_normal[,Targets],CNV_matrix = CNV_matrix)
  MA_normal[,Targets] <- MA_normal_corrected
  
  MA_cancer_corrected <- DataCorrection(MA = MA_cancer[,Targets],CNV_matrix = CNV_matrix)
  MA_cancer[,Targets] <- MA_cancer_corrected
}
#ProcessedData <- list(MA_matrix= MA_matrix,MA_matrix_corrected=MA_matrix_corrected,CNV_matrix=CNV_matrix,TFs=TFs)

# standardize data
MA_normal = MA_normal - matrix(1,nrow(MA_normal),1) %*% colMeans(MA_normal)
MA_normal = MA_normal / sqrt( (1/nrow(MA_normal)) * matrix(1,nrow(MA_normal),1) %*% colSums(MA_normal^2))

##########################################################################
##### 2nd step: Inferring the regulatory network from normal tissues #####
##########################################################################
cat(paste0("\n\t Inferring the gene regulatory network running ",Method," on normal tissues."))

if (Method == "hLICORN"){
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
}

if (Method == "Lasso"){
  # connect TFs to their targets
  GRNnetwork <- Lasso.bigraph(numericalExpression = MA_normal,TFlist = TFs,
                              subsamples,weakness,maxScore,min.ratio,n.lambda,mc.cores)

  # connect TFs to other TFs
  GRNnetwork2 <- c()
  for (i in (1:length(TFs))){
    GRNnetwork2 <- rbind(GRNnetwork2,Lasso.bigraph(numericalExpression = MA_normal[,TFs],TFlist = TFs[-i],
                               subsamples,weakness,maxScore,min.ratio,n.lambda,mc.cores))
  }
  
  # collapse networks
  GRNnetwork2[,1] <- str_replace_all(GRNnetwork2[,1],as.character(GRNnetwork2[,1]),paste0(GRNnetwork2[,1],"_TF"))
  GRNnetwork <- rbind(GRNnetwork,GRNnetwork2)
}

if (Method == "Scoop"){
  # data
  cat("Not available yet.")
  
}

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
InterTargets <- intersect(colnames(MA_cancer),colnames(MA_cancer_TFs)) 
InterTFs <- colnames(MA_cancer_TFs)[grep("_TF",colnames(MA_cancer_TFs))]
MA_cancer_TFs[,InterTargets] <- as.matrix(MA_cancer[,InterTargets],nrow(MA_cancer),length(InterTargets)) 
MA_cancer_TFs[,InterTFs] <- as.matrix(MA_cancer[,str_replace_all(InterTFs,"_TF","")],nrow(MA_cancer),length(InterTFs))

if (CNV_correction == TRUE){
  MA_cancer_TFs <- DataCorrection(MA_matrix = MA_cancer_TFs,CNV_matrix = CNV_matrix)
}
ProcessedData <- list(MA_matrix= MA_matrix,MA_matrix_corrected=MA_matrix_corrected,CNV_matrix=CNV_matrix,TFs=TFs)

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
Scores_log <- as.matrix(Score)
#Scores_log <- scale(log10(as.matrix(Score)))

OptimizationLoop <- function(i){
  a = Adj_matrix
  b <- Scores_log[,i]
  e <- rbind(diag(ncol(Adj_matrix)),-diag(ncol(Adj_matrix)))
  f <- c(rep(0,ncol(Adj_matrix)),rep(-1,ncol(Adj_matrix)))
  lsei(a=a,b=b,
#       c=rep(1,ncol(Adj_matrix)),d=1,
       e=e,f=f)
}
Beta <- matrix(unlist(mclapply(X=1:ncol(Scores_log), FUN=OptimizationLoop)), ncol = ncol(Scores_log), byrow = FALSE)
rownames(Beta) <- colnames(Adj_matrix)
colnames(Beta) <- colnames(Scores_log)
Beta[which(Beta<0)] <- rep(0,length(which(Beta<0)))

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

stab.lasso  <- function(X, Y, lambda, subsamples, weakness, mc.cores) {
  # define the parameters for resamplings
  n <- nrow(X)
  p <- ncol(X)
  folds <- replicate(subsamples, sample(1:n, floor(n/2)), simplify=FALSE)
  blocs <- split(1:subsamples, 1:mc.cores)
  
  # function to run on each core
  bloc.stability <- function(subsets) {
    return(Reduce("+", lapply(subsets, function(s) {
      weights <- runif(p,weakness,1)
      lars.out <- lars(scale(X, FALSE, 1/weights), Y, intercept=FALSE, normalize=FALSE)
      return(predict(lars.out, s=lambda, type="coefficients", mode="lambda")$coefficients != 0)
    }))/(length(subsets)*length(blocs)))
  }
  return(Reduce("+", mclapply(blocs, bloc.stability, mc.cores=mc.cores)))
}

Lasso.bigraph <- function(numericalExpression,TFlist,
                          subsamples,weakness,maxScore,min.ratio,n.lambda,mc.cores){
  # define the matrices
  MA_TFs <- numericalExpression[,TFlist]
  if (length(intersect(colnames(numericalExpression),TFlist))==(ncol(numericalExpression)-1)){
    MA_Targets <- as.matrix(numericalExpression[,-which(colnames(numericalExpression) %in% TFlist)],ncol=1,nrow=nrow(numericalExpression))
    colnames(MA_Targets) <-  colnames(numericalExpression)[-which(colnames(numericalExpression)%in%TFlist)]
    rownames(MA_Targets) <- rownames(numericalExpression)
  } else {
    MA_Targets <- numericalExpression[,-which(colnames(numericalExpression) %in% TFlist)]
  }
  
  # set the grid of lambda
  lmax <- max(apply(abs(crossprod(as.matrix(MA_TFs),as.matrix(MA_Targets))), 1, max))
  lambda <- 10^seq(log10(lmax), log10(min.ratio * lmax), len = n.lambda)

  # compute the stability path for each target gene
  if (ncol(MA_Targets)>1) pb <- txtProgressBar(min = 0, max = ncol(MA_Targets), style = 3)
  all.stab.path <- lapply(1:ncol(MA_Targets), function(i) {
    if (ncol(MA_Targets)>1){
      setTxtProgressBar(pb,i)
    }
    stab.lasso(X=MA_TFs, Y = MA_Targets[, i],lambda =  lambda,subsamples =  subsamples,weakness =  weakness, mc.cores=mc.cores)
  }
  )

  # compute a confidence score for edges
  score.func <- function(stab.path){
   m <- as.matrix(stab.path)
    max <- colMaxs(m)
    concerned <- which(max>0)
    return(sparseVector(max[concerned],concerned,ncol(m)))
  }
  score.list <- sapply(all.stab.path,score.func)
  score.mat <- Matrix(0,ncol(MA_TFs),ncol(MA_Targets))

  for (i in 1:length(score.list)){
    score.mat[,i] <- score.list[[i]]
  }
  colnames(score.mat) <- colnames(MA_Targets)
  rownames(score.mat) <- colnames(MA_TFs)
  ind <- which(score.mat != 0, arr.ind=TRUE)
  Score <- data.frame(TFs = colnames(MA_TFs)[ind[, 1]], Targets = colnames(MA_Targets)[ind[, 2]], weight=c(score.mat[ind]))

  # keep the best scores
  Score <- Score[Score$weight>maxScore,c(1,2)]

  # refit the model with OLS
  Targets <- unique(Score$Targets)
  predictors <- lapply(lapply(Targets, function(target) which(Score$Targets %in% target)), function(pred.id) Score$TFs[pred.id])
  edges.refit <- vector("list",length(Targets))
  for (i in 1:length(Targets)){
    Y <- MA_Targets[,match(Targets[i],colnames(MA_Targets))]
    X <- MA_TFs[,predictors[[i]],drop=FALSE]
    X <- as.matrix(X,ncol(X),nrow=nrow(X))
    edges.refit[[i]] <- coefficients(lm.fit(X,Y))
  }

  # create the associated graph and list of edges
  g <- graph.edgelist(as.matrix(Score), directed=TRUE)
  status <- rep("target", length(V(g)$name))
  status[V(g)$name %in% Score[, 1]] <- "TF"
  status <- factor(status)
  V(g)$color <- as.numeric(status)
  E(g)$value <- unlist(edges.refit)
  E(g)$color <- sign(unlist(edges.refit)) + 3    

  dg <- get.data.frame(g)
  dg$status <- rep("activator", nrow(dg))
  dg$status[dg$value <= 0] <- "inhibitor"
  dg$status <- as.factor(dg$status)
  dg <- split(dg[, c(1,5)], dg$to)
  
  GRNnetwork <- matrix(0,ncol=3,nrow=length(dg))
  GRNnetwork[,1] <- names(dg)
  colnames(GRNnetwork) <- c("gene","activator","inhibitor")
  for (i in (1:length(dg))){
    g <- dg[[i]]
    GRNnetwork[i,2] <- paste(g$from[which(g$status=="activator")],collapse=" ")
    GRNnetwork[i,3] <- paste(g$from[which(g$status=="inhibitor")],collapse=" ")
  }
  GRNnetwork <- data.frame(GRNnetwork)
  return(GRNnetwork)
}

Scoop.lasso <- function(numericalExpression,TFlist){
  ################### to be done
  # load the data
  MA_TFs <- numericalExpression[,TFlist]
  if (length(intersect(colnames(numericalExpression),TFlist))==(ncol(numericalExpression)-1)){
    MA_Targets <- as.matrix(numericalExpression[,-which(colnames(numericalExpression) %in% TFlist)],ncol=1,nrow=nrow(numericalExpression))
    colnames(MA_Targets) <-  colnames(numericalExpression)[-which(colnames(numericalExpression)%in%TFlist)]
    rownames(MA_Targets) <- rownames(numericalExpression)
  } else {
    MA_Targets <- numericalExpression[,-which(colnames(numericalExpression) %in% TFlist)]
  }
  
  # get a comparable set of lambda
  lmax <- max(apply(abs(crossprod(as.matrix(MA_TFs),as.matrix(MA_Targets))), 1, max))
  lambda <- 10^seq(log10(lmax), log10(min.ratio * lmax), len = n.lambda)

  # 1st step: build the hierarchy
  distance <- as.dist(1-cor(MA_TFs)/2)
  TFhierarchy <- hclust(distance,method = "average")
}

## this function is useful to correct gene expression for CNV
DataCorrection <- function(MA,CNV_matrix){
  # 1st: we need to overlap the datasets
  OverlapSamples <- intersect(rownames(MA),rownames(CNV_matrix))
  OverlapGenes <- intersect(colnames(MA),colnames(CNV_matrix))
  cat(paste0("We have ",length(OverlapGenes)," genes and ",length(OverlapSamples)," samples with both MA and CNV data.\n"))
  
  CNV_matrix <- CNV_matrix[OverlapSamples,OverlapGenes]
  MA_matrix_reduced <- MA[OverlapSamples,OverlapGenes]
  
  # 2nd: correction using CNV
  cat("Correcting gene expression using CNV.\n")
  CorrectedExpression <- foreach(j=1:ncol(CNV_matrix),.combine='cbind') %do% {
    Predictions <-resid(lm(MA_matrix_reduced[,j]~as.numeric(CNV_matrix[,j]),na.action=na.exclude))
  }
  rownames(CorrectedExpression) <- rownames(MA_matrix_reduced)
  CorrectedExpression[is.na(CorrectedExpression)] <- MA_matrix_reduced[is.na(CorrectedExpression)]
  MA_matrix_corrected = MA
  MA_matrix_corrected[OverlapSamples,OverlapGenes] <- CorrectedExpression
  #MA_matrix_corrected <- scale(MA_matrix_corrected)
  return(MA_matrix_corrected)
}