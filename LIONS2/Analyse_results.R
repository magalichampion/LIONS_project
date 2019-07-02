### Analyse the different results

# libraries
library(stringr)
library(ggplot2)
library(igraph)

##### load the different data types #####
# paths
path <- "/Users/mchampion/Desktop/LIONS/" # Mac path
path = "/users/map5_2/mchampio/Bureau/LIONS/" # descartes path
TargetDirectory <- paste0(path,"LIONS/LIONS2/Results/")
DataDirectory <- paste0(path,"Data/TCGA")

# expression data
load(paste0(DataDirectory,"/ProcessedData_BLCA.Rdata"))
MA <- t(ProcessedData$MA_TCGA)

# CNV data
CNV <- t(ProcessedData$CNV_TCGA)
CNV_thr <- t(ProcessedData$CNV_TCGA_thresholded)
rownames(CNV_thr) <- paste0(rownames(CNV_thr),"-01")

# subtypes data
Subtypes <- read.table(paste0(path,"Data/Subtypes/Classif_TCGA_mars_2018.csv"),sep=",",header=TRUE)
samples <- as.character(Subtypes$ID)
samples <- substr(samples,1,nchar(samples)-1)
Subtypes <- Subtypes[,4]
Subtypes <- as.character(Subtypes)
names(Subtypes) <- samples

Samples <- intersect(names(Subtypes),intersect(rownames(CNV),rownames(MA)))
Genes <- intersect(colnames(CNV),colnames(MA))
Subtypes <- Subtypes[Samples]
MA <- MA[Samples,Genes]
CNV <- CNV[Samples,Genes]
CNV_thr <- CNV_thr[Samples,]

all_subtypes <-c("BasalSquamous","Luminal","LuminalInfiltrated","LuminalPapillary","Neuronal")

#### Combine all Beta ####
load(paste0(TargetDirectory,"Beta_",all_subtypes[1],".Rdata"))
Beta <- Beta[!(rowSums(Beta) == 0),]
Beta_mean <- matrix(0,nrow=nrow(Beta),ncol=length(all_subtypes))
Beta_sparse <- matrix(0,nrow=nrow(Beta),ncol=length(all_subtypes))
colnames(Beta_mean) = colnames(Beta_sparse) <- all_subtypes
Beta_mean[,1] <- rowMeans(Beta)
Beta_sparse[,1] <- rowSums(Beta>0)/ncol(Beta)
rownames(Beta_mean) = rownames(Beta_sparse) <- rownames(Beta)

for (i in (2:length(all_subtypes))){
  load(paste0(TargetDirectory,"Beta_",all_subtypes[i],".Rdata"))
  Beta <- Beta[!(rowSums(Beta) == 0),]
  Beta_t <- rowMeans(Beta)
  Beta_t2 <- rowSums(Beta>0)/ncol(Beta)
  genes <- intersect(names(Beta_t),rownames(Beta_mean))
  Beta_mean[genes,i] <- Beta_t[genes] 
  Beta_sparse[genes,i] <- Beta_t2[genes]
  
  Beta_t3 <- Beta_t[-which(names(Beta_t)%in%genes)]
  Beta_t4 <- Beta_t2[-which(names(Beta_t2)%in%genes)]
  Beta_mean2 <- matrix(0,ncol=length(all_subtypes),nrow=nrow(Beta_mean)+length(Beta_t3))
  colnames(Beta_mean2) <- colnames(Beta_mean)
  rownames(Beta_mean2) <- c(rownames(Beta_mean),names(Beta_t3))
  Beta_sparse2 <- matrix(0,ncol=length(all_subtypes),nrow=nrow(Beta_sparse)+length(Beta_t4))
  colnames(Beta_sparse2) <- colnames(Beta_sparse)
  rownames(Beta_sparse2) <- c(rownames(Beta_sparse),names(Beta_t4))
  
  Beta_mean2[rownames(Beta_mean),colnames(Beta_mean)] <- Beta_mean
  Beta_mean2[names(Beta_t3),i] <- Beta_t3
  Beta_mean <- Beta_mean2
  Beta_sparse2[rownames(Beta_sparse),colnames(Beta_sparse)] <- Beta_sparse
  Beta_sparse2[names(Beta_t4),i] <- Beta_t4
  Beta_sparse <- Beta_sparse2
}
write.table(Beta_mean,file=paste0(TargetDirectory,"Beta_Mean.csv"),sep=";")
write.table(Beta_sparse,file=paste0(TargetDirectory,"Beta_Sparse.csv"),sep=";")

#### work on one subtype ####
which_subtype <- all_subtypes[3]

# network
GRN_LicornOutput <- read.table(paste0(TargetDirectory,"Network_",which_subtype,".txt"), sep=";", strip.white=TRUE, header=TRUE)
GRN <- rbind(as.matrix(GRN_LicornOutput[, c(1,2)]),
             as.matrix(GRN_LicornOutput[, c(1,3)]))
GRN <- GRN[GRN[, 2] != "", ]
rep <- sapply(strsplit(GRN[, 2], " "), length)
GRN <- cbind(rep(GRN[, 1], rep), unlist(strsplit(GRN[, 2], " ")))
Adj_matrix <- as_adjacency_matrix(graph_from_edgelist(GRN, directed = TRUE),sparse=FALSE)
Adj_matrix <- Adj_matrix[!(rowSums(Adj_matrix) == 0), ]
Adj_matrix <- Adj_matrix[,!(colSums(Adj_matrix) == 0)]
Adj_matrix <- Adj_matrix[order(rownames(Adj_matrix)), ]
cat("\n\t The inferred network is made of",nrow(GRN),"edges, which connect",length(unique(GRN[,2])),"TFs to",length(unique(GRN[,1])),"target genes.")  
cat("\n\t Each target gene is associated with an averaged number of",mean(rowSums(Adj_matrix)),"TF. Conversely, a TF is associated with an averaged number of",mean(colSums(Adj_matrix)),"target genes.")

# scores
Score <- read.table(paste0(TargetDirectory,"Scores_",which_subtype,".txt"),sep=',')
colnames(Score) <- gsub(".", '-', colnames(Score), fixed = T)
Score <- Score[!(rowSums(Score) == 0), ]
Score <- Score[order(rownames(Score)), ]
Score <- as.matrix(Score)

# Beta
load(paste0(TargetDirectory,"Beta_",which_subtype,".Rdata"))
Beta <- Beta[!(rowSums(Beta) == 0),]

# mutation
Mutations <- read.table(paste0(path,"Data/TCGA/mutations_complet.csv"),header=TRUE,sep=";")
genes <- unique(Mutations[,2])
mutations <- matrix(0,ncol=ncol(Mutations)-1,nrow=length(genes))
colnames(mutations) <- colnames(Mutations)[-1]
rownames(mutations) <- genes
for (i in c(1:length(genes))){
  gene <- genes[i]
  I <- which(Mutations[,2]==gene)
  mutations[i,] <- mapply(Mutations[I[1],2:ncol(Mutations)],FUN=as.character)
}
mutations <- mutations[,-1]
colnames(mutations) <- gsub(".", '-', colnames(mutations), fixed = T)

##### Histograms #####
Hist <- hist(Score)
Hist$counts <- log(Hist$counts)
Hist$counts[which(Hist$counts<0)] <- rep(0,length(which(Hist$counts<0)))
plot(Hist,main="",xlab="Deregulation scores",ylab="Frequency",yaxt="n")
axis(side=2,at=c(0,2.3,4.6,6.9,9.2,11.5,13.8),labels=c("1","10","1e2","1e3","1e4","1e5","1e6"))

Hist <- hist(Beta)
Hist$counts <- log(Hist$counts)
Hist$counts[which(Hist$counts<0)] <- rep(0,length(which(Hist$counts<0)))
plot(Hist,main="",xlab="B",ylab="Frequency",yaxt="n")
axis(side=2,at=c(0,2.3,4.6,6.9,9.2,11.5,13.8),labels=c("1","10","1e2","1e3","1e4","1e5","1e6"))

#### Vérification des Beta (pour un échantillon) ####
beta <- Beta[,2]
ordre <- sort(beta,decreasing = T)
Incidence <- colSums(Adj_matrix>0)
Incidence[names(ordre)] # ok, ça marche!

##### CNV comparison #####
which_subtype <- all_subtypes[5]

Score <- read.table(paste0(TargetDirectory,"Scores_",which_subtype,".txt"),sep=',')
colnames(Score) <- gsub(".", '-', colnames(Score), fixed = T)
Score <- Score[!(rowSums(Score) == 0), ]
Score <- Score[order(rownames(Score)), ]
Score <- as.matrix(Score)

Input <- Score # Beta ou Score
Genes <- intersect(rownames(Input),colnames(CNV_thr))
Samples <- intersect(colnames(Input),rownames(CNV_thr))
CNV_t <- as.matrix(CNV_thr[Samples,Genes])
Input_CNV <- Input[Genes,Samples]

# plot the distribution
I <- which(Input_CNV<1e-2)
Input_m <- Input_CNV[-I]
CNV_m <- CNV_t[-I]
par(mfrow = c(1, 1))

F2 <- ecdf(Input_m[which(CNV_m==2)])
F1 <- ecdf(Input_m[which(CNV_m==1)])
F1m <- ecdf(Input_m[which(CNV_m==-1)])
F2m <- ecdf(Input_m[which(CNV_m==-2)])
F0 <- ecdf(Input_m[which(CNV_m==0)])

plot(F2,verticals=TRUE,do.points=FALSE,main="Cumulative empirical distribution",xlab="Score of deregulation",ylab="",xlim=c(0,1))
lines(F1,verticals=TRUE,do.points=FALSE,col="green")
lines(F1m,verticals=TRUE,do.points=FALSE,col="red")
lines(F2m,verticals=TRUE,do.points=FALSE,col="yellow")
lines(F0,verticals=TRUE,do.points=FALSE,col="blue")
legend(0.6,0.7,c("CNV=2","CNV=1","CNV=-1","CNV=-2","CNV=0"),lty=rep(1,5),col=c("black","green","red","yellow","blue"))

# Student test for testing the differences
Ttest <- rep(0,length(table(CNV_t)))
names(Ttest) <- names(table(CNV_t))
for (i in (1:length(table(CNV_t)))){
  test <- t.test(Input_CNV[which(CNV_t==names(table(CNV_t))[i])],Input_CNV[which(CNV_t==0)],var.equal=F,paired=F)
  Ttest[i] <- test$p.value
}
Ttest[-3]
p.adjust(Ttest[-3])

#### correlation with mutation data ####
which_subtype <- all_subtypes[1]
load(paste0(TargetDirectory,"Beta_",which_subtype,".Rdata"))
Beta <- Beta[!(rowSums(Beta) == 0),]

Input <- Beta
Samples <- intersect(colnames(Input),colnames(mutations))
Genes <- intersect(rownames(Input),rownames(mutations))
Input <- Input[Genes,Samples]
Mutations_tr <- mutations[Genes,Samples]
Mutations_thr <- Mutations_tr
Mutations_thr[(Mutations_thr=="NaN")] <- 0
Mutations_thr[!(Mutations_thr=="0")] <- 1
Mutations_thr <- mapply(Mutations_thr, FUN=as.numeric)
Mutations_thr <- matrix(data=Mutations_thr, ncol=ncol(Mutations_tr), nrow=nrow(Mutations_tr))
colnames(Mutations_thr) <- colnames(Mutations_tr)
rownames(Mutations_thr) <- rownames(Mutations_tr)

SignMut <- which(rowSums(Mutations_thr=="1")>0)
if (length(SignMut)>0){
  Mutations_thr<- Mutations_thr[SignMut,]
  Input <- Input[SignMut,]
}

No <- c()
PVal_Genes <- c()
for (i in (1:nrow(Mutations_thr))){
  Gene <- rownames(Mutations_thr)[i]
  Mutated <- which(Mutations_thr[i,]==1)
  if (length(Mutated)==0){
    No <- c(No,Gene)
  } else {
    C <- wilcox.test(Input[Gene,]~Mutations_thr[Gene,])
    PVal_Genes <- c(PVal_Genes,C$p.value)
  }
}
names(PVal_Genes) <- rownames(Mutations_thr)[which(is.na(match(rownames(Mutations_thr),No))==TRUE)]
PVal_Genes <- PVal_Genes[order(PVal_Genes,decreasing=FALSE)]
PVal_Genes_Adj <- p.adjust(PVal_Genes,method = "BH")
length(which(PVal_Genes_Adj<0.1))
length(which(PVal_Genes<0.1))

# plot
TF_CNV <- names(PVal_Genes)[which(PVal_Genes<0.1)]

for (i in (1:length(TF_CNV))){
  boxplot(Input[TF_CNV[i],]~Mutations_thr[TF_CNV[i],],main=paste0("TF ",TF_CNV[i]))
}

dfr <- matrix(0,ncol=3,nrow=length(TF_CNV)*ncol(Beta))
colnames(dfr) <- c("Beta","Gene","Mutation")
for (i in (1:length(TF_CNV))){
  gene <- TF_CNV[i]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),1] <- Input[gene,]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),2] <- rep(gene,ncol(Input))
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),3] <- Mutations_thr[gene,]
}
dfr <- data.frame(dfr)
dfr$Beta <- as.numeric(as.character(dfr$Beta))
p<-ggplot(dfr, aes(x=Gene,y=Beta, fill=Mutation)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("-2","-1","0","1","2"),values=c("red", "green","blue","yellow","black")) +ylim(0,0.15)

#### matrice de confusion beta != 0 vs. mutation
which_subtype <- all_subtypes[1]
load(paste0(TargetDirectory,"Beta_",which_subtype,".Rdata"))
Beta <- Beta[!(rowSums(Beta) == 0),]
Input <- Beta
Samples <- intersect(colnames(Input),colnames(mutations))
Genes <- intersect(rownames(Input),rownames(mutations))
Input <- Input[Genes,Samples]
Mutations_tr <- mutations[Genes,Samples]
Mutations_thr <- Mutations_tr
Mutations_thr[(Mutations_thr=="NaN")] <- 0
Mutations_thr[!(Mutations_thr=="0")] <- 1
Mutations_thr <- mapply(Mutations_thr, FUN=as.numeric)
Mutations_thr <- matrix(data=Mutations_thr, ncol=ncol(Mutations_tr), nrow=nrow(Mutations_tr))
colnames(Mutations_thr) <- colnames(Mutations_tr)
rownames(Mutations_thr) <- rownames(Mutations_tr)

beta <- Input
beta[which(beta>0)]<- 1

confusion <- matrix(0,ncol=2,nrow=2)
colnames(confusion) <- c("Beta_nnnull","Beta_null")
rownames(confusion) <- c("NnMut","Mut")  
confusion[2,2] <- length(which(beta-Mutations_thr==-1))
confusion[1,1] <- length(which(beta-Mutations_thr==1))
confusion[2,1] <- length(which(Mutations_thr==1))-confusion[2,2]
confusion[1,2] <- length(which(Mutations_thr==0))-confusion[1,1]
accuracy <- (confusion[1,1] + confusion[2,2])/(sum(confusion))
precision <- confusion[1,2] / (sum(confusion[,1]))  
rappel <- confusion[1,2]/sum(confusion[2,])
Fscore <- 2 * (precision*rappel)/(precision +rappel)

# tests aléatoires
alea <- function(v){
  s <- sample(length(v))
  v[s]
}
pvalue <- function(beta,Mutations_thr){
  alea1 <- alea(beta)
  alea2 <- alea(Mutations_thr)
  confusion <- matrix(0,ncol=2,nrow=2)
  colnames(confusion) <- c("Beta_nnnull","Beta_null")
  rownames(confusion) <- c("NnMut","Mut")  
  confusion[2,2] <- length(which(alea1-alea2==-1))
  confusion[1,1] <- length(which(alea1-alea2==1))
  confusion[2,1] <- length(which(alea2==1))-confusion[2,2]
  confusion[1,2] <- length(which(alea2==0))-confusion[1,1]
  accuracy <- (confusion[1,1] + confusion[2,2])/(sum(confusion))
  precision <- confusion[1,2] / (sum(confusion[,1]))  
  rappel <- confusion[1,2]/sum(confusion[2,])
  Fscore <- 2 * (precision*rappel)/(precision +rappel)
  return(c(accuracy,Fscore))
}
pvalue(beta,Mutations_thr)

results <- matrix(0,ncol=2,nrow=1e6)
colnames(results) <- c("accuracy","Fscore")
for (i in (1:1e6)){
  results[i,] <- pvalue(beta,Mutations_thr)
}
length(which(results[,1]>accuracy))/1e6
length(which(results[,2]<Fscore))/1e6

#### correlation with subtypes
which_subtype <- all_subtypes[5]
load(paste0(TargetDirectory,"Beta_",which_subtype,".Rdata"))
Beta <- Beta[!(rowSums(Beta) == 0),]

library(ComplexHeatmap) 
library(circlize)
# plotting
heatmap <- Heatmap(Beta,cluster_rows = T,cluster_columns = T,show_row_dend = F,
                   show_column_names=FALSE,
                   col=colorRamp2(c(min(Beta),median(Beta[Beta>0]),max(Beta)), c("blue","aquamarine","white")),
                   heatmap_legend_param = list(color_bar = "continuous"))
draw(heatmap)

# plot MA data by subtypes
p <- c()
for (i in (1:length(cibles))){
  if (length(cibles[[i]])>1){
    ModuleData <- t(MA)[cibles[[i]],]
    RegulatorData <- t(MA)[names(cibles)[i],]
    GeneClustering=hclust(dist(ModuleData), method = "complete", members = NULL)
  } else {
    ModuleData <- matrix(0,nrow=1,ncol=ncol(t(MA)))
    colnames(ModuleData) <- rownames(MA)
    rownames(ModuleData) <- cibles[[i]]
    ModuleData[1,] <- t(MA)[cibles[[i]],] 
    RegulatorData <- t(MA)[names(cibles)[i],]
    GeneClustering$order <- cibles[[i]]
  }
  
  SampleClustering=hclust(dist(RegulatorData), method = "complete", members = NULL)    

  ClustRegulatorData <- RegulatorData[SampleClustering$order]
  ClustModuleData <- ModuleData[GeneClustering$order,SampleClustering$order]
  ClustCombinedData <- rbind(ClustModuleData,ClustRegulatorData)
  rownames(ClustCombinedData) <- c(GeneClustering$order,names(cibles)[i])
  
  my_subtype <- rep("other",ncol(ClustCombinedData))  
  names(my_subtype) <- colnames(ClustCombinedData)
  my_subtype[which(colnames(ClustCombinedData) %in% colnames(Beta))] <- "subtype"

  data <- data.frame(ClustRegulatorData,my_subtype)
  boxplot(data$ClustRegulatorData~data$my_subtype,main=paste0("Subtype ",which_subtype, " TF ",names(cibles)[i]))
  C <- wilcox.test(ClustRegulatorData~my_subtype)
  p <- c(p,C$p.value)
  print(paste0("Pour le TF ",names(cibles)[i],", la p-value est de ",C$p.value))
}
names(p) <- names(cibles)
padjust <- p.adjust(p,method = "BH")
padjust

# plot the heatmap (if wanted)
heatmap_MA <- function(which_subtype,TF){
  load(paste0(TargetDirectory,"Beta_",which_subtype,".Rdata"))
  Beta <- Beta[!(rowSums(Beta) == 0),]
  
  GRN_LicornOutput <- read.table(paste0(TargetDirectory,"Network_",which_subtype,".txt"), sep=";", strip.white=TRUE, header=TRUE)
  GRN <- rbind(as.matrix(GRN_LicornOutput[, c(1,2)]),
               as.matrix(GRN_LicornOutput[, c(1,3)]))
  GRN <- GRN[GRN[, 2] != "", ]
  rep <- sapply(strsplit(GRN[, 2], " "), length)
  GRN <- cbind(rep(GRN[, 1], rep), unlist(strsplit(GRN[, 2], " ")))
  Adj_matrix <- as_adjacency_matrix(graph_from_edgelist(GRN, directed = TRUE),sparse=FALSE)
  Adj_matrix <- Adj_matrix[!(rowSums(Adj_matrix) == 0), ]
  Adj_matrix <- Adj_matrix[,!(colSums(Adj_matrix) == 0)]
  Adj_matrix <- Adj_matrix[order(rownames(Adj_matrix)), ]
  
  cibles <- names(which(Adj_matrix[,TF]==1))
  cibles <- str_replace_all(cibles,"_TF","")

  if (length(cibles)>1){
    ModuleData <- t(MA)[cibles,]
    RegulatorData <- t(MA)[TF,]
    GeneClustering=hclust(dist(ModuleData), method = "complete", members = NULL)
  } else {
    ModuleData <- matrix(0,nrow=1,ncol=ncol(t(MA)))
    colnames(ModuleData) <- rownames(MA)
    rownames(ModuleData) <- cibles
    ModuleData[1,] <- t(MA)[cibles,] 
    RegulatorData <- t(MA)[TF,]
    GeneClustering$order <- cibles
  }
  
  SampleClustering=hclust(dist(RegulatorData), method = "complete", members = NULL)    
  
  ClustRegulatorData <- RegulatorData[SampleClustering$order]
  ClustModuleData <- ModuleData[GeneClustering$order,SampleClustering$order]
  ClustCombinedData <- rbind(ClustModuleData,ClustRegulatorData)
  rownames(ClustCombinedData) <- c(GeneClustering$order,TF)
  
  my_subtype <- rep("other",ncol(ClustCombinedData))  
  names(my_subtype) <- colnames(ClustCombinedData)
  my_subtype[which(colnames(ClustCombinedData) %in% colnames(Beta))] <- "subtype"
  
  haRow = HeatmapAnnotation(test = my_subtype,name=paste0(which_subtype),
                            col = list(test=c("other"="white","subtype"="black")))
  
  heatmap <- Heatmap(ClustCombinedData, name = "Gene expression", column_title = paste('Module',TF),
                     cluster_rows=FALSE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,row_names_gp=gpar(col=c(rep("white",nrow(ModuleData)),"black"),fontsize=10),
                     column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Cibles",nrow(ModuleData)),"Regulators"),gap = unit(5, "mm"),
                     #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                     col=colorRamp2(c(min(ClustCombinedData), median(ClustCombinedData), max(ClustCombinedData)), c("green", "black", "red")),heatmap_legend_param = list(color_bar = "continuous"),bottom_annotation = haRow)
  draw(heatmap)
}

## Analyse différentielle
TFS <- read.table(paste0(DataDirectory,"/AllHumanTranscriptionFactor.txt"),sep="\t")
MA_TF <- MA[,TFS$V1]
MAcancerNorm = sweep(MA_TF, 2, colMeans(MA_TF), "-")
#Select the genes that are at least 2 fold above mean
OverExpressedGene=list()
#For each patient, select the genes that are >2, i.e. 2 fold higher than mean expression for that gene
for(j in 1:nrow(MAcancerNorm)){
  OverExpressedGene[[j]]=MAcancerNorm[j,MAcancerNorm[j,]>=2]
}
names(OverExpressedGene)=rownames(MAcancerNorm)
#Select the genes that are at least 2 fold below mean
UnderExpressedGene=list()
#For each patient, select the genes that are <2, i.e. 2 fold lower than mean expression for that gene
for(j in 1:nrow(MAcancerNorm)){
  UnderExpressedGene[[j]]=MAcancerNorm[j,MAcancerNorm[j,]<=-2]
}
names(UnderExpressedGene)=rownames(MAcancerNorm)

FoldChanges <- MAcancerNorm
FoldChanges[abs(FoldChanges)<2] <- 0
FoldChanges <- as.matrix(FoldChanges)
FoldChanges_reduced <- FoldChanges[,-which(colSums(FoldChanges)==0)]

ResultsTFs <- matrix(0,ncol=length(unique(Subtypes)),nrow=ncol(FoldChanges_reduced))
rownames(ResultsTFs) <- colnames(FoldChanges_reduced)
colnames(ResultsTFs) <- unique(Subtypes)
for (i in (1:length(unique(Subtypes)))){
  I = names(Subtypes)[which(Subtypes==unique(Subtypes)[i])]
  subtype = colSums(abs(FoldChanges_reduced[I,])>0)/length(I)
  ResultsTFs[,i] <- subtype
}
total = rowSums(ResultsTFs)
ResultsTFs <- ResultsTFs[order(total,decreasing=TRUE),]
write.csv2(ResultsTFs,file=paste0(TargetDirectory,"DeregulatedTFs_FoldChanges.csv"))

##### Differential expression analysis using limma
library(limma)
TFs <- read.table(paste0(DataDirectory,"/AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1
TFs <- as.character(TFs)
TFs <- intersect(TFs,colnames(MA))
Results <- matrix(0,ncol=2*length(all_subtypes),nrow=length(TFs))
rownames(Results) <- TFs
colnames(Results) <- rep(sort(unique(Subtypes)),each=2)
for (i in (1:length(unique(Subtypes)))){
  # define the subtypes
  sub <- sort(unique(Subtypes))[i]
  sub2 <- sort(unique(Subtypes))[-which(sort(unique(Subtypes))==sub)]
  
  Sub_tmp <- Subtypes
  Sub_tmp[Sub_tmp==sub] <- "sub"
  Sub_tmp[Sub_tmp %in% sub2] <- "other"
  Sub_tmp <- as.factor(Sub_tmp)
  
  # define the design matrix
  design <- model.matrix(~ 0 + Sub_tmp)
  colnames(design) <- levels(Sub_tmp)
  
  # fit the linear model
  fit <- lmFit(t(MA[,TFs]),design)
  
  # define the contrasts
  cont.matrix <- makeContrasts(test_sub=other - sub,levels=design)
  fit.cont <- contrasts.fit(fit, cont.matrix)
  
  # ebayes
  fit.cont <- eBayes(fit.cont)
  
  # summary of tests
  summa.fit <- decideTests(fit.cont)
  summary(summa.fit)
  
  # output genes
  gene_table <- topTable(fit.cont,coef="test_sub",sort.by="p",n="Inf")
  Results[rownames(gene_table),((i*2)-1):(2*i)] <- as.matrix(gene_table[c("logFC","adj.P.Val")],ncol=ncol(gene_table),nrow=nrow(gene_table))
}
write.csv2(Results,file=paste0(TargetDirectory,"DeregulatedTFs_FoldChanges.csv"))

# Venn diagram
library(VennDiagram)
Beta_sparse <- read.csv2(paste0(TargetDirectory,"/Beta_sparse.csv"),dec=".")
rownames(Beta_sparse) <- Beta_sparse$X
Beta_sparse <- Beta_sparse[,-1]

Diff_genes <- read.csv2(paste0(TargetDirectory,"/DeregulatedTFs_FoldChanges.csv"))
rownames(Diff_genes) <- Diff_genes$X
Diff_genes <- Diff_genes[,-1]

par(mfrow=c(1,1))
threshold <- 0.5
for (i in (1:ncol(Beta_sparse))){
  I1 <- rownames(Beta_sparse)[(which(Beta_sparse[,i]>threshold))]
  I2 <- rownames(Diff_genes)[which(Diff_genes[,2*i]<0.01)]
  g <- draw.pairwise.venn(length(I1),length(I2), length(intersect(I1,I2)), c("Diff exprimés","Dérégulés "),col=c("green","red"),scaled=FALSE,euler.d=TRUE,ind=FALSE)
  require(gridExtra)
  grid.arrange(gTree(children=g), top=paste0("Venn Diagram for ",colnames(Beta_sparse)[i], " subtype"))
}

# quelques analyses sur les cibles
Target_analysis <- function(which_subtype,threshold,TF,Beta=0,Adj_matrix,Score,DataBase){
  # TF : précise le TF étudié (obligatoire)
  # threshold : à partir de quel moment on considère le TF comme dérégulé
  # which_subtype : sur quel sous-type travaille t-on
  # Beta, Adj_matrix, Score : on peut directement charger Beta, Adj_matrix et Score pour gagner du temps
  # DataBase : base de données à charger pour l'analyse d'enrichissement
  
  if (Beta==0){
    load(paste0(TargetDirectory,"Beta_",which_subtype,".Rdata"))
    Beta <- Beta[!(rowSums(Beta) == 0),]
    
    GRN_LicornOutput <- read.table(paste0(TargetDirectory,"Network_",which_subtype,".txt"), sep=";", strip.white=TRUE, header=TRUE)
    GRN <- rbind(as.matrix(GRN_LicornOutput[, c(1,2)]),
                 as.matrix(GRN_LicornOutput[, c(1,3)]))
    GRN <- GRN[GRN[, 2] != "", ]
    rep <- sapply(strsplit(GRN[, 2], " "), length)
    GRN <- cbind(rep(GRN[, 1], rep), unlist(strsplit(GRN[, 2], " ")))
    Adj_matrix <- as_adjacency_matrix(graph_from_edgelist(GRN, directed = TRUE),sparse=FALSE)
    Adj_matrix <- Adj_matrix[!(rowSums(Adj_matrix) == 0), ]
    Adj_matrix <- Adj_matrix[,!(colSums(Adj_matrix) == 0)]
    Adj_matrix <- Adj_matrix[order(rownames(Adj_matrix)), ]
    
    Score <- read.table(paste0(TargetDirectory,"Scores_",which_subtype,".txt"),sep=',')
    colnames(Score) <- gsub(".", '-', colnames(Score), fixed = T)
    Score <- Score[!(rowSums(Score) == 0), ]
    Score <- Score[order(rownames(Score)), ]
    Score <- as.matrix(Score)
  }  
  
  if ((rowSums(Beta>0)[TF]/ncol(Beta))<threshold){
    print(paste0("On ne considère pas le TF comme dérégulé car il est non nul dans ",round(rowSums(Beta>0)[TF]/ncol(Beta),2)*100,"% des cas"))
  } else {
    print(paste0("Le TF est dérégulé car il est non nul dans ",round(rowSums(Beta>0)[TF]/ncol(Beta),2)*100,"% des cas"))
    
    # def des cibles
    cibles <- names(which(Adj_matrix[,TF]==1))
    cibles_TF <- str_replace_all(cibles,"_TF","")
    
    # cibles dérégulés
    Score_cibles <- Score[cibles,]
    score <- sort(rowMeans(Score_cibles),decreasing = T)
    
    # GSEA pour les cibles
    Results <- EnrichmentAnalysis(Data=cibles_TF,DataBase=DataBase)
    GSEA <- Results$ResultThreshold
    
    res <- list(cibles = cibles,scores=score,GSEA)
    names(res) <- c("cibles","score_cibles","GSEA")
    return(res)
  }
}

library(limma)
library(GSA)
DataBase <- CreateDataBase() # run the code directly
Target_analysis(which_subtype=all_subtypes[1],threshold=0.5,TF="PPARG",DataBase=DataBase)
