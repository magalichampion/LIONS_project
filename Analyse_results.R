### Analyse the different results

# libraries
library(stringr)
library(ggplot2)
library(igraph)

##### load the different data types #####
# paths
path <- "/Users/Magali/Desktop/recherche/LIONS/"
TargetDirectory <- paste0(path,"results/")

# network
GRN_LicornOutput <- read.table(paste0(TargetDirectory,"Network_TCGA2.txt"), sep=";", strip.white=TRUE, header=TRUE)
GRN <- rbind(as.matrix(GRN_LicornOutput[, c(1,2)]),
             as.matrix(GRN_LicornOutput[, c(1,3)]))
GRN <- GRN[GRN[, 2] != "", ]
rep <- sapply(strsplit(GRN[, 2], " "), length)
GRN <- cbind(rep(GRN[, 1], rep), unlist(strsplit(GRN[, 2], " ")))
Adj_matrix <- as_adjacency_matrix(graph_from_edgelist(GRN, directed = TRUE),sparse=FALSE)
Adj_matrix <- Adj_matrix[!(rowSums(Adj_matrix) == 0), ]
Adj_matrix <- Adj_matrix[,!(colSums(Adj_matrix) == 0)]
Adj_matrix <- Adj_matrix[order(rownames(Adj_matrix)), ]

# scores
Score <- read.table(paste0(TargetDirectory,"Scores_TCGA2_Corrected.txt"),sep=',')
colnames(Score) <- gsub(".", '-', colnames(Score), fixed = T)
#colnames(Score) <- str_replace_all(colnames(Score),"X","")
Score <- Score[!(rowSums(Score) == 0), ]
Score <- Score[order(rownames(Score)), ]
Score <- as.matrix(Score)

# Beta
load(paste0(TargetDirectory,"Beta_TCGA2_Corrected.Rdata"))
Beta <- Beta[!(rowSums(Beta) == 0),]

# subtypes data
#load(paste0(path,"Data/Subtypes/CIT_BLCA_Subgroup.RData"))
#Subgroups <- CIT_BLCA_Subgroup
#Samples <- sub('[A-Z]{3}$', "", names(Subgroups))
#names(Subgroups) <- Samples
Subgroups_full <- read.table(paste0(path,"Data/Subtypes/BLCA-TP.mergedcluster.csv"),sep=",",header=TRUE)
rownames(Subgroups_full) <- Subgroups_full$SampleName  
Subgroups <- Subgroups_full$mRNAseq_cHierarchical
names(Subgroups) <- rownames(Subgroups_full)

Subgroups <- as.character(Subgroups)
names(Subgroups) <- rownames(Subgroups_full)
Subgroups <- Subgroups[colnames(Beta)]
Subgroups[which(Subgroups=="1")] <- rep("p53",length(which(Subgroups==1)))
Subgroups[which(Subgroups=="2")] <- rep("luminal",length(which(Subgroups==2)))
Subgroups[which(Subgroups=="3")] <- rep("basal",length(which(Subgroups==3)))

# CNV
#load(paste0(path,"Data/copyNumber/CIT/CIT_BLCA_CNA.rda"))
#colnames(CIT_BLCA_CNA) <- str_replace_all(colnames(CIT_BLCA_CNA),"[A-Z]","")
#CNV <- CIT_BLCA_CNA
load(paste0(path,"Data/TCGA/ProcessedData_BLCA.Rdata"))
CNV <- ProcessedData$CNV_TCGA_thresholded
colnames(CNV) <- paste0(colnames(CNV),"-01")

# mutation
Mutations <- read.table(paste0(path,"Data/TCGA/Mutation_TCGA2_Corrected.csv"),header=TRUE,sep=",")
rownames(Mutations) <- Mutations[,1]
Mutations <- Mutations[,-c(1,2)]
colnames(Mutations) <- gsub(".", '-', colnames(Mutations), fixed = T)
Mutations <- Mutations[-1,]

# gene expression
#MA_cancer <- read.table(paste0(path,"Data/expression/TumorExpressionData_AllGenes.txt"))
#rownames(MA_cancer) <- str_replace_all(rownames(MA_cancer),"[A-Z]","") # anonymous data
#VarMax <- 0.75
#vars <- apply(MA_cancer, 2, var)
#MA_cancer <- MA_cancer[,order(vars, decreasing=TRUE)[1:floor(VarMax*ncol(MA_cancer))]]

# TFs
TFs <- read.table(paste0(path,"Data/expression/AllHumanTranscriptionFactor.txt"))
TFs <- TFs$V1 

#### Network ####
Adj_matrix_bin <- matrix(0,ncol=length(union(colnames(Adj_matrix),rownames(Adj_matrix))),nrow=length(union(colnames(Adj_matrix),rownames(Adj_matrix))))
colnames(Adj_matrix_bin) = rownames(Adj_matrix_bin) <- union(colnames(Adj_matrix),rownames(Adj_matrix))
Adj_matrix_bin[rownames(Adj_matrix),colnames(Adj_matrix)] <- Adj_matrix
Adj_matrix_bin <- t(Adj_matrix_bin) + Adj_matrix_bin
net1 <- graph.adjacency(Adj_matrix_bin)
lay1 <- layout.kamada.kawai(net1)
png("Network.png")
plot(net1,main="", layout=lay1)
dev.off()

##### Histograms #####
Hist <- hist(Score)
Hist$counts <- log(Hist$counts)
plot(Hist,main="",xlab="Deregulation scores",ylab="Frequency",yaxt="n")
axis(side=2,at=c(0,2.3,4.6,6.9,9.2,11.5,13.8),labels=c("1","10","1e2","1e3","1e4","1e5","1e6"))

Hist <- hist(Beta)
Hist$counts <- log(Hist$counts)
plot(Hist,main="",xlab="Deregulation scores",ylab="Frequency",yaxt="n")
axis(side=2,at=c(0,2.3,4.6,6.9,9.2,11.5,13.8),labels=c("1","10","1e2","1e3","1e4","1e5","1e6"))

##### CNV comparison #####
Comp <- "Beta" # you can also choose Score as inputs
if (Comp=="Beta"){
  Input <- Beta
} else {
  Input <- Score  
  rownames(Input) <- str_replace_all(rownames(Input),"_TF","")
}
Genes <- intersect(rownames(Input),rownames(CNV))
Samples <- intersect(colnames(Input),colnames(CNV))
CNV <- as.matrix(CNV[Genes,Samples])
Input <- Input[Genes,Samples]

# plot the distribution
I <- which(Input<0.05)
if (length(I)>0){
  Input_m <- Input[-I]
  CNV_m <- CNV[-I]
} else {
  Input_m <- Input
  CNV_m <- CNV
}
par(mfrow = c(1, 1))
I2 <- which(CNV_m==2)
I1 <- which(CNV_m==1)  
I1m <- which(CNV_m==-1)  
I2m <- which(CNV_m==-2)
I0 <- which(CNV_m==0)

F2 <- ecdf(Input_m[I2])
F1 <- ecdf(Input_m[I1])
F1m <- ecdf(Input_m[I1m])
F2m <- ecdf(Input_m[I2m])
F0 <- ecdf(Input_m[I0])

plot(F2,verticals=TRUE,do.points=FALSE,main="Cumulative empirical distribution",xlab="Score of deregulation",ylab="")
lines(F1,verticals=TRUE,do.points=FALSE,col="green")
lines(F1m,verticals=TRUE,do.points=FALSE,col="red")
lines(F2m,verticals=TRUE,do.points=FALSE,col="yellow")
lines(F0,verticals=TRUE,do.points=FALSE,col="blue")
legend(0.8,0.4,c("CNV=2","CNV=1","CNV=-1","CNV=-2","CNV=0"),lty=rep(1,5),col=c("black","green","red","yellow","blue"))

dfr <- matrix(0,ncol=4,nrow=(nrow(Input)*ncol(Input)+nrow(Input_corrected)*ncol(Input_corrected)))
colnames(dfr) <- c("Score","Gene","CNV","Correction")
for (i in (1:nrow(Input))){
  gene <- rownames(Input)[i]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),1] <- Input[gene,]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),2] <- rep(gene,ncol(Input))
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),3] <- CNV[gene,]
}
for (i in ((nrow(Input)+1):(nrow(Input)+nrow(Input_corrected)))){
  gene <- rownames(Input_corrected)[i-nrow(Input)]
  dfr[(ncol(Input_corrected)*(i-1)+1):(ncol(Input_corrected)*i),1] <- Input_corrected[gene,]
  dfr[(ncol(Input_corrected)*(i-1)+1):(ncol(Input_corrected)*i),2] <- rep(gene,ncol(Input_corrected))
  dfr[(ncol(Input_corrected)*(i-1)+1):(ncol(Input_corrected)*i),3] <- CNV[gene,]
}
dfr[,4] <- c(rep("Non corrected",nrow(Input)*ncol(Input)),rep("Corrected",nrow(Input_corrected)*ncol(Input_corrected)))

I <- which(is.na(dfr[,3])==TRUE)
if (length(I)>0){
  dfr <- dfr[-I,]
}
dfr <- data.frame(dfr)
dfr$Score <- as.numeric(as.character(dfr$Score))

p<-ggplot(dfr, aes(x=Correction,y=Score, fill=CNV)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("-2","-1","0","1","2"),values=c("green", "red","blue","yellow","black"))

# Student test for testing the differences
t.test(Input_m[I2],Input_m[I0],var.equal=F,paired=F)

# List of significantly associated genes
No <- c()
PVal_Genes <- c()
for (i in (1:nrow(Input))){
  Gene <- rownames(Input)[i]
  if (dim(table(CNV[Gene,]))==1){
    No <- c(No,Gene)
  } else {
    I <- which(is.na(CNV[Gene,])==TRUE)
    if (length(I)>0){
      if (length(I)==ncol(Input)){
        # no data
        No <- c(No,Gene)
      } else {
        C <- kruskal.test(Input[Gene,-I]~CNV[Gene,-I])
        PVal_Genes <- c(PVal_Genes,C$p.value)
      } 
    } else {
      C <- kruskal.test(Input[Gene,]~CNV[Gene,])
      PVal_Genes <- c(PVal_Genes,C$p.value)
    }
  }
}
names(PVal_Genes) <- rownames(Input)[which(is.na(match(rownames(Input),No))==TRUE)]
PVal_Genes <- PVal_Genes[order(PVal_Genes,decreasing=FALSE)]

Adj_pval <- p.adjust(PVal_Genes,method="BH")
TF_CNV <- names(PVal_Genes)[which(p.adjust(PVal_Genes,method="BH")<0.05)]

for (i in (1:length(TF_CNV))){
  boxplot(Input[TF_CNV[i],]~CNV[TF_CNV[i],],main=paste0("TF ",TF_CNV[i]))
}

dfr <- matrix(0,ncol=3,nrow=length(TF_CNV)*ncol(Input))
colnames(dfr) <- c("Beta","Gene","CNV")
for (i in (1:length(TF_CNV))){
  gene <- TF_CNV[i]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),1] <- Input[gene,]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),2] <- rep(gene,ncol(Input))
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),3] <- CNV[gene,]
}
I <- which(is.na(dfr[,3])==TRUE)
dfr <- dfr[-I,]
dfr <- data.frame(dfr)
dfr$Beta <- as.numeric(as.character(dfr$Beta))
p<-ggplot(dfr, aes(x=Gene,y=Beta, fill=CNV)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("-2","-1","0","1","2"),values=c("red", "green","blue","yellow","black"))

# distribution of these genes
load(paste0(TargetDirectory,"TF_CNV_To_correct.Rdata"))
Input_TF <- Input[TF_CNV,]
CNV_TF <- CNV[TF_CNV,]
I <- which(Input_TF==0)
par(mfrow = c(1, 1))
I2 <- which(CNV_TF==2)
I1 <- which(CNV_TF==1)  
I1m <- which(CNV_TF==-1)  
I2m <- which(CNV_TF==-2)
I0 <- which(CNV_TF==0)

if (length(I)>0){
  Input_m <- Input_TF[-I]
} else {
  Input_m <- Input_TF
}
F2 <- ecdf(Input_m[I2])
F1 <- ecdf(Input_m[I1])
F1m <- ecdf(Input_m[I1m])
F2m <- ecdf(Input_m[I2m])
F0 <- ecdf(Input_m[I0])

plot(F2,verticals=TRUE,do.points=FALSE,main="Cumulative empirical distribution",xlab="Score of deregulation",ylab="")
lines(F1,verticals=TRUE,do.points=FALSE,col="green")
lines(F1m,verticals=TRUE,do.points=FALSE,col="red")
lines(F2m,verticals=TRUE,do.points=FALSE,col="yellow")
lines(F0,verticals=TRUE,do.points=FALSE,col="blue")
legend(0.8,0.4,c("CNV=2","CNV=1","CNV=-1","CNV=-2","CNV=0"),lty=rep(1,5),col=c("black","green","red","yellow","blue"))

t.test(Input_m[I2],Input_m[I0],var.equal=F,paired=F)

#### Beta - table ####
Beta_small <- Beta[,names(Subgroups)]
ResultsTFs <- matrix(0,ncol=length(unique(Subgroups)),nrow=nrow(Beta_small))
rownames(ResultsTFs) <- rownames(Beta_small)
colnames(ResultsTFs) <- unique(Subgroups)
for (i in (1:length(unique(Subgroups)))){
  I = names(Subgroups)[which(Subgroups==unique(Subgroups)[i])]
  subtype = rowMeans(Beta_small[,I])
  ResultsTFs[,i] <- subtype
}
total = rowMeans(Beta_small)
ResultsTFs <- cbind(ResultsTFs,total)
ResultsTFs <- ResultsTFs[order(total,decreasing=TRUE),]
TargetDirectory <- paste0(path,"results/")
write.csv2(ResultsTFs,file=paste0(TargetDirectory,"DeregulatedTFs_TCGA2_Corrected.csv"))

#### expression diff ####
MA_cancer_norm <- MA_cancer[names(Subgroups),]
MA_cancer_norm = MA_cancer_norm - matrix(1,nrow(MA_cancer_norm),1) %*% colMeans(MA_cancer_norm)
MA_cancer_norm = MA_cancer_norm / sqrt( (1/nrow(MA_cancer_norm)) * matrix(1,nrow(MA_cancer_norm),1) %*% colSums(MA_cancer_norm^2))

library(limma)
Design <- matrix(0,ncol=length(unique(Subgroups)),nrow=length(Subgroups))
colnames(Design) <- unique(Subgroups)
Names <- c()
for (i in (1:length(unique(Subgroups)))){
  I <- which(Subgroups==unique(Subgroups)[i])
  Design[I,i] <- rep(1,length(I))
  Names <- c(Names,names(I))
}
rownames(Design) <- Names

# for each of the subtypes 
DEG <- list()
MA_cancer_norm <- MA_cancer_norm[rownames(Design),]
for (i in (1:length(unique(Subgroups)))){
  Design2 <- cbind(Design[,i],rowSums(Design[,-i]))
  #colnames(Design) <- c("Pap","Lum","Basal","tcga")
  fit <- lmFit(t(MA_cancer_norm),Design2)
  #constrat.matrix <- makeContrasts("Pap-Lum","Pap-Basal","Pap-tcga","Lum-Basal","Lum-tcga","Basal-tcga",levels=Design)
  #fit2 <- contrasts.fit(fit,constrat.matrix)
  fit <- eBayes(fit)
  Expressed <- topTable(fit,coef=1,sort.by = "logFC" ,number=ncol(MA_cancer_norm),p.value = 0.05,lfc = 0.1375)
  deg <- Expressed[,1]
  names(deg) <- rownames(Expressed)
  tfs <- intersect(TFs,names(deg))
  deg <- deg[tfs]
  deg <- deg[order(abs(deg),decreasing=TRUE)]
  DEG <- c(DEG,list(deg))
}
names(DEG) <- unique(Subgroups)

# the same for dereg scores (exactly the same as our computed averaged score - except p-values)
DSG <- list()
Beta <- Beta[,rownames(Design)]
for (i in (1:length(unique(Subgroups)))){
  Design2 <- cbind(Design[,i],rowSums(Design[,-i]))
  #colnames(Design) <- c("Pap","Lum","Basal","tcga")
  fit <- lmFit(Beta,Design2)
  #constrat.matrix <- makeContrasts("Pap-Lum","Pap-Basal","Pap-tcga","Lum-Basal","Lum-tcga","Basal-tcga",levels=Design)
  #fit2 <- contrasts.fit(fit,constrat.matrix)
  fit <- eBayes(fit)
  Expressed <- topTable(fit,coef=1,sort.by = "logFC" ,number=nrow(Beta),p.value = 0.05,lfc = 0.263)
  deg <- Expressed[,1]
  names(deg) <- rownames(Expressed)
  tfs <- intersect(TFs,names(deg))
  deg <- deg[tfs]
  deg <- deg[order(abs(deg),decreasing=TRUE)]
  DSG <- c(DSG,list(deg))
}
names(DSG) <- unique(Subgroups)

# the same for CNV
CNV <- CIT_BLCA_CNA[,names(Subgroups)]
DCNVG <- list()
for (i in (1:length(unique(Subgroups)))){
  Design2 <- cbind(Design[,i],rowSums(Design[,-i]))
  #colnames(Design) <- c("Pap","Lum","Basal","tcga")
  fit <- lmFit(CNV,Design2)
  #constrat.matrix <- makeContrasts("Pap-Lum","Pap-Basal","Pap-tcga","Lum-Basal","Lum-tcga","Basal-tcga",levels=Design)
  #fit2 <- contrasts.fit(fit,constrat.matrix)
  fit <- eBayes(fit)
  Expressed <- topTable(fit,coef=1,sort.by = "logFC" ,number=nrow(CNV),p.value = 0.05,lfc = 0.1375)
  deg <- Expressed[,1]
  names(deg) <- rownames(Expressed)
  tfs <- intersect(TFs,names(deg))
  deg <- deg[tfs]
  deg <- deg[order(abs(deg),decreasing=TRUE)]
  DCNVG <- c(DCNVG,list(deg))
}
names(DCNVG) <- unique(Subgroups)

# Venn Diagram
library(VennDiagram)
for (i in (1:length(unique(Subgroups)))){
  VennDiagram::draw.triple.venn(area1 = length(DEG[[i]]), area2 = length(DSG[[i]]), area3 = length(DCNVG[[i]]),
                   n12 = length(intersect(names(DEG[[i]]),names(DSG[[i]]))),n23=length(intersect(names(DSG[[i]]),names(DCNVG[[i]]))),n13=length(intersect(names(DEG[[i]]),names(DCNVG[[i]]))),
                   n123 = length(intersect(intersect(names(DEG[[i]]),names(DSG[[i]])),names(DCNVG[[i]]))),
                   category = c("Differentially expressed genes","Differentially perturbed genes","Differentially CNV altered genes"),
                   lty="blank",
                   fill=c("skyblue","pink1","mediumorchid")
)
}

#### correlation with mutation data
Samples <- intersect(colnames(Beta),colnames(Mutations))
Input <- Beta
Input <- Input[,Samples]
Input <- Input[!(rowSums(Input) == 0),]
Genes <- intersect(rownames(Input),rownames(Mutations))
Input <- Input[Genes,]
Mutations <- Mutations[Genes,Samples]
Mutation_threshold <- as.matrix(Mutations)
Mutation_thr <- Mutation_threshold

Mutation_thr[which(is.na(Mutation_threshold)==TRUE)] <- rep(0,length(which(is.na(Mutation_threshold)==TRUE)))
Mutation_thr[which(Mutation_threshold=="NaN")] <- rep(0,length(which(Mutation_threshold=="NaN")))
Mutation_thr[-which(Mutation_thr=="0")] <- rep(1,length(Mutation_thr)-length(which(Mutation_thr=="0")))

No <- c()
PVal_Genes <- c()
for (i in (1:nrow(Mutations))){
  Gene <- rownames(Mutations)[i]
  if (dim(table(Mutation_thr[Gene,]))==1){
    No <- c(No,Gene)
  } else {
      C <- wilcox.test(Input[Gene,]~Mutation_thr[Gene,])
      PVal_Genes <- c(PVal_Genes,C$p.value)
  }
}
names(PVal_Genes) <- rownames(Mutation_thr)[which(is.na(match(rownames(Mutation_thr),No))==TRUE)]
PVal_Genes <- PVal_Genes[order(PVal_Genes,decreasing=FALSE)]

# restrict to genes with enough mutation data
SignMut <- which(rowSums(Mutation_thr=="1")>4)
PVal_Genes2 <- c()
for (i in (1:length(SignMut))){
  Gene <- names(SignMut)[i]
  I <- which(Beta[Gene,]<0.05)
  Beta_small <- Beta
  Beta_small[Gene,I] <- rep(0,length(I))
  C <- wilcox.test(Beta_small[Gene,]~Mutation_thr[Gene,])
  PVal_Genes2 <- c(PVal_Genes2,C$p.value)
}
names(PVal_Genes2) <- names(SignMut)
PVal_Genes2 <- PVal_Genes2[!is.nan(PVal_Genes2)]
PVal_Genes2 <- PVal_Genes2[order(PVal_Genes2,decreasing=FALSE)]

Adj_pval <- p.adjust(PVal_Genes,method="BH")
TF_CNV <- names(PVal_Genes)[which(p.adjust(PVal_Genes,method="BH")<0.05)]

for (i in (1:length(TF_CNV))){
  boxplot(Input[TF_CNV[i],]~Mutation_thr[TF_CNV[i],],main=paste0("TF ",TF_CNV[i]))
}

dfr <- matrix(0,ncol=3,nrow=length(TF_CNV)*ncol(Input))
colnames(dfr) <- c("Beta","Gene","CNV")
for (i in (1:length(TF_CNV))){
  gene <- TF_CNV[i]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),1] <- Input[gene,]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),2] <- rep(gene,ncol(Input))
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),3] <- CNV[gene,]
}
I <- which(is.na(dfr[,3])==TRUE)
dfr <- dfr[-I,]
dfr <- data.frame(dfr)
dfr$Beta <- as.numeric(as.character(dfr$Beta))
p<-ggplot(dfr, aes(x=Gene,y=Beta, fill=CNV)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("-2","-1","0","1","2"),values=c("red", "green","blue","yellow","black"))

MutatedGenes <- names(PVal_Genes[which(PVal_Genes<0.05)])
dfr <- matrix(0,ncol=3,nrow=length(MutatedGenes)*ncol(Input))
colnames(dfr) <- c("Beta","Gene","Mutation")
for (i in (1:length(MutatedGenes))){
  gene <- MutatedGenes[i]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),1] <- Input[gene,]
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),2] <- rep(gene,ncol(Input))
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),3] <- Mutation_thr[gene,]
}
dfr <- data.frame(dfr)
dfr$Beta <- as.numeric(as.character(dfr$Beta))
p<-ggplot(dfr, aes(x=Gene,y=Beta, fill=Mutation_thr)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("0","1"),values=c("red", "green"))

#### correlation with subtypes
Beta_sub <- Beta[,names(Subgroups)]
PVal <- matrix(0,nrow=length(unique(Subgroups)),ncol=nrow(Beta_sub))
for (i in (1:length(unique(Subgroups)))){
  Sub <- Subgroups
  sub <- which(Subgroups==unique(Subgroups)[i])
  Sub[sub] <- rep(1,length(sub))
  Sub[-sub] <- rep(0,length(Sub)-length(sub))
  for (j in (1:nrow(Beta_sub))){
    gene <- rownames(Beta_sub)[j]
    c <- wilcox.test(Beta_sub[gene,]~Sub)
    PVal[i,j] <- c$p.value
  }
  PVal[i,] <- p.adjust(PVal[i,],method="BH")
}
colnames(PVal) <- rownames(Beta_sub)
rownames(PVal) <- unique(Subgroups)

for (i in (1:length(unique(Subgroups)))){
  Sub <- Subgroups
  sub <- which(Subgroups==unique(Subgroups)[i])
  Sub[sub] <- rep(1,length(sub))
  Sub[-sub] <- rep(0,length(Sub)-length(sub))
  
  ModuleData <- Beta_sub[which(PVal[i,]<0.05),]
  SampleClustering=hclust(dist(t(ModuleData)), method = "complete", members = NULL)    
  GeneClustering=hclust(dist(ModuleData), method = "complete", members = NULL)
  ClustModuleData <- ModuleData[GeneClustering$order,SampleClustering$order]

  # create annotations 
  Sub_data <- matrix(0,nrow=1,ncol=ncol(ClustModuleData))
  colnames(Sub_data) <- colnames(ClustModuleData)
  rownames(Sub_data) <- "Subtype"
  Sub_data[1,] <- Sub[colnames(ClustModuleData)]
  Genes <- data.frame(t(Sub_data),No=rep(" ",ncol(ClustModuleData)))
  ColAnnotation <- c("1" = "black", "0" = "white"," "="white")
  ColAnnotation <- rep(list(ColAnnotation),ncol(Genes))
  names(ColAnnotation) <- colnames(Genes)
  haRow = HeatmapAnnotation(df=Genes ,name="test",
                            col = ColAnnotation,which="column",show_legend=FALSE)
  
  # plotting
  heatmap <- Heatmap(log(ClustModuleData), name = "Gene expression", column_title = 'Module', cluster_rows=FALSE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,row_names_gp=gpar(col=c(rep("white",nrow(ClustModuleData))),fontsize=10),
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Module Genes",nrow(ClustModuleData))),gap = unit(5, "mm"),
                   #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                   col=colorRamp2(c(min(log(ClustModuleData)),-10,0), c("blue","aquamarine","white")),heatmap_legend_param = list(color_bar = "continuous"),top_annotation=haRow)
  draw(heatmap)
}
  
#### heatmap for Beta 
library(ComplexHeatmap)
library(circlize)
Beta_clus <- c()
Sub_clus <- c()
for (j in (1:length(unique(Subgroups)))){
  I <- which(as.character(Subgroups)==unique(as.character(Subgroups))[j])
  EnrichSmall <- Beta_small[,I]
  SampleClustering <- hclust(dist(t(EnrichSmall)), method = "complete", members = NULL)    
  EnrichSmall <- EnrichSmall[,SampleClustering$order]
  Beta_clus <- cbind(Beta_clus,EnrichSmall)
  Sub_small <- Subgroups[I]
  Sub_clus <- c(Sub_clus,as.character(Sub_small[SampleClustering$order]))
}
#GeneClustering <- hclust(dist(Beta_clus), method = "complete", members = NULL)
#Beta_clus <- Beta_clus[GeneClustering$order,]
Subtypes <- data.frame(Sub_clus)

ha = HeatmapAnnotation(df = Subtypes)

heatmap <- Heatmap(Beta_clus, name = "Score", cluster_rows=TRUE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,show_row_names=TRUE,row_names_gp=gpar(col=c(rep("black",nrow(Beta_clus))),fontsize=8),
                     column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Transcription factors",nrow(Beta_clus))),gap = unit(5, "mm"),
                     #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                     col=colorRamp2(c(0,0.1,max(Beta_clus)), c("white","aquamarine","blue")),heatmap_legend_param = list(color_bar = "continuous"),top_annotation=ha)
draw(heatmap)  

# top varying beta genes
library(matrixStats)
var <- rowVars(Beta)
names(var) <- rownames(Beta)
var <- var[order(var,decreasing = TRUE)]
VarGenes <- names(var)[1:(10/100*length(var))]
Beta_small <- Beta[VarGenes,]
