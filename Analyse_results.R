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
GRN_LicornOutput <- read.table(paste0(TargetDirectory,"Network_Corrected_First.txt"), sep=";", strip.white=TRUE, header=TRUE)
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
Score <- read.table(paste0(TargetDirectory,"Scores_Corrected_First.txt"),sep=',')
colnames(Score) <- str_replace_all(colnames(Score),"X","")
Score <- Score[!(rowSums(Score) == 0), ]
Score <- Score[order(rownames(Score)), ]
Score <- as.matrix(Score)

# Beta
load(paste0(TargetDirectory,"Beta_Corrected_First.Rdata"))
Beta <- Beta[!(rowSums(Beta) == 0),]

# subtypes data
load(paste0(path,"Data/Subtypes/CIT_BLCA_Subgroup.RData"))
Subgroups <- CIT_BLCA_Subgroup
Samples <- sub('[A-Z]{3}$', "", names(Subgroups))
Subgroups <- as.character(Subgroups)
names(Subgroups) <- Samples

# CNV
load(paste0(path,"Data/copyNumber/CIT/CIT_BLCA_CNA.rda"))
colnames(CIT_BLCA_CNA) <- str_replace_all(colnames(CIT_BLCA_CNA),"[A-Z]","")

# gene expression
MA_cancer <- read.table(paste0(path,"Data/expression/TumorExpressionData_AllGenes.txt"))
rownames(MA_cancer) <- str_replace_all(rownames(MA_cancer),"[A-Z]","") # anonymous data
VarMax <- 0.75
vars <- apply(MA_cancer, 2, var)
MA_cancer <- MA_cancer[,order(vars, decreasing=TRUE)[1:floor(VarMax*ncol(MA_cancer))]]

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
Genes <- intersect(rownames(Input),rownames(CIT_BLCA_CNA))
Samples <- intersect(colnames(Input),colnames(CIT_BLCA_CNA))
CNV <- as.matrix(CIT_BLCA_CNA[Genes,Samples])
Input <- Input[Genes,Samples]

# plot the distribution
I <- which(Input==0)
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
  dfr[(ncol(Input)*(i-1)+1):(ncol(Input)*i),3] <- CNV_non_corrected[gene,]
}
for (i in ((nrow(Input)+1):(nrow(Input)+nrow(Input_corrected)))){
  gene <- rownames(Input_corrected)[i-nrow(Input)]
  dfr[(ncol(Input_corrected)*(i-1)+1):(ncol(Input_corrected)*i),1] <- Input_corrected[gene,]
  dfr[(ncol(Input_corrected)*(i-1)+1):(ncol(Input_corrected)*i),2] <- rep(gene,ncol(Input_corrected))
  dfr[(ncol(Input_corrected)*(i-1)+1):(ncol(Input_corrected)*i),3] <- CNV[gene,]
}
dfr[,4] <- c(rep("Non corrected",nrow(Input)*ncol(Input)),rep("Corrected",nrow(Input_corrected)*ncol(Input_corrected)))

I <- which(is.na(dfr[,3])==TRUE)
dfr <- dfr[-I,]
dfr <- data.frame(dfr)
dfr$Score <- as.numeric(as.character(dfr$Score))
p<-ggplot(dfr, aes(x=Correction,y=Score, fill=CNV)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("-2","-1","0","1","2"),values=c("red", "green","blue","yellow","black"))

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
Beta <- Beta[,names(Subgroups)]
ResultsTFs <- matrix(0,ncol=length(unique(Subgroups)),nrow=nrow(Beta))
rownames(ResultsTFs) <- rownames(Beta)
colnames(ResultsTFs) <- unique(Subgroups)
for (i in (1:length(unique(Subgroups)))){
  I = names(Subgroups)[which(Subgroups==unique(Subgroups)[i])]
  subtype = rowMeans(Beta[,I])
  ResultsTFs[,i] <- subtype
}
total = rowMeans(Beta)
ResultsTFs <- cbind(ResultsTFs,total)
ResultsTFs <- ResultsTFs[order(total,decreasing=TRUE),]
TargetDirectory <- paste0(path,"results/")
write.csv2(ResultsTFs,file=paste0(TargetDirectory,"DeregulatedTFsCIT_Licorn.csv"))

#### expression diff ####
MA_cancer_norm <- MA_cancer[names(Subgroups),]
MA_cancer_norm = MA_cancer_norm - matrix(1,nrow(MA_cancer_norm),1) %*% colMeans(MA_cancer_norm)
MA_cancer_norm = MA_cancer_norm / sqrt( (1/nrow(MA_cancer_norm)) * matrix(1,nrow(MA_cancer_norm),1) %*% colSums(MA_cancer_norm^2))

library(limma)
Design <- matrix(0,ncol=length(unique(Subgroups)),nrow=nrow(MA_cancer_norm))
colnames(Design) <- unique(Subgroups)
rownames(Design) <- rownames(MA_cancer_norm)
for (i in (1:length(unique(Subgroups)))){
  I <- which(Subgroups==unique(Subgroups)[i])
  Design[I,i] <- rep(1,length(I))
}

# for each of the subtypes 
DEG <- list()
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
Beta <- Beta[,names(Subgroups)]
for (i in (1:length(unique(Subgroups)))){
  Design2 <- cbind(Design[,i],rowSums(Design[,-i]))
  #colnames(Design) <- c("Pap","Lum","Basal","tcga")
  fit <- lmFit(Beta,Design2)
  #constrat.matrix <- makeContrasts("Pap-Lum","Pap-Basal","Pap-tcga","Lum-Basal","Lum-tcga","Basal-tcga",levels=Design)
  #fit2 <- contrasts.fit(fit,constrat.matrix)
  fit <- eBayes(fit)
  Expressed <- topTable(fit,coef=1,sort.by = "logFC" ,number=nrow(Beta),p.value = 0.05,lfc = 0.1375)
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