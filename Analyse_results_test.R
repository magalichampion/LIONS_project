##### In this code, we analyse the LIONS results ####
path <- "/Users/Magali/Desktop/recherche/LIONS/"
TargetDirectory <- paste0(path,"results/")

# load the results
load(paste0(path,"results/Results2.Rdata"))
Beta <- Results$Beta
Beta <- Beta[!(rowSums(Beta) == 0),]

# on regarde la répartition de Beta par sous-types
load(paste0(path,"Data/Subtypes/CIT_BLCA_Subgroup.RData"))
Subgroups <- CIT_BLCA_Subgroup
Samples <- sub('[A-Z]{3}$', "", names(Subgroups))
Subgroups <- as.character(Subgroups)
names(Subgroups) <- Samples

Subgroups_full <- read.table(paste0(path,"Data/Subtypes/BLCA-TP.mergedcluster.csv"),sep=",",header=TRUE)
rownames(Subgroups_full) <- Subgroups_full$SampleName  
Subgroups <- Subgroups_full$mRNAseq_cHierarchical
names(Subgroups) <- rownames(Subgroups_full)
Subgroups <- Subgroups[colnames(Beta)]

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
write.csv2(ResultsTFs,file=paste0(TargetDirectory,"DeregulatedTFsTCGA_Licorn.csv"))

# define our own subtypes
library(ConsensusClusterPlus)
results <- ConsensusClusterPlus(Beta,maxK=20,reps=100,pItem=0.8,pFeature=1,distance="euclidian",clusterAlg="km",plot="png")
NrClusters <- 4
ClustResults <- results[[NrClusters]]
Subtypes_Beta <- ClustResults$consensusClass
Subtypes_Beta <- as.character(Subtypes_Beta)
names(Subtypes_Beta) <- names(ClustResults$consensusClass)

# another way to define clusters
d <-  dist(t(Beta),method="euclidian")
fit <- hclust(d,method="ward.D")
plot(fit)
groups <- cutree(fit,k=4)
rect.hclust(fit,k=4,border="red")

load(paste0(path,"Data/Subtypes/OurSubtypes.Rdata"))

library(ComplexHeatmap)
library(circlize)
Subtypes <- Subgroups # to fill
Beta_clus <- c()
Sub_clus <- c()
for (j in (1:length(unique(Subtypes)))){
  I <- names(Subtypes)[which(as.character(Subtypes)==unique(Subtypes)[j])]
  EnrichSmall <- Beta[,I]
  if (length(I)>1){
    SampleClustering <- hclust(dist(t(EnrichSmall)), method = "complete", members = NULL)    
    EnrichSmall <- EnrichSmall[,SampleClustering$order]
    Beta_clus <- cbind(Beta_clus,EnrichSmall)
    Sub_small <- Subtypes[I]
    Sub_clus <- c(Sub_clus,as.character(Sub_small[SampleClustering$order]))
  } else {
    Beta_clus <- cbind(Beta_clus,EnrichSmall)
    colnames(Beta_clus)[ncol(Beta_clus)] <- I
    Sub_clus <- c(Sub_clus,paste0("",j))
  }
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

# another clustering view
library(mixOmics)
cim(Beta)

# two annotations
Beta_clus <- c()
Sub_clus <- c()
for (j in (1:length(unique(Subtypes_Beta)))){
  I <- names(Subtypes_Beta)[which(as.character(Subtypes_Beta)==unique(Subtypes_Beta)[j])]
  EnrichSmall <- Beta[,I]
  if (length(I)>1){
    for (t in (1:length(unique(Subgroups)))){
      Sub <- Subgroups[intersect(names(Subtypes_Beta),I)]
      I2 <- names(Sub[which(Sub==unique(Subgroups)[t])])
      EnrichSmall2 <- EnrichSmall[,I2]
      
      if (length(I2)>1){
        SampleClustering <- hclust(dist(t(EnrichSmall2)), method = "complete", members = NULL)    
        EnrichSmall2 <- EnrichSmall2[,SampleClustering$order]
        Beta_clus <- cbind(Beta_clus,EnrichSmall2)
        Sub_small <- Subtypes_Beta[I2]
        Sub_clus <- c(Sub_clus,as.character(Sub_small[SampleClustering$order]))
      } else {
        if (length(I2)==1){
          Beta_clus <- cbind(Beta_clus,EnrichSmall2)
          colnames(Beta_clus)[ncol(Beta_clus)] <- I2
          Sub_clus <- c(Sub_clus,paste0("",j))
        }
      }
    }
  } else {
    Beta_clus <- cbind(Beta_clus,EnrichSmall)
    colnames(Beta_clus)[ncol(Beta_clus)] <- I
    Sub_clus <- c(Sub_clus,paste0("",j))
  }
}

Subtypes <- data.frame(OurSubtypes=Sub_clus,Subtypes=Subgroups[colnames(Beta_clus)])
ha2 <- HeatmapAnnotation(df=Subtypes, 
                col= list(OurSubtypes=c("1"="#FF0000FF","2"="#FF3300FF","3"="#FF6600FF","4"= "#FF9900FF","5"="#FFCC00FF"),
                Subtypes=c("1"="blue","2"="green3","3"="cyan")))
heatmap <- Heatmap(Beta_clus, name = "Score", cluster_rows=TRUE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,show_row_names=TRUE,row_names_gp=gpar(col=c(rep("black",nrow(Beta_clus))),fontsize=8),
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Transcription factors",nrow(Beta_clus))),gap = unit(5, "mm"),
                   #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                   col=colorRamp2(c(0,0.1,max(Beta_clus)), c("white","aquamarine","blue")),heatmap_legend_param = list(color_bar = "continuous"),top_annotation=ha2)
draw(heatmap)

threshold <- 0.01
Beta[which(Beta<threshold)] <- rep(0,length(which(Beta<threshold)))
Deregulated_TFs <- rowSums(Beta)[which(rowSums(Beta>0)>0)]
Beta <- Beta[names(Deregulated_TFs),]

ResultsTFs_Clust <- matrix(0,ncol=NrClusters,nrow=nrow(Beta))
rownames(ResultsTFs_Clust) <- rownames(Beta)
colnames(ResultsTFs_Clust) <- paste0("Cluster",c(1:7))
for (i in (1:NrClusters)){
  I = names(Subtypes_Beta)[which(Subtypes_Beta==i)]
  if (length(I)>1){
    subtype = rowSums(abs(Beta[,I])>0)/length(I)
  } else { 
    Nonnull <- which(abs(Beta[,I])>0)
    subtype <- rep(0,nrow(Beta))
    names(subtype) <- rownames(Beta)
    subtype[Nonnull] <- rep(1,length(Nonnull))
  }
  ResultsTFs_Clust[,i] <- subtype
}
total = rowSums(ResultsTFs_Clust)
ResultsTFs_Clust <- ResultsTFs_Clust[order(total,decreasing=TRUE),]
write.csv2(ResultsTFs_Clust,file=paste0(TargetDirectory,"DeregulatedTFsCIT_Licorn_OurSubtypes.csv"))

# code de Thomas  ############ DOES NOT WORK NOW
# load data (here we use the matrix beta as data)
Beta

# convert into boolean
Beta_boo <- abs(Beta)>0

# Convert boolean matrices into graphs (adjacency lists, actually)
graph.beta = matrix2graph(Beta_boo)

# We will remove the alterations that affect less than 2 samples
filterlevel=2
graph.beta_restricted = Filter(function(l) length(l)>=filterlevel,graph.beta)

#Compute the domination relationship among genes:
Subgroups_temp=Subtypes_Beta[colnames(Beta)]
reds = Subgroups_temp=="1" # try with this subtype
names(reds) = names(Subgroups_temp)
dom.genes=domination.graph(graph.beta_restricted,reds)

# For the analysis, we only keep the genes that are not dominated.
graph.beta.master=graph.beta_restricted[sapply(dom.genes$dominators,function(x) length(x)==0)]

# Another step in preparing the graph is merging the alterations that
# are "identical", i.e. affect exactly the same samples :
merged=merge.identical.nodes(graph.beta.master)
#look how they were grouped
#View(merged$family.data)

# the resulting graph is :
graph=merged$graph

# Now we invoke the core of this package : first see how the solution tree grows with the threshold
growth=treeGrowth(graph.beta.master,reds)

plot(growth$threshold,growth$nodes)
points(growth$threshold,growth$leaves,col='blue')

# The explotian seems to start at about thresold=0.0006
# To see  better, extract the convex bottom:
cb=convex.bottom(growth$threshold,growth$nodes)

lines(growth$threshold[cb],growth$nodes[cb],col='red')

# We choose to set
threshold = # to be set
  # 0.008502720 for luminal1
  # no significant results for luminal2
  # 0.0014251844 for tcga4
  # 0.001087472 for basal-like
  sol=solutions(graph,reds,threshold = threshold)

# and finally we can include pathway information in these results :
data("kegg_pathways")

sol=pathway_cover(sol,gene_pathways,pathwayname,merged$families)

### GSEA for some genes
library("GSA")
library("limma")

source("/Users/Magali/Desktop/recherche/PancancerAnalysis/amarettodevel/PostProcessing/EnrichmentAnalysis/EA_createData.R")
source("/Users/Magali/Desktop/recherche/PancancerAnalysis/amarettodevel/Postprocessing/EnrichmentAnalysis/EA_run.R")

DataBase <- CreateDataBase() 
GeneList <- c("NFKB1","HOXA10","ZNF567","TAF9B","ZNF420","MEIS2","SERTAD2","HOXB2","HOXD10","HOXC8","ZBTB10","TRIM29","ZNF23","RUVBL1","ZFAT","NOTCH2","CHuRC1","PITX3","ASB12","ZNF610")      # what you want
GeneList <- c("ZFP36L2","SOX15","SERTAD2","MEIS2","CHURC1","ZNF211","HOXA10","ATM","KLF9","BTBD11","ZNF655","RUVBL1","HOXA5","ZN22","ZNFX11","BATF")
GeneList <- c("SERTAD2","HOXA10","MEIS2")
GeneList <- c("MEIS2","HOXB2","NOTCH3","TRERF1","HOXD10","ETV6","ELK3","ELF3","FOXQ1","GRHL3","HOXA10","ZBTB10","GATA3","CHURC1","TRIM22","GRHL1","ZFP36L2","SOX15","ZNF655","JUN","NFKB1","NR2F6","PPARG")
GeneList <- c("FOXJ1","MEIS2","ETV6","GATA3","HOXC8","PITX3","NFKB1")
GeneList <- c("HOXA10","PLOR3F","KLF3","EOMES")
Results = EnrichmentAnalysis(Data=GeneList,DataBase=DataBase)

### Expression differentielle
MA <- read.table("/Users/Magali/Desktop/recherche/LIONS/results/Expression.txt")
load(paste0(path,"Data/Subtypes/CIT_BLCA_Subgroup.RData"))
Subgroups <- CIT_BLCA_Subgroup
names(Subgroups) <- str_replace_all(names(Subgroups),"[A-Z]","")
MA <- MA[names(Subgroups),intersect(TFs,colnames(MA))]

MAcancerNorm = sweep(MA, 2, colMeans(MA), "-")
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
ResultsTFs <- matrix(0,ncol=length(unique(Subgroups)),nrow=ncol(FoldChanges_reduced))
rownames(ResultsTFs) <- colnames(FoldChanges_reduced)
colnames(ResultsTFs) <- unique(Subgroups)
for (i in (1:length(unique(Subgroups)))){
  I = names(Subgroups)[which(Subgroups==unique(Subgroups)[i])]
  subtype = colSums(abs(FoldChanges_reduced[I,])>0)/length(I)
  ResultsTFs[,i] <- subtype
}
total = rowSums(ResultsTFs)
ResultsTFs <- ResultsTFs[order(total,decreasing=TRUE),]
write.csv2(ResultsTFs,file=paste0(path,"/results/DeregulatedTFs_FoldChanges.csv"))

library(limma)
Design <- matrix(0,ncol=length(unique(Subgroups)),nrow=nrow(MA_cancer))
colnames(Design) <- unique(Subgroups)
rownames(Design) <- rownames(MA_cancer)
for (i in (1:length(unique(Subgroups)))){
  I <- which(Subgroups==unique(Subgroups)[i])
  Design[I,i] <- rep(1,length(I))
}
Design2 <- cbind(Design[,3],Design[,1]+Design[,2])
fit <- lmFit(t(MA_cancer),Design2)
fit <- eBayes(fit)
topTable(fit,coef=2)

#### TCGA codes
CNV <- read.csv2(paste0(path,"Data/TCGA/CNV_TCGA.csv"),sep=",",header = TRUE,row.names=1)
samples <- intersect(rownames(CNV),colnames(Beta))
CNV <- t(CNV[samples,])
gene_names <- read.csv2(paste0(path,"Data/TCGA/gene_names.csv"),sep=",",header=FALSE)
gene_names <- as.matrix(gene_names)
I <- match(gene_names[,2],rownames(CNV))
rownames(CNV)[I] <- gene_names[,1]

CNV <- read.csv2(paste0(path,"Data/TCGA/blca_tcga_blca_tcga_gistic.csv"),sep=",",header=TRUE,row.names=1)
colnames(CNV) <- str_replace_all(colnames(CNV),"\\.","-")

load("/Users/Magali/Desktop/recherche/LIONS/Data/copyNumber/CIT/CIT_BLCA_CNA.rda")
colnames(CIT_BLCA_CNA) <- str_replace_all(colnames(CIT_BLCA_CNA),"[A-Z]","")
#Genes <- unique(str_replace_all(rownames(Scores_log),"_TF",""))
#Genes <- intersect(Genes,rownames(CIT_BLCA_CNA))
Genes <- intersect(rownames(CIT_BLCA_CNA),rownames(Beta))
CNV <- CIT_BLCA_CNA[Genes,intersect(colnames(Beta),colnames(CIT_BLCA_CNA))]
#I = rownames(Scores_log)[grep("_TF",rownames(Scores_log))]
#I <- str_replace_all(I,"_TF","")
#rownames(CNV)[which(rownames(CNV) %in% I)] <- paste0(rownames(CNV)[which(rownames(CNV) %in% I)],"_TF")
#Scores <- Scores_log[rownames(CNV),intersect(colnames(Scores_log),colnames(CIT_BLCA_CNA))]
Beta_short <- Beta[Genes,intersect(colnames(Beta),colnames(CIT_BLCA_CNA))]
CNV_transformed <- CNV
CNV_transformed <- as.matrix(CNV_transformed,nrow=nrow(CNV),ncol=ncol(CNV))
colnames(CNV_transformed) <- colnames(CNV)
rownames(CNV_transformed) <- rownames(CNV)

I_genes <- which(rowMeans(is.na(CNV_transformed))>0.1)
I_samples <- which(colMeans(is.na(CNV_transformed))>0.1)
if (length(I_genes)>0){
  CNV_transformed <- CNV_transformed[-I_genes,]
}
if (length(I_samples)>0){
  CNV_transformed <- CNV_transformed[,-I_samples]
}
CNV_transformed <- CNV_transformed[intersect(rownames(Beta),rownames(CNV_transformed)),intersect(colnames(Beta),colnames(CNV_transformed))]
Beta_short <- Beta[intersect(rownames(Beta),rownames(CNV_transformed)),intersect(colnames(Beta),colnames(CNV_transformed))]

# transfo 
CNV_transformed[which(abs(CNV_transformed)<2)] <- rep("0",length(which(abs(CNV_transformed)<2)))
CNV_transformed[which(CNV_transformed=="2")] <- rep("1",length(which(CNV_transformed=="2")))
CNV_transformed[which(CNV_transformed=="-2")] <- rep("1",length(which(CNV_transformed=="-2")))

CNV_transformed[which(abs(CNV_transformed)>0)] <- rep("1",length(which(abs(CNV_transformed)>0)))

# threshold for genes
# extract amplified and deleted genes
AMPfile = paste0(path,'Data/TCGA/amp_genes.conf_99.txt')
AMPtable=read.csv(AMPfile,sep="\t",head=TRUE,na.strings="NA") 
AMPgenes=c()
if (ncol(AMPtable)>1){
  for (i in 2:ncol(AMPtable)) {
    tmpGenes=as.vector(AMPtable[4:nrow(AMPtable),i])
    tmpGenes[tmpGenes==""]=NA
    tmpGenes=na.omit(tmpGenes)
    AMPgenes=c(AMPgenes,tmpGenes)         
  }
  AMPgenes=unique(AMPgenes)
  AMPgenes=AMPgenes[order(AMPgenes)]
} 
AMPgenes <- intersect(AMPgenes,rownames(Beta))
# Same for DEL genes
DELfile=paste0(path,'Data/TCGA/del_genes.conf_99.txt')
DELtable=read.csv(DELfile,sep="\t",head=TRUE,na.strings="NA")     
DELgenes=c()
if (ncol(DELtable)>1){
  for (i in 2:ncol(DELtable)) {
    tmpGenes=as.vector(DELtable[4:nrow(DELtable),i])
    tmpGenes[tmpGenes==""]=NA
    tmpGenes=na.omit(tmpGenes)
    DELgenes=c(DELgenes,tmpGenes)         
  }
  DELgenes=unique(DELgenes)
  DELgenes=DELgenes[order(DELgenes)]
}
DELgenes <- intersect(DELgenes,rownames(Beta))
AlteredGenes <- c(AMPgenes,DELgenes)

# filter genes and samples
CNVgenes <- names(which(rowMeans(CNV_transformed=="1",na.rm = TRUE)>0.10))
CNVgenes <- rownames(CNV_transformed)
#Beta_genes <- names(which(rowMeans(Beta)>0.1))

#Varying <- rowVars(Beta)
#names(Varying) <- rownames(Beta)
#Beta_genes <- names(Varying[which(Varying>threshold)])

# code for p-values
TopGenes <- CNVgenes
No <- c()
PVal_Genes <- c()
for (i in (1:length(TopGenes))){
  Gene <- TopGenes[i]
  if (dim(table(CNV_transformed[Gene,]))==1){
    No <- c(No,Gene)
  } else {
    I <- which(CNV_transformed[Gene,]=="NaN")
    if (length(I)>0){
      if (length(I)==ncol(Beta_short)){
        # no data
        No <- c(No,Gene)
      } else {
        C <- wilcox.test(Beta_short[Gene,-I]~CNV_transformed[Gene,-I])
        PVal_Genes <- c(PVal_Genes,C$p.value)
      } 
    } else {
      C <- wilcox.test(Beta_short[Gene,]~CNV_transformed[Gene,])
      PVal_Genes <- c(PVal_Genes,C$p.value)
    }
  }
}
names(PVal_Genes) <- TopGenes[which(is.na(match(TopGenes,No))==TRUE)]
PVal_Genes <- PVal_Genes[order(PVal_Genes,decreasing=FALSE)]
Adj_pval <- p.adjust(PVal_Genes,method="BH")

TF_CNV <- names(PVal_Genes)[which(p.adjust(PVal_Genes,method="BH")<0.05)]

No <- c()
PVal_Genes2 <- c()
for (i in (1:length(TopGenes))){
  Gene <- TopGenes[i]
  if (dim(table(CNV_transformed[Gene,]))==1){
    No <- c(No,Gene)
  } else {
    I <- which(CNV_transformed[Gene,]=="NaN")
    if (length(I)>0){
      if (length(I)==ncol(Beta_short)){
        # no data
        No <- c(No,Gene)
      } else {
        C <- kruskal.test(Beta_short[Gene,-I]~CNV_transformed[Gene,-I])
        PVal_Genes2 <- c(PVal_Genes2,C$p.value)
      } 
    } else {
      C <- kruskal.test(Beta_short[Gene,]~CNV_transformed[Gene,])
      PVal_Genes2 <- c(PVal_Genes2,C$p.value)
    }
  }
}
names(PVal_Genes2) <- TopGenes[which(is.na(match(TopGenes,No))==TRUE)]
PVal_Genes2 <- PVal_Genes2[!is.na(PVal_Genes2)]
Adjust <- p.adjust(PVal_Genes2,method="BH")
Adjust <- Adjust[order(Adjust,decreasing = FALSE)]
TF_CNV <- names(Adjust)[which(Adjust<0.05)]

TopGenes <- CNVgenes
No <- c()
PVal_Genes <- c()
for (i in (1:length(TopGenes))){
  Gene <- TopGenes[i]
  if (dim(table(CNV_transformed[Gene,]))==1){
    No <- c(No,Gene)
  } else {
    I <- which(CNV_transformed[Gene,]=="NaN")
    I0 <- which(Beta_short[Gene,]==0)
    if (length(I)>0){
      if (length(I)==ncol(Beta_short)){
        # no data
        No <- c(No,Gene)
      } else {
        if (dim(table(CNV_transformed[Gene,-c(I,I0)]))==1){
          No <- c(No,Gene)
        } else {
          C <- kruskal.test(Beta_short[Gene,-c(I,I0)]~CNV_transformed[Gene,-c(I,I0)])
          PVal_Genes <- c(PVal_Genes,C$p.value)
        }
      } 
    } else {
      if (dim(table(CNV_transformed[Gene,-I0]))<2){
        No <- c(No,Gene)
      } else {
        C <- kruskal.test(Beta_short[Gene,-I0]~CNV_transformed[Gene,-I0])
        PVal_Genes <- c(PVal_Genes,C$p.value)
      }
    }
  }
}
names(PVal_Genes) <- TopGenes[which(is.na(match(TopGenes,No))==TRUE)]
PVal_Genes <- PVal_Genes[order(PVal_Genes,decreasing=FALSE)]
Adj_pval <- p.adjust(PVal_Genes,method="BH")

TF_CNV <- names(PVal_Genes)[which(p.adjust(PVal_Genes,method="BH")<0.05)]

#corr_Samples <- c()
#for (i in (1:ncol(Beta))){
#  corr_Samples <- c(corr_Samples,cor(Beta[i,],CNV_transformed[i,],method="pearson",use = "pairwise.complete.obs"))
#}
#names(corr_Samples) <- colnames(Beta)
#corr_Genes <- corr_Genes[order(abs(corr_Genes),decreasing=TRUE)]
#corr_Samples <- corr_Samples[order(abs(corr_Samples),decreasing=TRUE)]

library(circlize)
library(ComplexHeatmap)
library(purrr)
library(dplyr)

Beta_small <- Beta_short[TF_CNV,]
Annot <- t(CNV_transformed[TF_CNV,])
Annot[which(is.na(Annot)==TRUE)] <- rep("0",length(which(is.na(Annot)==TRUE)))
Gene <- "ASH2L"
Beta_clus <- c()
annot <- c()
for (i in (1:length(unique(Annot[,Gene])))){
  I <- names(which(Annot[,Gene]==unique(Annot[,Gene])[i]))
  Beta_clus <- cbind(Beta_clus,Beta_small[,I])  
  annot <- c(annot,Annot[I,Gene])
}
annot <- data.frame(annot)
ha = HeatmapAnnotation(df = annot)
Heatmap(Beta_clus,cluster_columns = FALSE,cluster_rows = FALSE,show_column_dend=FALSE,show_column_names=FALSE,show_row_names=TRUE,col=colorRamp2(c(0,0.1,max(Beta_small)), c("white","aquamarine","blue")),heatmap_legend_param = list(color_bar = "continuous"),top_annotation=ha)

for (i in (1:length(TF_CNV))){
  boxplot(Beta_short[TF_CNV[i],]~CNV_transformed[TF_CNV[i],],main=paste0("TF ",TF_CNV[i]))
}
dfr <- matrix(0,ncol=3,nrow=length(TF_CNV)*ncol(Beta_short))
colnames(dfr) <- c("Beta","Gene","CNV")
for (i in (1:length(TF_CNV))){
  gene <- TF_CNV[i]
  dfr[(ncol(Beta_short)*(i-1)+1):(ncol(Beta_short)*i),1] <- Beta_short[gene,]
  dfr[(ncol(Beta_short)*(i-1)+1):(ncol(Beta_short)*i),2] <- rep(gene,ncol(Beta_short))
  dfr[(ncol(Beta_short)*(i-1)+1):(ncol(Beta_short)*i),3] <- CNV_transformed[gene,]
}
dfr <- data.frame(dfr)
dfr$Beta <- as.numeric(as.character(dfr$Beta))
p<-ggplot(dfr, aes(x=Gene,y=Beta, fill=CNV)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("-2","-1","0","1","2"),values=c("red", "green","blue","yellow","black"))

#### comparer avec et sans correction
load(paste0(path,"results/ResultsTCGA_Licorn.Rdata"))
Beta <- Results$Beta
Beta <- Beta[intersect(rownames(Beta),rownames(CNV_transformed)),intersect(colnames(Beta),colnames(CNV_transformed))]

load(paste0(path,"results/ResultsTCGA_Licorn_Corrected.Rdata"))
Beta_corrected <- Results$Beta
Beta_corrected <- Beta_corrected[intersect(rownames(Beta_corrected),rownames(CNV_transformed)),intersect(colnames(Beta_corrected),colnames(CNV_transformed))]

Ecart <- Beta - Beta_corrected
PVal <- c()
for (i in (1:nrow(Beta))){
  gene <- rownames(Beta)[i]
  cnv <- rep("Nothing",ncol(Beta))
  cnv[which(abs(CNV_transformed[gene,])>1)] <- rep("Altered",length(which(abs(CNV_transformed[gene,])>1)))
  if (length(table(cnv))==1){
    pval <- 1
    PVal <- c(PVal,pval)
  } else {
    pval <- wilcox.test(Ecart[gene,]~cnv)
    PVal <- c(PVal,pval$p.value)
  }
}
names(PVal) = rownames(Beta)
PVal <- PVal[order(PVal)] 

# reassocier les sous-types aux marqueurs de chaque sous-types
Genes <- c("CYP2J2","FGFR3","ERBB2","ERBB3","FOXA1","GATA3","GPX2","KRT18","KRT19","KRT20","KRT7","KRT8","PPARG","XBP1","CD44","CDH3","KRT1","KRT14","KRT16","KRT5","KRT6A","KRT6B","KRT6C","ACTG2","CNN1","MYH11","MFAP4","PGM5","FLNC","ACTC1","DES","PCP4")
Genes <- intersect(Genes,colnames(MA_cancer))
MA <- t(MA_cancer[,Genes])

library(ComplexHeatmap)
library(circlize)
Beta_clus <- c()
Sub_clus <- c()
for (j in (1:length(unique(Subtypes_Beta)))){
  I <- names(Subtypes_Beta)[which(as.character(Subtypes_Beta)==unique(Subtypes_Beta)[j])]
  EnrichSmall <- MA[,I]
  if (length(I)>1){
    SampleClustering <- hclust(dist(t(EnrichSmall)), method = "complete", members = NULL)    
    EnrichSmall <- EnrichSmall[,SampleClustering$order]
    Beta_clus <- cbind(Beta_clus,EnrichSmall)
    Sub_small <- Subtypes_Beta[I]
    Sub_clus <- c(Sub_clus,as.character(Sub_small[SampleClustering$order]))
  } else {
    Beta_clus <- cbind(Beta_clus,EnrichSmall)
    colnames(Beta_clus)[ncol(Beta_clus)] <- I
    Sub_clus <- c(Sub_clus,paste0("",j))
  }
}
#GeneClustering <- hclust(dist(Beta_clus), method = "complete", members = NULL)
#Beta_clus <- Beta_clus[GeneClustering$order,]
Subtypes <- data.frame(Sub_clus)

ha = HeatmapAnnotation(df = Subtypes)
heatmap <- Heatmap(Beta_clus, name = "Score", cluster_rows=FALSE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,show_row_names=TRUE,row_names_gp=gpar(col=c(rep("black",nrow(Beta_clus))),fontsize=8),
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Transcription factors",nrow(Beta_clus))),gap = unit(5, "mm"),
                   #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                   heatmap_legend_param = list(color_bar = "continuous"),top_annotation=ha)
draw(heatmap)

# distributions according to the CNV status
par(mfrow = c(1, 1))
I0 <- which(Beta_short==0)
Beta_very_small <- Beta_short[-I0]
CNV_very_small <- CNV_transformed[-I0]
I2 <- which(CNV_very_small==2)
F2 <- ecdf(Beta_very_small[I2])
plot(F2,verticals=TRUE,do.points=FALSE,main="Cumulative empirical distribution",xlab="Score of deregulation",ylab="")
I1 <- which(CNV_very_small==1)
F1 <- ecdf(Beta_very_small[I1])
lines(F1,verticals=TRUE,do.points=FALSE,col="green")
I1m <- which(CNV_very_small==-1)
F1bis <- ecdf(Beta_very_small[I1m])
lines(F1bis,verticals=TRUE,do.points=FALSE,col="red")
I2m <- which(CNV_very_small==-2)
F2bis <- ecdf(Beta_very_small[I2m])
lines(F2bis,verticals=TRUE,do.points=FALSE,col="yellow")
I <- which(CNV_very_small==0)
F0 <- ecdf(Beta_very_small[I])
lines(F0,verticals=TRUE,do.points=FALSE,col="blue")

C <- wilcox.test(Beta_very_small[c(I,I2)]~CNV_very_small[c(I,I2)])

# targets of ASHL2
library(stringr)
library(ggplot2)
Targets <- c("BRF2","DDHD2","LSM1","PROSC","WHSC1L1_TF")
score_targets <- Scores_log[Targets,colnames(Beta)]
CNV <- read.csv2(paste0(path,"results/ASH2LGistic.csv"),sep=",",header=TRUE,row.names=1)
colnames(CNV) <- str_replace_all(colnames(CNV),"\\.","-")
CNV_Targets <- CNV
CNV_Targets <- as.matrix(CNV_Targets,nrow=nrow(CNV),ncol=ncol(CNV))
colnames(CNV_Targets) <- colnames(CNV)
rownames(CNV_Targets) <- rownames(CNV)
I_genes <- which(rowMeans(is.na(CNV_Targets))>0.1)
I_samples <- which(colMeans(is.na(CNV_Targets))>0.1)
if (length(I_genes)>0){
  CNV_Targets <- CNV_Targets[-I_genes,]
}
if (length(I_samples)>0){
  CNV_Targets <- CNV_Targets[,-I_samples]
}
CNV_Targets <- CNV_Targets[intersect(rownames(score_targets),rownames(CNV_Targets)),intersect(colnames(score_targets),colnames(CNV_Targets))]
score_targets <- score_targets[intersect(rownames(score_targets),rownames(CNV_Targets)),intersect(colnames(score_targets),colnames(CNV_Targets))]

CNV_Targets[which(abs(CNV_Targets)>0)] <- rep("1",length(which(abs(CNV_Targets)>0)))
dfr2 <- matrix(0,ncol=3,nrow=length(Targets)*ncol(score_targets))
colnames(dfr2) <- c("score","Gene","CNV")
for (i in (1:length(Targets))){
  gene <- Targets[i]
  dfr2[(ncol(score_targets)*(i-1)+1):(ncol(score_targets)*i),1] <- score_targets[gene,]
  dfr2[(ncol(score_targets)*(i-1)+1):(ncol(score_targets)*i),2] <- rep(gene,ncol(score_targets))
  dfr2[(ncol(score_targets)*(i-1)+1):(ncol(score_targets)*i),3] <- CNV_Targets[gene,]
}
dfr2 <- data.frame(dfr2)
dfr2$score <- as.numeric(as.character(dfr2$score))
p<-ggplot(dfr2, aes(x=Gene,y=score, fill=CNV)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("0","1"),values=c("red", "green"))

# same analysis with mutations
Mutations <- read.csv2(paste0(path,"Data/TCGA/blca_tcga_blca_tcga_mutations.csv"),sep=",",header=TRUE,row.names=1)
colnames(Mutations) <- str_replace_all(colnames(Mutations),"\\.","-")
Mutations_transformed <- Mutations
Mutations_transformed <- as.matrix(Mutations_transformed,nrow=nrow(Mutations),ncol=ncol(Mutations))
Mutations_transformed[which(Mutations_transformed=="NaN")] <- rep("0",length(which(Mutations_transformed=="NaN")))
Mutations_transformed[which(is.na(Mutations_transformed)==TRUE)] <- rep("0",length(which(is.na(Mutations_transformed)==TRUE)))

Genes <- intersect(rownames(Beta),rownames(Mutations_transformed))
Mutations_transformed <- Mutations_transformed[intersect(rownames(Beta),rownames(Mutations_transformed)),intersect(colnames(Beta),colnames(Mutations_transformed))]

I <- which(Mutations_transformed=="0")
Tot <- c(1:length(Mutations_transformed))
Tot <- Tot[-I]
Mutations_transformed[Tot] <- rep("1",length(Tot))
Beta <- Beta[intersect(rownames(Beta),rownames(Mutations_transformed)),intersect(colnames(Beta),colnames(Mutations_transformed))]

mutatedgenes <- names(which(rowMeans(Mutations_transformed=="0")<0.99))

# code for p-values
TopGenes <- mutatedgenes
No <- c()
PVal_Genes <- c()
for (i in (1:length(TopGenes))){
  Gene <- TopGenes[i]
  if (dim(table(Mutations_transformed[Gene,]))==1){
    No <- c(No,Gene)
  } else {
    I <- which(Mutations_transformed[Gene,]=="NaN")
    if (length(I)>0){
      if (length(I)==ncol(Beta)){
        # no data
        No <- c(No,Gene)
      } else {
        C <- wilcox.test(Beta[Gene,-I]~Mutations_transformed[Gene,-I])
        PVal_Genes <- c(PVal_Genes,C$p.value)
      } 
    } else {
      C <- wilcox.test(Beta[Gene,]~Mutations_transformed[Gene,])
      PVal_Genes <- c(PVal_Genes,C$p.value)
    }
  }
}
names(PVal_Genes) <- TopGenes[which(is.na(match(TopGenes,No))==TRUE)]
PVal_Genes <- PVal_Genes[order(PVal_Genes,decreasing=FALSE)]
Adj_pval <- p.adjust(PVal_Genes,method="fdr")

TF_mut <- names(PVal_Genes)[which(PVal_Genes<0.05)]
dfr <- matrix(0,ncol=3,nrow=length(TF_mut)*ncol(Beta))
colnames(dfr) <- c("Beta","Gene","mutation")
for (i in (1:length(TF_mut))){
  gene <- TF_mut[i]
  dfr[(ncol(Beta)*(i-1)+1):(ncol(Beta)*i),1] <- Beta[gene,]
  dfr[(ncol(Beta)*(i-1)+1):(ncol(Beta)*i),2] <- rep(gene,ncol(Beta))
  dfr[(ncol(Beta)*(i-1)+1):(ncol(Beta)*i),3] <- Mutations_transformed[gene,]
}
dfr <- data.frame(dfr)
dfr$Beta <- as.numeric(as.character(dfr$Beta))
p<-ggplot(dfr, aes(x=Gene,y=Beta, fill=mutation)) +
  geom_boxplot(position=position_dodge(1))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 72, vjust = 1, size = 12, hjust = 1))
p + scale_fill_manual(breaks=c("0","1"),values=c("red", "green"))


#### try if glmnet is better than lsei
fit <- glmnet(Adj_matrix,Scores_log[,1],lambda=0,upper.limits=1,lower.limits=0) # rk, upper.limits ne change rien
Beta_glmnet <- as.vector(fit$beta)
names(Beta_glmnet) <- colnames(Adj_matrix)

e <- rbind(diag(ncol(Adj_matrix)),-diag(ncol(Adj_matrix)))
f <- c(rep(0,ncol(Adj_matrix)),rep(-1,ncol(Adj_matrix)))
Beta_lsei <- lsei(a=Adj_matrix,b=Scores_log[,1],e=e,f=f)

Beta <- matrix(0,ncol=ncol(Score),nrow=ncol(Adj_matrix))
for (i in (1:ncol(Score))){
  fit <- glmnet(Adj_matrix,log(Score[,i]),lambda=0,upper.limits=1,lower.limits=0)
  beta <- as.vector(fit$beta)
  names(beta) <- colnames(Adj_matrix)
  Beta[,i] <- beta
}
colnames(Beta) <- colnames(Score)
rownames(Beta) <- colnames(Adj_matrix)

# keep genes with significant CNV on expression
Input <- Score  
rownames(Input) <- str_replace_all(rownames(Input),"_TF","")
Genes <- intersect(rownames(Input),rownames(CIT_BLCA_CNA))
Samples <- intersect(colnames(Input),colnames(CIT_BLCA_CNA))
CNV <- as.matrix(CIT_BLCA_CNA[Genes,Samples])
MA_matrix <- as.matrix(t(MA_cancer))[Genes,Samples]

CNVdrivers=c()
pval <- c()
for (i in 1:length(rownames(CNV))) {
  #res=lm(MA_matrix[i,]~CNV[i,])
  #res.summary=summary(res)
  #if (res$coefficients[2]>0 & res.summary$coefficients[2,4]<0.05 & res.summary$r.squared>0.25) {
  #  CNVdrivers=c(CNVdrivers,rownames(CNV)[i])
  #}
  test = kruskal.test(MA_matrix[i,]~CNV[i,])
  pval <- c(pval,test$p.value)
}
pval <- p.adjust(pval,method="fdr")
names(pval) <- rownames(CNV)
CNVdrivers <- names(which(pval<0.05)) 

# correction
MA_cancer_TFs <- matrix(0,ncol=(nrow(Adj_matrix)+ncol(Adj_matrix)),nrow=nrow(MA_cancer))
dim(MA_cancer_TFs)
rownames(MA_cancer_TFs) = rownames(MA_cancer)
colnames(MA_cancer_TFs) <- c(colnames(Adj_matrix),rownames(Adj_matrix))
InterTargets <- intersect(colnames(MA_cancer_TFs),rownames(Adj_matrix))
length(InterTargets)
InterTFs <- intersect(colnames(MA_cancer_TFs),colnames(Adj_matrix))
length(InterTFs)
MA_cancer_TFs[,InterTFs] <- as.matrix(MA_cancer[,InterTFs],nrow(MA_cancer),length(InterTFs))
InterTargets_TFs <- InterTargets[grep("_TF",InterTargets)]
length(InterTargets_TFs)
MA_cancer_TFs[,InterTargets_TFs] <- as.matrix(MA_matrix_corrected[,str_replace_all(InterTargets_TFs,"_TF","")],nrow(MA_cancer),length(InterTargets_TFs))
MA_cancer_TFs[,intersect(InterTargets,colnames(MA_cancer))] <- as.matrix(MA_matrix_corrected[,intersect(InterTargets,colnames(MA_cancer))],nrow(MA_cancer),length(InterTargets))
MA_cancer_TFs[1:10,1:10]
write.table(MA_cancer_TFs,file=paste0(TargetDirectory,"Expression2_partially_CorrectedBis.txt"),sep=" ",quote=FALSE)
