##### In this code, we analyse the LIONS results ####
path <- "/Users/Magali/Desktop/recherche/LIONS/"
TargetDirectory <- paste0(path,"results/")

# load the results
load(paste0(path,"results/ResultsCIT_Licorn.Rdata"))
Beta <- Results$Beta

# set a threshold
threshold <- 0.01
Beta[which(Beta<threshold)] <- rep(0,length(which(Beta<threshold)))
Deregulated_TFs <- rowSums(Beta)[which(rowSums(Beta>0)>0)]
Beta <- Beta[names(Deregulated_TFs),]

# on regarde la répartition de Beta par sous-types
load(paste0(path,"Data/Subtypes/CIT_BLCA_Subgroup.RData"))
Subgroups <- CIT_BLCA_Subgroup
Samples <- sub('[A-Z]{3}$', "", names(Subgroups))
Subgroups <- as.character(Subgroups)
names(Subgroups) <- Samples
  
Beta <- Beta[,names(Subgroups)]
Deregulated_TFs <- rowSums(Beta)[which(rowSums(Beta>0)>0)]
Beta <- Beta[names(Deregulated_TFs),]

ResultsTFs <- matrix(0,ncol=length(unique(Subgroups)),nrow=nrow(Beta))
rownames(ResultsTFs) <- rownames(Beta)
colnames(ResultsTFs) <- unique(Subgroups)
for (i in (1:length(unique(Subgroups)))){
  I = names(Subgroups)[which(Subgroups==unique(Subgroups)[i])]
  subtype = rowSums(abs(Beta[,I])>0)/length(I)
  ResultsTFs[,i] <- subtype
}
total = rowSums(ResultsTFs)
ResultsTFs <- ResultsTFs[order(total,decreasing=TRUE),]
TargetDirectory <- paste0(path,"results/")
write.csv2(ResultsTFs,file=paste0(TargetDirectory,"DeregulatedTFsCIT_Licorn.csv"))

# define our own subtypes
library(ConsensusClusterPlus)
results <- ConsensusClusterPlus(Beta,maxK=20,reps=100,pItem=0.8,pFeature=1,distance="euclidian",clusterAlg="km",plot="png")
NrClusters <- 7
ClustResults <- results[[NrClusters]]
Subtypes_Beta <- ClustResults$consensusClass
Subtypes_Beta <- as.character(Subtypes_Beta)
names(Subtypes_Beta) <- names(ClustResults$consensusClass)

load(paste0(TargetDirectory,"NewSubtypes.Rdata"))

library(ComplexHeatmap)
library(circlize)
Beta_clus <- c()
Sub_clus <- c()
for (j in (1:length(unique(Subtypes_Beta)))){
  I <- names(Subtypes_Beta)[which(as.character(Subtypes_Beta)==unique(Subtypes_Beta)[j])]
  EnrichSmall <- Beta[,I]
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
heatmap <- Heatmap(Beta_clus, name = "Score", cluster_rows=TRUE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,show_row_names=TRUE,row_names_gp=gpar(col=c(rep("black",nrow(Beta_clus))),fontsize=8),
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Transcription factors",nrow(Beta_clus))),gap = unit(5, "mm"),
                   #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                   col=colorRamp2(c(0,0.1,max(Beta_clus)), c("white","aquamarine","blue")),heatmap_legend_param = list(color_bar = "continuous"),top_annotation=ha)
draw(heatmap)

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
                col= list(OurSubtypes=c("1"="#FF0000FF","2"="#FF3300FF","3"="#FF6600FF","4"= "#FF9900FF","5"="#FFCC00FF","6"="#FFFF00FF","7"="#FFFF80FF"),
                Subtypes=c("luminal-1"="blue","luminal-2"="green3","basal-like"="cyan","tcga-4"="purple")))
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
