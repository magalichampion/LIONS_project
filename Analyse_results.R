##### In this code, we analyse the LIONS results ####
path <- "/Users/Magali/Desktop/recherche/LIONS/"

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
ha2 <- HeatmapAnnotation(df=data.frame(Subgroups[colnames(Beta_clus)]))
heatmap <- Heatmap(Beta_clus, name = "Score", cluster_rows=TRUE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,show_row_names=TRUE,row_names_gp=gpar(col=c(rep("black",nrow(Beta_clus))),fontsize=8),
                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Transcription factors",nrow(Beta_clus))),gap = unit(5, "mm"),
                   #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                   col=colorRamp2(c(0,0.1,max(Beta_clus)), c("white","aquamarine","blue")),heatmap_legend_param = list(color_bar = "continuous"),top_annotation=ha)
draw(heatmap)

ResultsTFs_Clust <- matrix(0,ncol=NrClusters,nrow=nrow(Beta))
rownames(ResultsTFs_Clust) <- rownames(Beta)
colnames(ResultsTFs_Clust) <- paste0("Cluster",c(1:6))
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
