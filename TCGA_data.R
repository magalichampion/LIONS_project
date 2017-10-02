#### Load and Process TCGA data ####
#### These codes were created by O. Gevaert (Stanford) ####

##############################################################################################################################
# dataType can be stddata or analyses
# dataFileTag options for stddata are: mRNAseq_Preprocess.Level_3, 
# dataFileTag options for analyses are: CopyNumber_Gistic2.Level_4, Pathway_Paradigm_RNASeq_And_Copy_Number.Level_4, 
#                   Pathway_Paradigm_mRNA_And_Copy_Number.Level_4 and MutSigNozzleReport2.0.Level_4
# TCGA_acronym_uppercase: BRCA, LUAD, COADREAD, etc. 
##############################################################################################################################

##############################################################################################################################
#code for testing this script
##############################################################################################################################
# gdacURL= "http://gdac.broadinstitute.org/runs/"
# TCGA_acronym_uppercase="COADREAD"
# fileType= "tar.gz"
# downloadData=TRUE
# dataType="stddata"
# saveDir = "./"
# tmpData=get_firehoseData(downloadData,saveDir,TCGA_acronym_uppercase)
###############################################################################################################################

Download_CancerSite <- function(CancerSite,TargetDirectory,downloadData=TRUE) {    
  cat('Downloading MA and CNV data for:',CancerSite,'\n')
  Cancers=c('ACC','BLCA','BRCA','CESC','CHOL','COAD','COADREAD','ESCA','GBM','GBMLGG','HNSC','KICH','KIPAN','KIRC','KIRP','LAML',
            'LGG','LIHC','LUAD','LUSC','OV','PAAD','PCPG','READ','SARC','STAD','THCA','THYM','UCEC')
  if (!(CancerSite %in% Cancers)) {
    cat('This TCGA cancer site/type was not tested, continue at your own risk.\n')
  }
  
  # Creating the top level directory where all data will be stored.
  if (TargetDirectory=="./"){
    TargetDirectory=paste0(getwd(),"/")
  }
  dir.create(TargetDirectory,showWarnings=FALSE)
  
  # Settings
  TCGA_acronym_uppercase=toupper(CancerSite)
  
  # get RNA seq data (GBM does not have much RNAseq data.)
  dataType='stddata'    
  dataFileTag='mRNAseq_Preprocess.Level_3'     
  
  #special case for GBM and OV, not enough RNAseq data, so using the microarray data instead
  if (CancerSite=="GBM") {                  
    dataFileTag=c('Merge_transcriptome__agilentg4502a_07_1__unc_edu__Level_3__unc_lowess_normalization_gene_level__data','Merge_transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data')                     
  } else if(CancerSite=="OV") {                   
    dataFileTag='Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data'        
  }    
  cat('Searching MA data for:',CancerSite,"\n")
  if (length(dataFileTag)==1) {      
    MAdirectory=get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataFileTag=dataFileTag)        
  } else {        
    MAdirectory=c()      
    for (i in 1:length(dataFileTag)) {
      MAdirectory=c(MAdirectory,get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataFileTag=dataFileTag[i]))     
    }        
  }
  
  # get CNV GISTIC data.
  dataType='analyses'
  dataFileTag='CopyNumber_Gistic2.Level_4'
  cat('Searching CNV data for:',CancerSite,'\n')
  CNVdirectory=get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataType=dataType,dataFileTag=dataFileTag)    
  return(list(MAdirectory=MAdirectory,CNVdirectory=CNVdirectory))
}

get_firehoseData <- function(downloadData=TRUE,saveDir = "./",TCGA_acronym_uppercase = "LUAD",dataType="stddata",dataFileTag = "mRNAseq_Preprocess.Level_3",
                             FFPE=FALSE,fileType= "tar.gz",  gdacURL= "http://gdac.broadinstitute.org/runs/",untarUngzip=TRUE,printDisease_abbr=FALSE){  
  
  # Cases Shipped by BCR  # Cases with Data*  Date Last Updated (mm/dd/yy)
  cancers <- c("Adrenocortical carcinoma [ACC]	\n","Bladder Urothelial Carcinoma [BLCA] \n",
               "Breast invasive carcinoma [BRCA] \n","Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] \n",
               "Cholangiocarcinoma [CHOL] \n",	"Colon adenocarcinoma [COAD] \n",
               "Colorectal adenocarcinoma [COADREAD] \n","Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]	\n",
               "Esophageal carcinoma [ESCA] \n","Glioblastoma multiforme [GBM] \n","Glioma [GBMLGG] \n",
               "Head and Neck squamous cell carcinoma [HNSC]	\n","Kidney Chromophobe [KICH]	\n","Pan-kidney cohort [KIPAN] n",
               "Kidney renal clear cell carcinoma [KIRC]	\n","Kidney renal papillary cell carcinoma [KIRP]	\n",
               "Acute Myeloid Leukemia [LAML] \n","Brain Lower Grade Glioma [LGG] \n",
               "Liver hepatocellular carcinoma [LIHC]	\n","Lung adenocarcinoma [LUAD]	\n",
               "Lung squamous cell carcinoma [LUSC] \n","Mesothelioma [MESO] \n",
               "Ovarian serous cystadenocarcinoma [OV]	\n","Pancreatic adenocarcinoma [PAAD]	\n",
               "Pheochromocytoma and Paraganglioma [PCPG] \n","Prostate adenocarcinoma [PRAD] \n",
               "Rectum adenocarcinoma [READ]	\n","Sarcoma [SARC]	\n","Skin Cutaneous Melanoma [SKCM]	\n",
               "Stomach adenocarcinoma [STAD] \n","Testicular Germ Cell Tumors [TGCT] \n","Thyroid carcinoma [THCA]	\n",
               "Thymoma [THYM] \n","Uterine Corpus Endometrial Carcinoma [UCEC]	\n","Uterine Carcinosarcoma [UCS]	 \n",
               "Uveal Melanoma [UVM] \n");
  
  cancers_acronyms <-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","COADREAD","DLBC","ESCA","GBM","GBMLGG","HNSC","KICH","KIPAN","KIRC",
                       "KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
                       "SKCM","STAD","STES","TGCT","THCA","THYM","UCEC","UCS","UVM")
  
  if(printDisease_abbr){      
    return(cat("here are the possible TCGA database disease acronyms. \nRe-run this function with printDisease_abbr=FALSE to then run an actual query.\n\n",cancers));  
  }
  if (TCGA_acronym_uppercase %in% cancers_acronyms){
    gdacURL_orig <- gdacURL
    urlData <- getURL(gdacURL)
    urlData <- strsplit2(urlData,paste(dataType,"__",sep=""))
    
    #first column is junk
    urlData <- urlData[,2:dim(urlData)[2]]
    urlData <- strsplit2(urlData,"/")
    urlData <- urlData[,1]
    
    # some of these data are secured
    if (length(grep(".bak",urlData))>0){
      I = grep(".bak",urlData)
      urlData[I] <- rep("Secured",length(I))
    }
    #see the strptime codes here: http://stat.ethz.ch/R-manual/R-devel/library/base/html/strptime.html
    #I am transferring this text to dates: that way, we can programmatically find the latest one.
    #the POSIXct class in R is the date-time class
    urlData <- as.POSIXct(strptime(urlData, "%Y_%m_%d"),tz="America/Los_Angeles")
    #want to remove time zone: do as.Date
    dateData <- as.Date(urlData[which(!is.na(urlData))])
    lastDate <- dateData[match( summary(dateData)[which(names(summary(dateData))=="Max.")], dateData)]
    
    #sub back _ symbols
    lastDate <- gsub("-","_",as.character(lastDate))
    
    lastDateCompress <- gsub("_","",lastDate)
    #need last "/" for it to find this page.
    
    gdacURL <- paste(gdacURL,dataType,"__",lastDate,"/data/",TCGA_acronym_uppercase,"/",lastDateCompress,"/",sep="")
    
    #now get full dataset name. we just have the http link to the last page we select this dataset from right now.
    urlData <- getURL(gdacURL)
    
    #regular expressions: need \ to have R recognize any " or \ that's actually in our text
    urlData <- strsplit2(urlData,"href=\\\"")
    
    while(length(grep("was not found",urlData))>0) {
      cat(file.path("\tNOTE: the TCGA run dated ",lastDate, "for ", dataType," for disease ",TCGA_acronym_uppercase," isn't available     for download yet.\n"))
      cat("\tTaking the run dated just before this one.\n")
      dateData <-  dateData[-which(dateData==(summary(dateData)[which(names(summary(dateData))=="Max.")]))]
      lastDate <- dateData[match( summary(dateData)[which(names(summary(dateData))=="Max.")], dateData)]
      #sub back _ symbols
      lastDate <- gsub("-","_",as.character(lastDate))
      lastDateCompress <- gsub("_","",lastDate)
      #need last "/" for it to find this page.
      gdacURL <- paste(gdacURL_orig,dataType,"__",lastDate,"/data/",TCGA_acronym_uppercase,"/",lastDateCompress,"/",sep="")
      
      #now get full dataset name. we just have the http link to the last page we select this dataset from right now.
      urlData <- getURL(gdacURL)
      #regular expressions: need \ to have R recognize any " or \ that's actually in our text
      urlData <- strsplit2(urlData,"href=\\\"")
      
      #did we reach the end of the dates - ie only one left? leave the loop then.
      if(length(dateData)<=1){        
        break    
      }  
    } 
    
    #cat("Using data from date ",lastDate,"\n")
    if(length(grep("was not found",urlData))>0){  
      #this disease may not even be in the analyses directory yet.
      stop( paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase,". No data was downloaded.\n"))
    }
    
    #remove any FFPE datasets, or only keep those depending on user inputs.
    if (FFPE) { 
      urlData <- urlData[grep("FFPE",urlData)]    
      if(length(urlData)==0){        
        stop("\nNo FFPE data found for this query. Try FFPE=FALSE.\n")        
      }      
    } else {    
      #we DON'T want FFPE data.
      #but if no FFPE data to begin with: don't subset on this.
      if(length(grep("FFPE",urlData))>0){        
        urlData <- urlData[-grep("FFPE",urlData)]        
      }
      if(length(urlData)==0){        
        stop("\nNo non-FFPE data found for this query. Try FFPE=TRUE.\n")        
      }
    }
    #now get full dataset name.
    fileName <- urlData[grep(dataFileTag,urlData)]
    
    if(length(fileName)==0){      
      warnMessage <- paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase,"     for data type ",dataFileTag ,".No data was downloaded.\n")
      warning(warnMessage)
      return(NA)      
    }
    #some redundancy..but that' OK because we'll add back on the unique tar.gz file tag.
    #first file is one we want - not md5 file.
    fileName <- strsplit2(fileName,"tar.gz")[1,1]
    fileName <- paste(fileName,fileType,sep="")
    
    #final download url
    gdacURL <- paste(gdacURL,fileName,sep="")
    
    # Beed the savedir when we don't download !!!!!!!!
    saveDir <- paste(saveDir,"gdac_",lastDateCompress,'/',sep="")
    
    if(downloadData){        
      cat("\tDownloading",dataFileTag,"data, version:",lastDate,"\n")                
      cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
      
      # create dirs
      dir.create(saveDir,showWarnings=FALSE)
      
      # download file        
      setwd(saveDir)                
      download.file(gdacURL,fileName,quiet=TRUE,mode="wb")
      
      #this assumes a tar.gz file.
      if(fileType=="tar.gz" && untarUngzip) {                    
        cat("\tUnpacking data.\n")
        tarfile=paste0(saveDir,fileName)
        untar(fileName)
        
        #remove tarred file
        fileToRemove <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
        file.remove(paste0(saveDir,fileToRemove))
        
      } else if(untarUngzip) {        
        warning("File expansion/opening only built in for tar.gz files at the moment.\n")        
      }
      
      finalDir <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
      finalDir <- strsplit2(finalDir,fileType)        
      #must remove LAST period (ie laster character) only. 
      finalDir <- substr(finalDir,start=0,stop=(nchar(finalDir)-1))
      finalDir <- paste0(saveDir,finalDir)        
      cat("\tFinished downloading",dataFileTag,"data to",finalDir,"\n")
      
    } else {
      #just spit out the command you need
      cat("download data url is :\n ",gdacURL,'\n')
      finalDir <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
      finalDir <- strsplit2(finalDir,fileType)
      
      #must remove LAST period (ie laster character) only. 
      finalDir <- substr(finalDir,start=0,stop=(nchar(finalDir)-1))
      finalDir <- paste0(saveDir,finalDir)    
    }
    DownloadedFile=paste0(finalDir,'/')
    return(DownloadedFile)
  } else{
    return(cat(paste0("No data correspond to cancer ",TCGA_acronym_uppercase,".\n")))
  }
}

#code to process TCGA GISTIC (CNV) and gene expression (rna-seq or microarray) data 
#downloaded via code in TCGAdownload.R

#top wrapper code that analyzes both data types via one call:
Preprocess_CancerSite <- function(CancerSite,DataSetDirectories) {

  # Processing MA data, special case for OV and GBM where no RNA seq data is available
  if (CancerSite=="OV" || CancerSite=="GBM") { 
    MAstring='transcriptome__agilent'
  } else {
    MAstring='mRNAseq_RSEM_normalized_log2.txt'
  }               
  cat("Loading mRNA data.\n")
  if (length(DataSetDirectories$MAdirectory)>1) {              
    cat("\tFound multiple MA data sets.\n")
    DataSets=list()
    GeneLists=list()
    SampleLists=list()                
    MetaBatchData=data.frame()
    for (i in 1:length(DataSetDirectories$MAdirectory)) {        
      cat("\tProcessing data set",i,"\n")
      MAfiles=dir(DataSetDirectories$MAdirectory[i])
      MatchedFile=grep(MAstring,MAfiles)        
      if (length(MatchedFile)>0) {        
        DataSets[[i]]=Preprocess_MAdata(CancerSite,DataSetDirectories$MAdirectory[i],MAfiles[MatchedFile])
        GeneLists[[i]]=rownames(DataSets[[i]])
        SampleLists[[i]]=colnames(DataSets[[i]])
        
        # growing a batch data object
        currentBatch=matrix(i,length(colnames(DataSets[[i]])),1)                
        currentBatchData=data.frame(ArrayName=colnames(DataSets[[i]]),SampleName=colnames(DataSets[[i]]),Batch=currentBatch)
        MetaBatchData=rbind(MetaBatchData,currentBatchData)
      } else {
        cat("MA file not found for this cancer.\n")
      }           
    }
    # combine data sets with Combat. 
    cat("Combining data sets.\n")
    # overlap genes
    OverlapProbes=Reduce(intersect,GeneLists)
    OverlapSamples=Reduce(intersect,SampleLists)    
    if (length(OverlapSamples)>0) {
      cat('There is overlap between samples. No solution yet.\n')           
    }        
    for (i in 1:length(DataSetDirectories$MAdirectory)) {
      DataSets[[i]]=DataSets[[i]][OverlapProbes,]
    }
    # combat on data sets. 
    MA_TCGA=Reduce(cbind,DataSets)
    MA_TCGA=TCGA_BatchCorrection_MolecularData(MA_TCGA,MetaBatchData,MinPerBatch)    
    
  } else {
    
    MAfiles=dir(DataSetDirectories$MAdirectory)
    MatchedFile=grep(MAstring,MAfiles)        
    if (length(MatchedFile)>0) {   
      
      MA_TCGA=Preprocess_MAdata(CancerSite,DataSetDirectories$MAdirectory,MAfiles[MatchedFile])
      
    } else {
      
      stop("MA file not found for this cancer.\n")
    }           
  }         
  
  # Processing CNV data
  cat("Loading CNV data.\n")
  cat("\tProcessing GISTIC output.\n")
  CGH_Data=TCGA_Load_GISTICdata(DataSetDirectories$CNVdirectory) 
  CNVgenes=c(CGH_Data$AMPgenes,CGH_Data$DELgenes)
  CNVgenes=gsub('\\[','',CNVgenes)
  CNVgenes=gsub('\\]','',CNVgenes)
  
  # selecting the segmented GISTIC data
  if (length(CNVgenes)>1){
    CGH_Data$CGH_Data_Segmented=CGH_Data$CGH_Data_Segmented[unique(CNVgenes),]   
  } else {
    CGH_Data_Segmented_temp = matrix(0,1,ncol(CGH_Data$CGH_Data_Segmented))
    colnames(CGH_Data_Segmented_temp) = colnames(CGH_Data$CGH_Data_Segmented)
    rownames(CGH_Data_Segmented_temp) = CNVgenes
    CGH_Data_Segmented_temp[1,] <- CGH_Data$CGH_Data_Segmented[unique(CNVgenes),]  
    CGH_Data$CGH_Data_Segmented = CGH_Data_Segmented_temp
  }
  SampleNames=colnames(CGH_Data$CGH_Data_Segmented)
  if (CancerSite =='LAML') {
    colnames(CGH_Data$CGH_Data_Segmented)=paste(SampleNames,'-03',sep='')
  } else {
    colnames(CGH_Data$CGH_Data_Segmented)=paste(SampleNames,'-01',sep='')
  }        
  if (length(CNVgenes)>1){
    cat("\tBatch correction.\n")
    CNV_TCGA=TCGA_BatchCorrection_MolecularData(CGH_Data$CGH_Data_Segmented,BatchData,MinPerBatch=5)    
    Genes=rownames(CNV_TCGA)
  } else {
    cat("\tBatch correction can't be done with less than two CNV alterated genes.\n")
    CNV_TCGA = CGH_Data$CGH_Data_Segmented
    Genes = rownames(CNV_TCGA)
  }
  
  SplitGenes=strsplit2(Genes,'\\|')
  rownames(CNV_TCGA)=SplitGenes[,1]      
  CNV_TCGA=TCGA_GENERIC_MergeData(unique(rownames(CNV_TCGA)),CNV_TCGA)        
  
  # Overlapping all data sets (for MET and CNV, making sure all genes and samples exist in the MA data)
  cat('Summarizing:\n')
  cat('\tFound',length(rownames(MA_TCGA)),'genes and',length(colnames(MA_TCGA)),'samples for MA data.\n')    
  cat('\tFound',length(rownames(CNV_TCGA)),'genes and',length(colnames(CNV_TCGA)),'samples for GISTIC data before overlapping.\n')    
  cat('\tFound',length(rownames(MET_TCGA)),'genes and',length(colnames(MET_TCGA)),'samples for MethylMix data before overlapping.\n')    
  
  cat('Overlapping genes and samples for GISTIC and MethylMix.\n')
  OverlapGenes=Reduce(intersect,list(rownames(MA_TCGA),rownames(CNV_TCGA)))
  OverlapSamples=Reduce(intersect,list(colnames(MA_TCGA),colnames(CNV_TCGA)))
  if (length(OverlapGenes)==1 || length(OverlapSamples)==1){
    CNV_TCGA_temp <- matrix(0,length(OverlapGenes),length(OverlapSamples))
    rownames(CNV_TCGA_temp) <- OverlapGenes
    colnames(CNV_TCGA_temp) <- OverlapSamples
    if (length(OverlapGenes)==1){
      CNV_TCGA_temp[1,] <- CNV_TCGA[,OverlapSamples]
    } else {
      CNV_TCGA_temp[,1] <- CNV_TCGA[OverlapGenes,]
    }
    CNV_TCGA = CNV_TCGA_temp
  } else {
    CNV_TCGA=CNV_TCGA[OverlapGenes,OverlapSamples]
  }
  cat('\tFound',length(OverlapGenes),'overlapping genes and',length(OverlapSamples),'samples for GISTIC data.\n')    
  
  OverlapGenes=Reduce(intersect,list(rownames(MA_TCGA),rownames(MET_TCGA)))
  OverlapSamples=Reduce(intersect,list(colnames(MA_TCGA),colnames(MET_TCGA)))
  if (length(OverlapGenes)==1 || length(OverlapSamples)==1){
    MET_TCGA_temp <- matrix(0,length(OverlapGenes),length(OverlapSamples))
    rownames(MET_TCGA_temp) <- OverlapGenes
    colnames(MET_TCGA_temp) <- OverlapSamples
    if (length(OverlapGenes)==1){
      MET_TCGA_temp[1,] <- MET_TCGA[,OverlapSamples]
    } else {
      MET_TCGA_temp[,1] <- MET_TCGA[OverlapGenes,]
    }
    MET_TCGA = MET_TCGA_temp
  } else {
    MET_TCGA=MET_TCGA[OverlapGenes,OverlapSamples]
  }
  cat('\tFound',length(OverlapGenes),'overlapping genes and',length(OverlapSamples),'samples for MethylMix data.\n')
  
  return(list(MA_TCGA=MA_TCGA,CNV_TCGA=CNV_TCGA,MET_TCGA=MET_TCGA))
}

TCGA_Load_GISTICdata <- function (GisticDirectory) {
  
  # loading the data in R
  GenesFile=paste(GisticDirectory,'all_data_by_genes.txt',sep="")
  CGH_Data=read.csv(GenesFile,sep="\t",row.names=1,header=TRUE,na.strings="NA")
  
  # Fix the sample names
  SampleNames=colnames(CGH_Data)
  SampleNames=gsub('\\.','-',SampleNames) # which way do we want to go . or -???
  colnames(CGH_Data)=SampleNames
  
  # Process sample groups 
  Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(CGH_Data))
  
  # CLeaning up sample names
  CGH_Data=TCGA_GENERIC_CleanUpSampleNames(CGH_Data,12)
  
  # remove the first two columns
  CGH_Data=CGH_Data[,-1]
  CGH_Data=CGH_Data[,-1]
  CGH_Data=as.matrix(CGH_Data)
  
  # loading the thresholded data. 
  GenesFile=paste(GisticDirectory,'all_thresholded.by_genes.txt',sep="")
  CGH_Data_Thresholded=read.csv(GenesFile,sep="\t",row.names=1,header=TRUE,na.strings="NA")
  
  # Fix the sample names
  SampleNames=colnames(CGH_Data_Thresholded)
  SampleNames=gsub('\\.','-',SampleNames) # which way do we want to go . or -???
  colnames(CGH_Data_Thresholded)=SampleNames
  
  CGH_Data_Thresholded=TCGA_GENERIC_CleanUpSampleNames(CGH_Data_Thresholded,12)
  CGH_Data_Thresholded=CGH_Data_Thresholded[,-1]
  CGH_Data_Thresholded=CGH_Data_Thresholded[,-1]
  CGH_Data_Thresholded=as.matrix(CGH_Data_Thresholded)
  
  # load the amp and del genes
  AMPfile=paste(GisticDirectory,'amp_genes.conf_99.txt',sep="")
  AMPtable=read.csv(AMPfile,sep="\t",header=TRUE,na.strings="NA")     
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
  
  # Same for DEL genes
  DELfile=paste(GisticDirectory,'del_genes.conf_99.txt',sep="")
  DELtable=read.csv(DELfile,sep="\t",header=TRUE,na.strings="NA")     
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
  
  cat("There are",length(AMPgenes),"AMP genes and",length(DELgenes),"DEL genes.\n")
  
  return(list(CGH_Data_Segmented=CGH_Data,CGH_Data_Thresholded=CGH_Data_Thresholded,AMPgenes=AMPgenes,DELgenes=DELgenes))
}

Preprocess_MAdata <- function(CancerSite,Directory,File) {    
  MinPerBatch=5   
  
  cat("\tMissing value estimation.\n")
  MA_TCGA=TCGA_Load_MolecularData(paste(Directory,File,sep=''))        
  Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MA_TCGA))        
  if (CancerSite =='LAML') {
    MA_TCGA=MA_TCGA[,Samplegroups$PeripheralBloodCancer]
  } else {
    MA_TCGA=MA_TCGA[,Samplegroups$Primary]
  }          
  cat("\tBatch correction.\n")
  MA_TCGA=TCGA_BatchCorrection_MolecularData(MA_TCGA,BatchData,MinPerBatch)
  
  cat("\tProcessing gene ids and merging.\n")
  Genes=rownames(MA_TCGA)
  SplitGenes=strsplit2(Genes,'\\|')
  rownames(MA_TCGA)=SplitGenes[,1]        
  MA_TCGA=MA_TCGA[!rownames(MA_TCGA) %in% '?',]        
  MA_TCGA=TCGA_GENERIC_MergeData(unique(rownames(MA_TCGA)),MA_TCGA)  
  
  return(MA_TCGA=MA_TCGA)
}

TCGA_Load_MolecularData <- function (Filename) {     
  MET_Data=read.csv(Filename,sep="\t",row.names=1,header=TRUE,na.strings=c("NA","null"))
  
  if (rownames(MET_Data)[1]=='Composite Element REF') {
    cat("Removing first row with text stuff.\n")
    MET_Data=MET_Data[-1,]         
    Genes=rownames(MET_Data)
    MET_Data=apply(MET_Data,2,as.numeric)
    rownames(MET_Data)=Genes
  }
  
  SampleNames=colnames(MET_Data)
  SampleNames=gsub('\\.','-',SampleNames) # which way do we want to go . or -???
  colnames(MET_Data)=SampleNames
  
  # more than 10% = removal
  MissingValueThreshold=0.1
  
  # removing clones with too many missing values
  NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x))) /ncol(MET_Data)
  cat("Removing",sum(NrMissingsPerGene>MissingValueThreshold),"genes with more than 10% missing values.\n")
  if (sum(NrMissingsPerGene>MissingValueThreshold)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThreshold,]
  
  # removing patients with too many missings values     
  NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x))) /nrow(MET_Data)
  cat("Removing",sum(NrMissingsPerSample>MissingValueThreshold),"patients with more than 10% missing values.\n")
  if (sum(NrMissingsPerSample>MissingValueThreshold)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThreshold]
  
  # knn impute using Tibshirani's method     
  if (length(colnames(MET_Data))>1) {
    k=15
    KNNresults=impute.knn(as.matrix(MET_Data),k) # this does not work, gives storage error. ???
    MET_Data_KNN=KNNresults$data        
    # cleaning up sample names
    MET_Data_KNN_Clean=TCGA_GENERIC_CleanUpSampleNames(MET_Data_KNN,15)
    return(MET_Data_KNN_Clean)        
  } else {
    # when only 1 sample,need to make a matrix again
    #MET_Data=as.matrix(MET_Data)
    MET_Data_Clean=TCGA_GENERIC_CleanUpSampleNames(MET_Data,15)
    return(MET_Data_Clean)    
  }
  
}

TCGA_GENERIC_CleanUpSampleNames <-function(GEN_Data,IDlength=12) {     
  SampleNames=colnames(GEN_Data)
  SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
  if (length(SampleNamesShort)!=length(unique(SampleNamesShort))) {
    # remove the doubles           
    Counts=table(SampleNamesShort)
    Doubles=rownames(Counts)[which(Counts>1)]
    
    cat("Removing doubles for",length(Doubles),"samples.\n")
    for(i in 1:length(Doubles)) {                         
      CurrentDouble=Doubles[i]          
      pos=grep(CurrentDouble,SampleNames)
      #GEN_Data[1:10,pos]
      #cor(GEN_Data[,pos])
      GEN_Data=GEN_Data[,-pos[2:length(pos)]]     
      SampleNames=colnames(GEN_Data) # need to update samplenames because pos is relative to this
    }
    SampleNames=colnames(GEN_Data)
    SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
    
    # now set the samplenames
    colnames(GEN_Data)=SampleNamesShort
  } else {
    colnames(GEN_Data)=SampleNamesShort     
  }     
  return(GEN_Data)
}

TCGA_GENERIC_GetSampleGroups <-function(SampleNames) {
  
  # First replace any . with - so the sample groups are uniform. 
  SampleGroups=list()
  
  #1: Primary Tumor
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]01[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$Primary=SampleNames[Matches==1]
  
  #2: Recurrent tumor
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]02[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$Recurrent=SampleNames[Matches==1]
  
  #3: Primary blood derived cancer - peripheral blood
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]03[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$PeripheralBloodCancer=SampleNames[Matches==1]     
  
  #10: Blood derived normal
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]10[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$BloodNormal=SampleNames[Matches==1]     
  
  #11: Solid tissue derived normal
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]11[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$SolidNormal=SampleNames[Matches==1]     
  
  #20 Cellines
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]20[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$CellLines=SampleNames[Matches==1]    
  
  #06 Cellines
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]06[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$Metastatic=SampleNames[Matches==1]         
  
  return(SampleGroups)
}

TCGA_BatchCorrection_MolecularData <- function (GEN_Data,BatchData,MinInBatch) {
  
  # Remove samples with batch number 0
  if (length(-which(BatchData[,3]==0))>0) {
    BatchData=BatchData[-which(BatchData[,3]==0),]
  }
  
  # remove batches that are too small     
  #MinInBatch=5
  PresentSamples=is.element(BatchData[,1],colnames(GEN_Data))
  
  # changes this April 2014, such that the BatchDataSelected only deals with samples in the current GEN_Data
  BatchDataSelected=BatchData[PresentSamples,] 
  if (sum(PresentSamples) != length(colnames(GEN_Data))) BatchDataSelected=BatchData[-which(PresentSamples==FALSE),]
  BatchDataSelected$Batch <- factor(BatchDataSelected$Batch)
  
  NrPerBatch=table(BatchDataSelected$Batch)
  SmallBatches=NrPerBatch<MinInBatch
  BatchesToBeRemoved=names(SmallBatches)[which(SmallBatches==TRUE)]
  SamplesToBeRemoved=as.character(BatchDataSelected[which(BatchDataSelected$Batch %in% BatchesToBeRemoved),1])
  
  if (length(colnames(GEN_Data))-length(which(colnames(GEN_Data) %in% SamplesToBeRemoved))>5) { # just checking if we have enough samples after removing the too small batches
    if (length(which(colnames(GEN_Data) %in% SamplesToBeRemoved))>0) {
      cat("Removing",length(which(colnames(GEN_Data) %in% SamplesToBeRemoved)),"samples because their batches are too small.\n")
      GEN_Data=GEN_Data[,-which(colnames(GEN_Data) %in% SamplesToBeRemoved)]
    }          
    # batch correction with Combat, incorporate check for only 1 batch
    BatchCheck=TCGA_GENERIC_CheckBatchEffect(GEN_Data,BatchData)
    
    if (is.list(BatchCheck)) {
      GEN_Data_Corrected=TCGA_GENERIC_BatchCorrection(GEN_Data,BatchData)
      BatchCheck=TCGA_GENERIC_CheckBatchEffect(GEN_Data_Corrected,BatchData)
      return(GEN_Data_Corrected)
    } else {
      cat("Only one batch, no batch correction possible.\n")
      return(GEN_Data)
    }
  } else {
    cat("The nr of samples becomes to small, no batch correction possible.\n")
    return(GEN_Data)
  }
}

TCGA_GENERIC_BatchCorrection <-function(GEN_Data,BatchData) {
  
  # combat
  # select only samples with batch, others get deleted
  WithBatchSamples=is.element(colnames(GEN_Data),BatchData[,1])
  if (length(which(WithBatchSamples==FALSE))>0) GEN_Data=GEN_Data[,-which(WithBatchSamples==FALSE)]
  
  # select only the batch data that is present in the current data set, remove others (remember, the batch data is for all of TCGA)
  PresentSamples=is.element(BatchData[,1],colnames(GEN_Data))
  BatchDataSelected=BatchData
  if (sum(PresentSamples) != length(colnames(GEN_Data))) BatchDataSelected=BatchData[-which(PresentSamples==FALSE),]
  BatchDataSelected$Batch <- factor(BatchDataSelected$Batch)
  BatchDataSelected$ArrayName <- factor(BatchDataSelected$ArrayName)
  
  # reordening samples (not really necessary as Combat does this too)
  order <- match(colnames(GEN_Data),BatchDataSelected[,1])
  BatchDataSelected=BatchDataSelected[order,]
  BatchDataSelected$Batch <- factor(BatchDataSelected$Batch)
  
  # running combat
  CombatResults=ComBat_NoFiles(GEN_Data,BatchDataSelected)
  
  GEN_Data_Corrected=CombatResults[,-1]
  class(GEN_Data_Corrected) <- "numeric"
  return(GEN_Data_Corrected)
}

ComBat_NoFiles <- function(dat, saminfo, type='txt', write=FALSE, covariates='all', par.prior=TRUE, filter=FALSE, skip=0, prior.plots=FALSE) {
  #debug: expression_xls='exp.txt'; sample_info_file='sam.txt'; type='txt'; write=T; covariates='all'; par.prior=T; filter=F; skip=0; prior.plots=T
  cat('Reading Sample Information File\n')
  #saminfo <- read.table(sample_info_file, header=T, sep='\t',comment.char='')
  if(sum(colnames(saminfo)=="Batch")!=1){return('ERROR: Sample Information File does not have a Batch column!')}
  expression_xls='exp.txt'
  cat('Reading Expression Data File\n')
  #      if(type=='csv'){
  #           dat <- read.csv(expression_xls,header=T,row.names=1,as.is=T)
  #           #print(dat[1:2,])
  #           #  dat <- dat[,trim.dat(dat)]  
  #           #print(colnames(dat))
  #           #colnames(dat)=scan(expression_xls,what='character',nlines=1,sep=',',quiet=T)[1:ncol(dat)]
  #           #print(colnames(dat))
  #      }  else {
  #dat <- read.table(expression_xls,header=T,comment.char='',fill=T,sep='\t', as.is=T)
  if (nrow(dat)>1){
    dat <- dat[,trim.dat(dat)]
  } else {
    Trim <- trim.dat(dat)
    dat_temp <- matrix(0,1,ncol=length(Trim))
    colnames(dat_temp) <- colnames(dat)[Trim]
    rownames(dat_temp) <- rownames(dat)
    dat_temp[1,] <- dat[,Trim]
    dat = dat_temp
  }
  #colnames(dat)=scan(expression_xls,what='character',nlines=1,sep='\t',quiet=T)[1:ncol(dat)]
  #      }
  
  if (skip>0){
    geneinfo <- as.matrix(dat[,1:skip])
    dat <- dat[,-c(1:skip)]
  } else {
    geneinfo=NULL
  }
  
  if(filter){
    ngenes <- nrow(dat)
    col <- ncol(dat)/2
    present <- apply(dat, 1, filter.absent, filter)
    dat <- dat[present, -(2*(1:col))]
    if (skip>0){geneinfo <- geneinfo[present,]}
    cat('Filtered genes absent in more than',filter,'of samples. Genes remaining:',nrow(dat),'; Genes filtered:',ngenes-nrow(dat),'\n')
  }
  
  if(any(apply(dat,2,mode)!='numeric')){return('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option)')}
  
  tmp <- match(colnames(dat),saminfo[,1])
  if(any(is.na(tmp))){return('ERROR: Sample Information File and Data Array Names are not the same!')}
  tmp1 <- match(saminfo[,1],colnames(dat))
  saminfo <- saminfo[tmp1[!is.na(tmp1)],]      
  
  if(any(covariates != 'all')){saminfo <- saminfo[,c(1:2,covariates)]}
  design <- design.mat(saminfo)    
  
  batches <- list.batch(saminfo)
  n.batch <- length(batches)
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  
  ## Check for missing values
  NAs = any(is.na(dat))
  if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
  #print(dat[1:2,])
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))}else{B.hat=apply(dat,1,Beta.NA,design)} #Standarization Model
  grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
  if (!NAs){var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=TRUE)}
  
  stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
  if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}    
  s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
  
  ##Get regression batch effect parameters
  cat("Fitting L/S model and finding priors\n")
  batch.design <- design[,1:n.batch]
  if (!NAs){gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))}else{gamma.hat=apply(s.data,1,Beta.NA,batch.design)}
  delta.hat <- NULL
  for (i in batches){
    delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=TRUE))
  }
  
  ##Find Priors
  gamma.bar <- apply(gamma.hat, 1, mean)
  t2 <- apply(gamma.hat, 1, var)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  
  
  ##Plot empirical and parametric priors    
  if (prior.plots & par.prior){
    pdf(file='prior_plots.pdf')
    par(mfrow=c(2,2))
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main="Density Plot")
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
    qqnorm(gamma.hat[1,])    
    qqline(gamma.hat[1,], col=2)    
    
    tmp <- density(delta.hat[1,])
    invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
    tmp1 <- density(invgam)
    plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
    lines(tmp1, col=2)
    qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')    
    lines(c(0,max(invgam)),c(0,max(invgam)),col=2)    
    title('Q-Q Plot')
    dev.off()
  }
  
  ##Find EB batch adjustments    
  gamma.star <- delta.star <- NULL
  if(par.prior){
    cat("Finding parametric adjustments\n")
    for (i in 1:n.batch){
      temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }else{
    cat("Finding nonparametric adjustments\n")
    for (i in 1:n.batch){
      temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }
  
  ### Normalize the Data ###
  cat("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
  if(write) {
    output_file <- paste(expression_xls,'Adjusted','.txt',sep='_')
    #print(geneinfo[1:2])
    #print(bayesdata[1:2,1:4])
    #cat(c(colnames(geneinfo),colnames(dat),'\n'),file=output_file,sep='\t')
    #suppressWarnings(write.table(cbind(geneinfo,formatC(as.matrix(bayesdata), format = "f")), file=output_file, sep="\t", quote=F,row.names=F,col.names=F,append=T))
    outdata <- cbind(ProbeID=rownames(dat), bayesdata)
    write.table(outdata, file=output_file, sep="\t")
    cat("Adjusted data saved in file:",output_file,"\n")
  } else {
    return(cbind(rownames(dat),bayesdata))
  }  
}

# filters data based on presence/absence call
filter.absent <- function(x,pct){
  present <- TRUE
  col <- length(x)/2
  pct.absent <- (sum(x[2*(1:col)]=="A") + sum(x[2*(1:col)]=="M"))/col
  if(pct.absent > pct){present <- FALSE}
  present
}

# Next two functions make the design matrix (X) from the sample info file 
build.design <- function(vec, des=NULL, start=2){
  tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
  for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
  cbind(des,tmp)
}

design.mat <- function(saminfo){
  tmp <- which(colnames(saminfo) == 'Batch')
  tmp1 <- as.factor(saminfo[,tmp])
  cat("Found",nlevels(tmp1),'batches\n')
  design <- build.design(tmp1,start=1)
  ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
  cat("Found",ncov,'covariate(s)\n')
  if(ncov>0){
    for (j in 1:ncov){
      tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
      design <- build.design(tmp1,des=design)
    }
  }
  design
}

# Makes a list with elements pointing to which array belongs to which batch
list.batch <- function(saminfo){
  tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
  batches <- NULL
  for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
  batches
}

# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
  tmp <- strsplit(colnames(dat),'\\.')
  tr <- NULL
  for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
  tr
}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}

# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments
it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- apply(!is.na(sdat),1,sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat,g.bar,n,d.old,t2)
    sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=TRUE)
    d.new <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  #cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

#likelihood function used below
L <- function(x,g.hat,d.hat){prod(dnorm(x,g.hat,sqrt(d.hat)))}

# Monte Carlo integration function to find the nonparametric adjustments
int.eprior <- function(sdat,g.hat,d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]        
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x),length(g),n,byrow=TRUE)
    resid2 <- (dat-g)^2
    sum2 <- resid2%*%j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star,sum(g*LH)/sum(LH))
    d.star <- c(d.star,sum(d*LH)/sum(LH))
    #if(i%%1000==0){cat(i,'\n')}
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  adjust    
} 

#fits the L/S model in the presence of missing data values

Beta.NA = function(y,X){
  des=X[!is.na(y),]
  y1=y[!is.na(y)]
  B <- solve(t(des)%*%des)%*%t(des)%*%y1
  B
}

TCGA_GENERIC_MergeData <-function(NewIDListUnique,DataMatrix,MergeMethod) {
  
  NrUniqueGenes=length(NewIDListUnique)
  MergedData=matrix(0,NrUniqueGenes,length(colnames(DataMatrix)))
  for (i in 1:NrUniqueGenes) {
    currentID=NewIDListUnique[i]          
    tmpData=DataMatrix[which(rownames(DataMatrix) %in% currentID),]
    if (length(rownames(tmpData)) >1) {
      MergedData[i,]=colMeans(tmpData)
    } else {
      MergedData[i,]=tmpData
    }
  }
  rownames(MergedData)=NewIDListUnique
  colnames(MergedData)=colnames(DataMatrix)
  
  return(MergedData)
}

TCGA_GENERIC_CheckBatchEffect <-function(GEN_Data,BatchData) {
  # GEN_Data ;;;
  # Barch
  
  # first match the samples to the batch
  Order=match(colnames(GEN_Data),BatchData[,1])
  BatchDataSelected=BatchData[Order,]
  BatchDataSelected$Batch <- factor(BatchDataSelected$Batch) 
  
  # PCA analysis
  # alternatively use fast.prcomp from package gmodels, but tests do not show this is faster
  PCAanalysis=prcomp(t(GEN_Data))
  PCdata=PCAanalysis$x
  plot(PCdata[,1]~BatchDataSelected[,3])
  
  if (length(unique(BatchDataSelected$Batch))>1) {
    tmp=aov(PCdata[,1]~BatchDataSelected[,3])          
    return(list(Pvalues=summary(tmp),PCA=PCdata,BatchDataSelected=BatchDataSelected))
  } else {
    return(-1)
  }
}

Save_CancerSite <- function(CancerSite,TargetDirectory,DataSetDirectories,ProcessedData) {
  
  if (length(DataSetDirectories$MAdirectory)>1) {
    tmpMAdir=DataSetDirectories$MAdirectory[1]
  } else {
    tmpMAdir=DataSetDirectories$MAdirectory
  }
  
  MAdir=strsplit2(tmpMAdir,'/')
  tmpPos=grep('gdac',MAdir)
  MAversion=gsub('gdac_','',MAdir[tmpPos[1]])
  
  CNVdir=strsplit2(DataSetDirectories$CNVdirectory,'/')
  CNVversion=gsub('gdac_','',CNVdir[tmpPos[1]])
  
  METversion="MethylMix2015"
  
  SaveFile=paste(TargetDirectory,'TCGA_',CancerSite,"_ProcessedData_MA",MAversion,"_CNV",CNVversion,"_MET",METversion,".RData",sep='')
  save(file=SaveFile,ProcessedData)
  return(SaveFile)
}