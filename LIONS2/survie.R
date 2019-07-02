### Analyse de survie
path <- "/Users/mchampion/Desktop/LIONS/"
DataDirectory <- paste0(path,"Data/TCGA/")

library(stringr)
library(survival)

# load clinical data
SurvData <- read.csv(paste0(DataDirectory,"survie.csv"),sep=";")
SurvData$patient.days_to_death <- as.numeric(SurvData$patient.days_to_death)
SurvData$patient.days_to_last_followup <- as.numeric(SurvData$patient.days_to_last_followup)
SurvData$patient.vital_status <- as.character(SurvData$patient.vital_status)
rownames(SurvData) <- SurvData$patient.bcr_patient_barcode
rownames(SurvData) <- paste0(rownames(SurvData),"-01")
rownames(SurvData) <- str_to_upper(rownames(SurvData))

# creation de l'objet survie
Survie <- data.frame(days_to_event=pmax(SurvData$patient.days_to_death,SurvData$patient.days_to_last_followup,na.rm=T),event=SurvData$patient.vital_status)
rownames(Survie) <- rownames(SurvData)
Survie <- Survie[-which(Survie$days_to_event<=0),]
Survie$event <- as.character(Survie$event)
Survie$event[Survie$event=="alive"] <- rep(0,length(which(Survie$event=="alive"))) 
Survie$event[Survie$event=="dead"] <- rep(1,length(which(Survie$event=="dead")))
Survie$event <- as.numeric(Survie$event)

mydata.surv = Surv(time=Survie$days_to_event,event = Survie$event)

# survival estimation
surv <- survfit(Surv(days_to_event,event)~1,data=Survie)
plot(surv)

# survival analysis by data
load(paste0(DataDirectory,"/ProcessedData_BLCA.Rdata"))
MA <- t(ProcessedData$MA_TCGA)
CNV <- t(ProcessedData$CNV_TCGA)
CNV_thr <- t(ProcessedData$CNV_TCGA_thresholded)
rownames(CNV_thr) <- paste0(rownames(CNV_thr),"-01")

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

all_subtypes <-c("Basal_squamous","Luminal","Luminal_infiltrated","Luminal_papillary","Neuronal")
which.subtype <- all_subtypes[5]

MA_sub <- MA[which(Subtypes==which.subtype),]
Survie_sub <- Survie[which(Subtypes==which.subtype),]
# genes_basal <- c("SPOCD1","ZNF382","RCOR2","ATM","HABP4","IRX3","IFI16","TEAD2","NOTCH4","ZNF211","SNAI2","LMO3","MAFG","CRY1")
# genes_lum <- c("ZNF268","HES2","TBX2","PRDM8","TSHZ1")
# genes_lumInfil <- c("TSHZ1","ZNF354B","AR","HES2","HTATIP2","MAFG","ENO1")
# genes_lumPap <- c("RARB","RFX5","CBFA2T3","TBX18","TBX3","PTRF","TBX2","PPARG","NCOR2")
# genes_neuronal <- c("FAIM3","SMARCA2","RARB","ZNF235")
Exp_basal <- MA_sub[,genes_neuronal]
Survie_sub <- cbind(Survie_sub,Exp_basal)

for (i in (3:(ncol(Exp_basal)+2))){
  Survie_sub[Survie_sub[,i]<median(Survie_sub[,i]),i] <- 'LowExpressed'
  Survie_sub[!Survie_sub[,i]=="LowExpressed",i] <- 'HighlyExpressed'
}

for (i in (3:(ncol(Exp_basal)+2))){
  surv <- survfit(Surv(days_to_event,event)~Survie_sub[,i],data=Survie_sub)
  plot(surv,main=paste0("Survival analysis for gene ",colnames(Survie_sub)[i]),col=2:3,lty=2)
  legend(100,0.3,c("HighlyExpressed","LowExpressed"),col=2:3,lty = 2)
  
#  fit <- coxph(Surv(days_to_event, event) ~ c,data=Survie_sub) 
#  plot(survfit(fit, newdata=data.frame(Survie_sub[,i])),
#       xscale=365.25, xlab = "Years", ylab="Survival") 
  test <- survdiff(Surv(days_to_event,event)~Survie_sub[,i],data=Survie_sub)
  p.val <- 1 - pchisq(test$chisq, length(test$n) - 1)
  print(paste0("P-valeur du test vaut : ",p.val," pour le gene ",colnames(Survie_sub)[i]))
}
