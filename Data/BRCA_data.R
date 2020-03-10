#Files
#Metabric & OsloVal data tidy
#Data source:
#https://www.synapse.org/#!Synapse:syn1688369/wiki/27311
#https://www.synapse.org/#!Synapse:syn1688370/wiki/27309

#Synapse client
#As a sudo (in cmd)
#install pip
#pip install synapseclient
#synapse login -u me@example.com -p password --rememberMe
#synapse get syn1757053 #clinical data (rbin)
#synapse get syn1757063 #expression (rbin)
#synapse get syn1757057 #survival (rbin)
#synapse get syn1757055 #Disease free surv (rbin)

source("./Functions/functions.R")

metabric_files="./Data/Metabric/" #Directory where the files were downloaded" #Directory where the files were downloaded

load(paste0(metabric_files, "Complete_METABRIC_Clinical_Features_Data.rbin"))
load(paste0(metabric_files, "Complete_METABRIC_Expression_Data.rbin"))
load(paste0(metabric_files, "Complete_METABRIC_Clinical_Survival_Data_OS.rbin"))
load(paste0(metabric_files, "Complete_METABRIC_Clinical_Survival_Data__DSS.rbin"))


#Download of the Illumina probes metadata
download.file("https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/v3/HumanHt-12_V3_0_R3_11283641_A.zip","./Data/HumanHt-12_V3_0_R3_11283641_A.zip")
unzip("./Data/HumanHt-12_V3_0_R3_11283641_A.zip",exdir="./Data")
ANNOT=readBGX("./Data/HumanHt-12_V3_0_R3_11283641_A.bgx")
Expression=exprs(Complete_METABRIC_Expression_Data)


#############################
#Expressionset tidy function#
#############################

#Tidy METABRIC and OSLO data

expressiontidy2= function(Expression,ANNOT, CLINIC, SURV,BATCH){
  require(illuminaio)
  require(AnnotationDbi)
  require(limma)
  
  CLINIC$ER.Expr= gsub("\\+", "positive", CLINIC$ER.Expr)
  CLINIC$ER.Expr= gsub("\\-", "negative", CLINIC$ER.Expr)
  CLINIC$PR.Expr= gsub("\\+", "positive", CLINIC$PR.Expr)
  CLINIC$PR.Expr= gsub("\\-", "negative", CLINIC$PR.Expr)
  CLINIC$Her2.Expr= gsub("\\+", "positive", CLINIC$Her2.Expr)
  CLINIC$Her2.Expr= gsub("\\-", "negative", CLINIC$Her2.Expr)
  SURV[,2]= gsub("1", "deceased", SURV[,2])
  SURV[,2]= gsub("0", "living", SURV[,2])
  pData= data.frame(
    row.names= make.names(rownames(CLINIC)),
    sample_name= make.names(rownames(CLINIC)),
    alt_sample_name= NA,
    unique_patient_ID= NA,
    sample_type= "tumor",
    er= CLINIC$ER.Expr,
    pgr= CLINIC$PR.Expr,
    her2= CLINIC$Her2.Expr,
    tumor_size= CLINIC$size/10,
    T=NA,
    N=NA,
    age_at_initial_pathologic_diagnosis= CLINIC$age_at_diagnosis,
    grade=CLINIC$grade,
    dmfs_days=NA,
    dmfs_status=NA,
    days_to_tumor_recurrence=NA,
    recurrence_status=NA,
    days_to_death=SURV[,1],
    vital_status=SURV[,2],
    tissue=NA,
    treatment=CLINIC$Treatment,
    percent_normal_cells=NA,
    percent_stromal_cells=NA,
    percent_tumor_cells=NA,
    batch=BATCH,
    uncurated_author_metadata=NA,
    duplicates=NA)
  probes=rownames(Expression)
  ANNOT=ANNOT$probe[ANNOT$probe[,"Probe_Id"] %in% probes,]
  ANNOT= ANNOT[!is.na(ANNOT$Entrez_Gene_ID),]
  Expression=Expression[ANNOT[,"Probe_Id"],]
  fData= data.frame(probeset=ANNOT$Probe_Id,gene=ANNOT$Symbol, EntrezGene.ID=ANNOT$Entrez_Gene_ID,best_probe=TRUE,row.names=ANNOT$Probe_Id)
  colnames(Expression)=make.names(colnames(Expression))
  phenoData <- new("AnnotatedDataFrame",data=pData)
  featureData <- new("AnnotatedDataFrame",data=fData)
  ExpSet= ExpressionSet(assayData=Expression,phenoData=phenoData,featureData=featureData)
  return(ExpSet)
}





Expset=expressiontidy2(Expression,ANNOT,Complete_METABRIC_Clinical_Features_Data,Complete_METABRIC_Clinical_Survival_Data__DSS,BATCH="metabric") 

metabric=Expset

rm(list=c("Complete_METABRIC_Expression_Data","Complete_METABRIC_Clinical_Features_Data","Complete_METABRIC_Clinical_Survival_Data__DSS","Complete_METABRIC_Clinical_Survival_Data_OS","ANNOT","Expression","metabric_files"))

#Oslo files
#Synapse client
#As a sudo (in cmd)
#install pip
#pip install synapseclient
#synapse login -u me@example.com -p secret --rememberMe
#synapse get syn2188646 / syn1710393 #clinical data (txt/rbin)
#synapse get syn2188662 / syn1710395 #expression (txt/rbin)
#synapse get syn2188659 /syn1710396 #survival (txt/rbin)
#unzip archive.zip

oslo_files="./Data/Oslo/" # Directory where the files were downloaded

load(paste0(oslo_files, "OsloValClinicalData.rbin"))
load(paste0(oslo_files, "OsloValExpn.rbin"))
load(paste0(oslo_files, "OsloValSurvival.rbin"))


#Download of the Illumina probes metadata
download.file("https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/HumanHT-12_V4_0_R2_15002873_B.zip","./Data/HumanHT-12_V4_0_R2_15002873_B.zip")
unzip("./Data/HumanHT-12_V4_0_R2_15002873_B.zip",exdir="./Data") 
ANNOT=readBGX("./Data/HumanHT-12_V4_0_R2_15002873_B.bgx")
Expression=exprs(OsloValExpn)
Expset=expressiontidy2(Expression,ANNOT,OsloValClinicalData,OsloValSurvival,BATCH="oslo") 

oslo=Expset

rm(list=c("OsloValExpn","OsloValClinicalData","OsloValSurvival","ANNOT","Expression","oslo_files"))

#MetaGxBreast
library(MetaGxBreast)

#with Recurrence status according https://www.biorxiv.org/content/biorxiv/early/2016/05/12/052910.full.pdf
#"GSE58644" with error, partially solved with clinical data from geo
#CAL dataset with missmatch labels
datasets= c("UPP","UNT","UNC4","UCSF","TRANSBIG","STNO2","STK","PNC","NKI","NCI","CAL","MAINZ","EMC2","DFHCC","VDX","GSE25066","GSE58644","MSK")
platform= c("GPL97","GPL96","GPL885","NA","GPL96","NA","GPL96","GPL570","Agilent","custom","GPL4685","GPL96","GPL570","GPL570","GPL96","GPL96","GPL6244","GPL96")
affyplatforms= c("GPL96","GPL97","GPL570","GPL571","GPL4685","GPL6244")

downl= datasets[which(platform %in% affyplatforms)]
downl= downl[-which(downl=="UNT")] #UNT data is scaled and is not possible to get AIMS classification
downl= downl[-which(downl=="CAL")] #CAL data does not perform well in any signature, unknown bias

ALL= loadBreastEsets(loadString = downl, removeDuplicates = TRUE,
                     quantileCutoff = 0, rescale = FALSE, minNumberGenes = 0,
                     minNumberEvents = 0, minSampleSize = 0, removeRetracted = TRUE,
                     removeSubsets = TRUE, keepCommonOnly = FALSE, imputeMissing = FALSE)

esets=ALL$esets

#Correcting GSE58644 data
library(GEOquery)
gds <- getGEO("GSE58644")
gds <- gds[[1]]
original_time= as.numeric(pData(gds)$"time:ch1")*30.41 #(original data is in months and MetaGxBreast values are in days)

original_status= as.numeric(pData(gds)$"event:ch1")
original_status[original_status==1 & !is.na(original_status)]='recurrence'
original_status[original_status==0 & !is.na(original_status)]='norecurrence'


esets2= esets[["GSE58644"]]
MetaGx_time= pData(esets2)$dmfs_days
MetaGx_status= pData(esets2)$dmfs_status

identical(rownames(pData(esets2)),rownames(pData(gds))) #TRUE, Patients have same name and order 

pData(esets[["GSE58644"]])$dmfs_days=original_time
pData(esets[["GSE58644"]])$dmfs_status=original_status
fData(esets[["GSE58644"]])$gene=gsub(" ","",fData(esets[["GSE58644"]])$gene)


for(i in downl){
pData(esets[[i]])$status<- NA
pData(esets[[i]])[is.na(pData(esets[[i]])$recurrence_status),"status"]= pData(esets[[i]])[is.na(pData(esets[[i]])$recurrence_status),"dmfs_status"]
pData(esets[[i]])[!is.na(pData(esets[[i]])$recurrence_status),"status"]= pData(esets[[i]])[!is.na(pData(esets[[i]])$recurrence_status),"recurrence_status"]
pData(esets[[i]])$status= as.numeric(pData(esets[[i]])$status=="recurrence")
pData(esets[[i]])$time<- NA
pData(esets[[i]])[is.na(pData(esets[[i]])$recurrence_status),"time"]= pData(esets[[i]])[is.na(pData(esets[[i]])$recurrence_status),"dmfs_days"]
pData(esets[[i]])[!is.na(pData(esets[[i]])$recurrence_status),"time"]= pData(esets[[i]])[!is.na(pData(esets[[i]])$recurrence_status),"days_to_tumor_recurrence"]

ct=censor.time(pData(esets[[i]])$time, pData(esets[[i]])$status, time.cens = 5475 )
pData(esets[[i]])$time=ct$surv.time.cens
pData(esets[[i]])$status=ct$surv.event.cens

}

pData(metabric)$status= as.numeric(pData(metabric)$vital_status=="deceased")
pData(metabric)$time= as.numeric(as.character(pData(metabric)$days_to_death))
ct=censor.time(pData(metabric)$time, pData(metabric)$status, time.cens = 5475 )
pData(metabric)$time=ct$surv.time.cens
pData(metabric)$status=ct$surv.event.cens
esets$metabric=metabric


pData(oslo)$status= as.numeric(pData(oslo)$vital_status=="deceased")
pData(oslo)$time= as.numeric(as.character(pData(oslo)$days_to_death))
ct=censor.time(pData(oslo)$time, pData(oslo)$status, time.cens = 5475 )
pData(oslo)$time=ct$surv.time.cens
pData(oslo)$status=ct$surv.event.cens
esets$oslo=oslo

downl=c(downl,"metabric","oslo")

#################
#expandProbesets#
#################
#Expand probesets for multiple genes for a single probe

expandProbesets <- function (eset, sep = "///")
{
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  y<- lapply(as.character(fData(eset)$gene), function(x) strsplit(x, sep))
  eset <- eset[order(sapply(x, length)), ]
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  y<- lapply(as.character(fData(eset)$gene), function(x) strsplit(x, sep))
  idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
  idy <- unlist(sapply(1:length(y), function(i) rep(i, length(y[[i]]))))
  xx <- !duplicated(unlist(x))
  idx <- idx[xx]
  idy <- idy[xx]
  x <- unlist(x)[xx]
  y <- unlist(y)[xx]
  
  eset <- eset[idx, ]
  featureNames(eset) <- x
  fData(eset)$EntrezGene.ID <- x
  fData(eset)$gene <- y
  return(eset)
  
}



for(i in downl){

  #Drop duplicated genes (keep genes with highest variance)
  Dup=unique(fData(esets[[i]])[which(duplicated(fData(esets[[i]])$EntrezGene.ID)),"EntrezGene.ID"])
  Var= apply(exprs(esets[[i]]),1,var)
  drop=NULL
  for(j in Dup){
   pos=which(fData(esets[[i]])$EntrezGene.ID==j)
    drop= c(drop,pos[-which.max(Var[pos])])
    
  }
  esets[[i]]=esets[[i]][-drop,]

featureNames(esets[[i]]) <- fData(esets[[i]])$EntrezGene.ID

esets[[i]]=expandProbesets(esets[[i]])
 
#Apply robust linear scaling
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3283537/#bib61
Expression= apply(exprs(esets[[i]]),1,genefu::rescale,na.rm=TRUE,q=0.05)

data=Expression
annot= fData(esets[[i]])
colnames(annot)[2]="Gene.Symbol"
annot$probe=annot$EntrezGene.ID

#Perform molecular subtyping
AIMS<-  molecular.subtyping(sbt.model = "AIMS",data = t(exprs(esets[[i]])), annot = annot,do.mapping = TRUE)
pData(esets[[i]])$AIMS=AIMS$subtype
PAM50<- molecular.subtyping(sbt.model = "pam50",data = data, annot = annot,do.mapping = TRUE)
pData(esets[[i]])$pam50=PAM50$subtype
SCM<- molecular.subtyping(sbt.model = "scmgene",data = data,annot = annot,do.mapping = TRUE)
pData(esets[[i]])$scm=SCM$subtype
IntClust<- molecular.subtyping(sbt.model = "intClust",data = data,annot = annot,do.mapping = TRUE)
pData(esets[[i]])$intclust=IntClust$subtype
ONCOTYPE= oncotypedx(data=data,annot=annot,do.mapping=TRUE)
pData(esets[[i]])$oncotypedx=ONCOTYPE$score
#GGI= ggi()
ENDO=endoPredict(data=data,annot=annot,do.mapping=TRUE)
pData(esets[[i]])$endopredict=ENDO$score
MAMMA=gene70(data=data,annot=annot,do.mapping=TRUE)
pData(esets[[i]])$Mammaprint=MAMMA$score

Expression=t(Expression)
exprs(esets[[i]])=Expression

}


inALL=Reduce(intersect, lapply(esets,featureNames))
inALLpData=Reduce(intersect, lapply(esets,function(x) colnames(pData(x))))
for(i in downl){
  esets[[i]]=esets[[i]][inALL,]
  pData(esets[[i]])=pData(esets[[i]])[,inALLpData]
}
for(i in downl){
esets[[i]]=esets[[i]][,!is.na(Surv(pData(esets[[i]])$time,pData(esets[[i]])$status))]
}

save(esets,file="./Data/RNA_BRCA.rda") #Save in the target directory
