
#####Ovarian Cancer signature
#source("https://bioconductor.org/biocLite.R")
#biocLite("MetaGxOvarian")
library(MetaGxOvarian)
library(consensusOV)

source("./Functions/functions.R")
#with OS status and grade according https://www.biorxiv.org/content/biorxiv/early/2016/05/12/052910.full.pdf
datasets= c("PMID17290060","GSE9891","GSE8842","GSE51088","GSE49997","GSE32063","GSE32062","GSE30161","GSE30009","GSE26193","GSE18520","GSE17260","GSE14764","GSE13876","GSE26712","TCGAOVARIAN")
platform= c("GPL96", "GPL570","GPL5689","GPL7264","GPL2986","GPL6480","GPL6480","GPL570","GPL13728","GPL570","GPL570","GPL6480","GPL96","GPL7759","GPL96","GPL3921")
affyplatforms= c("GPL96","GPL97","GPL570","GPL571","GPL4685","GPL6244") 
downl= datasets[which(platform %in% affyplatforms)]

ALL= loadOvarianEsets(removeDuplicates = TRUE,
                     quantileCutoff = 0, rescale = FALSE, minNumberGenes = 0,
                     minNumberEvents = 0, minSampleSize = 0, removeRetracted = TRUE,
                     removeSubsets = TRUE, keepCommonOnly = FALSE, imputeMissing = FALSE)

esets=ALL$esets[downl]


for(i in downl){

RFS_time= as.numeric(pData(esets[[i]])$days_to_tumor_recurrence)
RFS_status= as.character(pData(esets[[i]])$recurrence_status)
RFS_status[RFS_status=="recurrence"]=1
RFS_status[RFS_status=="norecurrence"]=0
RFS_status=as.numeric(RFS_status)

OS_time= as.numeric(pData(esets[[i]])$days_to_death)
OS_status= as.character(pData(esets[[i]])$vital_status)
OS_status[OS_status=="deceased"]=1
OS_status[OS_status=="living"]=0
OS_status=as.numeric(OS_status)

esetClinic= data.frame(histology=pData(esets[[i]])$histological_type,grade=pData(esets[[i]])$summarygrade, age=pData(esets[[i]])$age_at_initial_pathologic_diagnosis,days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=pData(esets[[i]])$tumorstage,batch=i,row.names=rownames(pData(esets[[i]])))

pData(esets[[i]])=esetClinic

ct=censor.time(pData(esets[[i]])$days_to_os, pData(esets[[i]])$os_status, time.cens = 1825 )
pData(esets[[i]])$time=ct$surv.time.cens
pData(esets[[i]])$status=ct$surv.event.cens


esets[[i]]= esets[[i]][,pData(esets[[i]])$histology=="ser" & !is.na(pData(esets[[i]])$histology)]
esets[[i]]= esets[[i]][,pData(esets[[i]])$grade=="high" & !is.na(pData(esets[[i]])$grade)]
esets[[i]]= esets[[i]][,!is.na(Surv(pData(esets[[i]])$time,pData(esets[[i]])$status))]


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
colnames(fData(esets[[i]]))[1]="ID"
fData(esets[[i]])=fData(esets[[i]])[,c("ID","EntrezGene.ID","gene")]


 #Apply robust linear scaling
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3283537/#bib61
Expression= apply(exprs(esets[[i]]),1,genefu::rescale,na.rm=TRUE,q=0.05)
Expression=t(Expression)
exprs(esets[[i]])=Expression


esets[[i]]$Verhaak.subtypes=get.verhaak.subtypes (exprs(esets[[i]]),fData(esets[[i]])$EntrezGene.ID)$Verhaak.subtypes
esets[[i]]$Bentink.subtypes=get.bentink.subtypes (exprs(esets[[i]]),fData(esets[[i]])$EntrezGene.ID)$Bentink.subtypes
esets[[i]]$Helland.subtypes=get.helland.subtypes (exprs(esets[[i]]),fData(esets[[i]])$EntrezGene.ID)$Helland.subtypes
esets[[i]]$Konecny.subtypes=get.konecny.subtypes (exprs(esets[[i]]),fData(esets[[i]])$EntrezGene.ID)$Konecny.subtypes
esets[[i]]$Consensus.subtypes=get.consensus.subtypes (exprs(esets[[i]]),fData(esets[[i]])$EntrezGene.ID)$consensusOV.subtypes

}

library(TCGAbiolinks)
setwd("/home/mguerrero/Genetic_alg/App_FINAL/Data")


query <- GDCquery(project = "TCGA-OV",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification", 
                 sample.type=c("Primary solid Tumor"),
                 workflow.type = "HTSeq - Counts",
                 legacy=FALSE)


GDCdownload(query)

data <- GDCprepare(query, save = TRUE, 
                   save.filename = "RNA_OV.rda",
                   remove.files.prepared = F)

#Open expression matrix
load("./Data/RNA_OV.rda")
library(SummarizedExperiment)

data=data[,!duplicated(colData(data)$patient)]

#Voom normalization
library(limma)
library(edgeR)
Counts=assays(data)$HTSeq

#Filter low count genes
keep= filterByExpr(Counts)
data=data[keep,]
filteredCounts=Counts[keep,]
VoomCounts=voom(filteredCounts,plot=TRUE)
assays(data)$HTSeq=VoomCounts$E
colnames(VoomCounts$E)=colData(data)$patient


histology= "ser"
grade="high"

RFS_time= NA
RFS_status= NA


OS_status=colData(data)$"vital_status"
OS_status[colData(data)$"vital_status"=="dead"]= 1
OS_status[colData(data)$"vital_status"=="alive"]= 0
OS_status=as.numeric(OS_status)

OS_time= as.numeric(colData(data)$"days_to_death")
OS_time[is.na(OS_time)]=as.numeric(colData(data)$"days_to_last_follow_up")[is.na(OS_time)]

stageall=NA


esetClinic= data.frame(histology=histology,grade=grade, age=colData(data)$age_at_diagnosis/365,days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=stageall,batch="TCGA.RNAseq",row.names=rownames(colData(data)))


ct=censor.time(esetClinic$days_to_os, esetClinic$os_status, time.cens = 1825 )
esetClinic$time=ct$surv.time.cens
esetClinic$status=ct$surv.event.cens
rownames(esetClinic)=colData(data)$patient
phenoData <- new("AnnotatedDataFrame",data=esetClinic)

#Reformulate fData
library(hgu133plus2.db)

ID     <- rowData(data)$ensembl_gene_id
Symbol <- rowData(data)$external_gene_name
Entrez <- mapIds(hgu133plus2.db, keys=Symbol, c("ENTREZID"), keytype="SYMBOL")


featureData= new("AnnotatedDataFrame",data=data.frame(ID=ID,EntrezGene.ID=Entrez,gene=Symbol,row.names=ID))
eset7= ExpressionSet(assayData=VoomCounts$E,phenoData=phenoData,featureData=featureData)

eset7= eset7[,!is.na(Surv(eset7$time,eset7$status))]


#Drop NA genes
Out=which(is.na(fData(eset7)$EntrezGene.ID))
eset7=eset7[-Out,]

#Drop duplicated genes (keep genes with highest variance)
Dup=unique(fData(eset7)[which(duplicated(fData(eset7)$EntrezGene.ID)),"EntrezGene.ID"])
Var= apply(exprs(eset7),1,var)
drop=NULL
for(j in Dup){
  pos=which(fData(eset7)$EntrezGene.ID==j)
  drop= c(drop,pos[-which.max(Var[pos])])
}
eset7=eset7[-drop,]
featureNames(eset7) <- fData(eset7)$EntrezGene.ID
eset7=expandProbesets(eset7)
fData(eset7)=fData(eset7)[,c("ID","EntrezGene.ID","gene")]

 #Apply robust linear scaling
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3283537/#bib61
Expression= apply(eset7,1,genefu::rescale,na.rm=TRUE,q=0.05)
Expression=t(Expression)
exprs(eset7)=Expression


eset7$Verhaak.subtypes=get.verhaak.subtypes (exprs(eset7),fData(eset7)$EntrezGene.ID)$Verhaak.subtypes
eset7$Bentink.subtypes=get.bentink.subtypes (exprs(eset7),fData(eset7)$EntrezGene.ID)$Bentink.subtypes
eset7$Helland.subtypes=get.helland.subtypes (exprs(eset7),fData(eset7)$EntrezGene.ID)$Helland.subtypes
eset7$Konecny.subtypes=get.konecny.subtypes (exprs(eset7),fData(eset7)$EntrezGene.ID)$Konecny.subtypes
eset7$Consensus.subtypes=get.consensus.subtypes (exprs(eset7),fData(eset7)$EntrezGene.ID)$consensusOV.subtypes


downl=c(downl,"TCGA.RNAseq")
esets= append(esets,eset7)
names(esets)=downl

inALL=Reduce(intersect, lapply(esets,featureNames))
inALLpData=Reduce(intersect, lapply(esets,function(x) colnames(pData(x))))
for(i in downl){
  esets[[i]]=esets[[i]][inALL,]
  pData(esets[[i]])=pData(esets[[i]])[,inALLpData]
}


save(esets,file="/home/mguerrero/Genetic_alg/App_FINAL/Data/RNA_OV_all.rda")
