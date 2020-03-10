
#GSE39582 dataset

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582
gds <- getGEO("GSE39582", destdir = "./Data")
gds <- gds[[1]]
destdir = "./Data"
gds<- getGEO(filename=paste(destdir,"GSE39582_series_matrix.txt.gz",sep="/"))
eset1 <- gds

primarysite=pData(eset1)$"tumor.location:ch1"
primarysite[primarysite=="N/A"]=NA
primarysite[!is.na(primarysite)]="co"

age= as.numeric(pData(eset1)$"age.at.diagnosis (year):ch1")

esetClinic= data.frame(primarysite=primarysite,age=age, days_to_recurrence=as.numeric(pData(eset1)$"rfs.delay:ch1")*30,rfs_status=pData(eset1)$"rfs.event:ch1",days_to_os=as.numeric(pData(eset1)$"os.delay:ch1")*30,os_status=pData(eset1)$"os.event:ch1", stageall=pData(eset1)$"tnm.stage:ch1",batch="GSE39582",row.names=rownames(pData(eset1)))

pData(eset1)=esetClinic

#GSE14333 dataset

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14333
gds <- getGEO("GSE14333", destdir = "./Data")
gds <- gds[[1]]
destdir = "./Data"
gds<- getGEO(filename=paste(destdir,"GSE14333_series_matrix.txt.gz",sep="/"))
eset2 <- gds

GSE14333Clin=read.table(text = gsub(":", ";", pData(eset2)$"Location:ch1"),sep=";")
primarysite=as.character(GSE14333Clin[,1])
primarysite[primarysite=="Rectum"]="re"
primarysite[primarysite=="Colon"| primarysite=="Left"| primarysite=="Right"]="co"
primarysite[primarysite==""]=NA

days_to_recurrence=as.numeric(as.character(GSE14333Clin[,"V9"]))*30
rfs_status= abs(as.numeric(as.character(GSE14333Clin[,"V11"]))-1)

age= as.numeric(GSE14333Clin[,"V5"])

esetClinic= data.frame(primarysite=primarysite,age=age,days_to_recurrence=days_to_recurrence,rfs_status=rfs_status,days_to_os=NA,os_status=NA, stageall=as.numeric(GSE14333Clin[,"V3"]),batch="GSE14333",row.names=rownames(pData(eset2)))
pData(eset2)=esetClinic

#GSE17536 dataset

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17536
gds <- getGEO("GSE17536", destdir = "./Data")
gds <- gds[[1]]
destdir = "./Data"
gds<- getGEO(filename=paste(destdir,"GSE17536_series_matrix.txt.gz",sep="/"))
eset3 <- gds

rfs_status= as.numeric(pData(eset3)$"dfs_event (disease free survival; cancer recurrence):ch1"=="recurrence")
days_to_recurrence=as.numeric(pData(eset3)$"dfs_time:ch1")*30

os_status= as.numeric(pData(eset3)$"overall_event (death from any cause):ch1"=="death")
days_to_os=as.numeric(pData(eset3)$"overall survival follow-up time:ch1")*30

age=as.numeric(pData(eset3)$"age:ch1")

esetClinic= data.frame(primarysite=NA,age=age,days_to_recurrence=days_to_recurrence,rfs_status=rfs_status,days_to_os=days_to_os,os_status=os_status, stageall=as.numeric(pData(eset3)$"ajcc_stage:ch1"),batch="GSE17536",row.names=rownames(pData(eset3)))
pData(eset3)=esetClinic

#GSE33113 dataset

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33113
gds <- getGEO("GSE33113", destdir = "./Data")
gds <- gds[[1]]
destdir = "./Data"
gds<- getGEO(filename=paste(destdir,"GSE33113_series_matrix.txt.gz",sep="/"))
eset4 <- gds

eset4= eset4[,pData(eset4)$"tissue:ch1"=="primary tumor resection"]
rfs_status= as.numeric(pData(eset4)$"meta or recurrence within 3 years:ch1"=="yes")
days_to_recurrence=as.numeric(pData(eset4)$"time to meta or recurrence:ch1")
stageall= 2

age=as.numeric(gsub(",",".",pData(eset4)$"age at diagnosis:ch1"))

esetClinic= data.frame(primarysite=NA,age=age,days_to_recurrence=days_to_recurrence,rfs_status=rfs_status,days_to_os=NA,os_status=NA, stageall=stageall,batch="GSE33113",row.names=pData(eset4)$title)
sampleNames(eset4)=as.character(pData(eset4)$title)
pData(eset4)=esetClinic


#GSE37892 dataset

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37892
gds <- getGEO("GSE37892", destdir = "./Data")
gds <- gds[[1]]
destdir = "./Data"
gds<- getGEO(filename=paste(destdir,"GSE37892_series_matrix.txt.gz",sep="/"))
eset5 <- gds

surgerytime=as.Date(pData(eset5)$"date at surgery:ch1")
mettime=  pData(eset5)$"date at distant metastasis:ch1"
mettime[pData(eset5)$"date at distant metastasis:ch1"=="NA"]=NA
mettime=as.Date(mettime)
lastcontacttime=as.Date(pData(eset5)$"date at last contact:ch1")
rfs_status=as.numeric(!is.na(mettime))
days_to_recurrence= as.numeric( lastcontacttime-surgerytime)
days_to_recurrence[rfs_status==1]=as.numeric(lastcontacttime[rfs_status==1]-mettime[rfs_status==1])

age=as.numeric(pData(eset5)$"age at diagnosis:ch1")
esetClinic= data.frame(primarysite="co",age=age,days_to_recurrence=days_to_recurrence,rfs_status=rfs_status,days_to_os=NA,os_status=NA, stageall=pData(eset5)$"Stage:ch1",batch="GSE37892",row.names=rownames(pData(eset5)))
pData(eset5)=esetClinic


downl= c("GSE39582","GSE14333","GSE17536","GSE33113","GSE37892")
esets= list(eset1,eset2,eset3,eset4,eset5)
names(esets)= downl

#Add consensus cluster to clinical data

#Download "cms_labels_public_all.txt" from "https://www.synapse.org/#!Synapse:syn4978511"
#https://www.synapse.org/#!Synapse:syn2623706/wiki/67246
coadCLASS=read.table("./Data/cms_labels_public_all.txt",sep="\t",header=TRUE)

for(i in downl){

require(dplyr)

pData(esets[[i]])$sample=rownames(pData(esets[[i]]))
P=inner_join(pData(esets[[i]]), coadCLASS)
rownames(P)=P$sample
P$sample=NULL
P$dataset=NULL
colnames(P)[which(colnames(P)=="CMS_final_network_plus_RFclassifier_in_nonconsensus_samples")]="CMS_Final"

esets[[i]]= esets[[i]][,sampleNames(esets[[i]]) %in% rownames(P)]
pData(esets[[i]])=P

#Check gene names!

require(hgu133plus2.db)

ID     <- fData(esets[[i]])$ID
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("ENTREZID"), keytype="PROBEID")
fData(esets[[i]])$EntrezGene.ID=out
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("SYMBOL"), keytype="PROBEID")
fData(esets[[i]])$gene=out

#Drop NA genes
Out=which(is.na(fData(esets[[i]])$EntrezGene.ID))
esets[[i]]=esets[[i]][-Out,]

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


esets[[i]]=expandProbesets(esets[[i]])


 #Apply robust linear scaling
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3283537/#bib61
Expression= apply(exprs(esets[[i]]),1,genefu::rescale,na.rm=TRUE,q=0.05)
Expression=t(Expression)
if(anyNA(Expression)){
NAs=apply(Expression,1,anyNA)
Expression[NAs,]=0}
exprs(esets[[i]])=Expression

esets[[i]]$rfs_status=as.numeric(as.character(esets[[i]]$rfs_status))
esets[[i]]$os_status=as.numeric(as.character(esets[[i]]$os_status))

ct=censor.time(pData(esets[[i]])$days_to_recurrence, pData(esets[[i]])$rfs_status, time.cens = 1825 )
pData(esets[[i]])$time=ct$surv.time.cens
pData(esets[[i]])$status=ct$surv.event.cens

#Remove NAs in survival and stage 4 samples
esets[[i]]=esets[[i]][,!c(is.na(pData(esets[[i]])$time) | is.na(pData(esets[[i]])$status) | pData(esets[[i]])$stageall==4)]

#Reformulate fData
fData(esets[[i]])=fData(esets[[i]])[,c("ID","EntrezGene.ID","gene")]

}

rm("P")
detach("package:hgu133plus2.db", unload=TRUE)






#TCGA data

library(TCGAbiolinks)
setwd("./Data")


query <- GDCquery(project = "TCGA-COAD",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification", 
                 sample.type=c("Primary solid Tumor"),
                 workflow.type = "HTSeq - Counts",
                 legacy=FALSE)


GDCdownload(query)

data <- GDCprepare(query, save = TRUE, 
                   save.filename = "RNA_COAD.rda",
                   remove.files.prepared = F)


query <- GDCquery(project = "TCGA-READ",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  sample.type=c("Primary solid Tumor"),
                  workflow.type = "HTSeq - Counts",
                  legacy=FALSE)


GDCdownload(query)

data <- GDCprepare(query, save = TRUE, 
                   save.filename = "RNA_READ.rda",
                   remove.files.prepared = F)

#Open expression matrix
setwd('..')
load("./Data/RNA_READ.rda")
data1=data
load("./Data/RNA_COAD.rda")
data2=data


library(SummarizedExperiment)
data=cbind(data1,data2)

data=data[,!duplicated(colData(data)$patient)]

#Primary site
primarysite=colData(data)$name
primarysite[primarysite=="Colon Adenocarcinoma"]="co"
primarysite[primarysite=="Rectum Adenocarcinoma"]="re"

#Stages
stageall=colData(data)$tumor_stage

isna= is.na(colData(data)$tumor_stage)
notrep=colData(data)$tumor_stage=="not reported"
stageall[isna|notrep]=NA

stagei= colData(data)$tumor_stage=="stage i"
stageia= colData(data)$tumor_stage=="stage ia"
stageall[stagei|stageia]=1

stageii= colData(data)$tumor_stage=="stage ii"
stageiia= colData(data)$tumor_stage=="stage iia"
stageiib= colData(data)$tumor_stage=="stage iib"
stageiic= colData(data)$tumor_stage=="stage iic"
stageall[stageii|stageiia|stageiib|stageiic]=2

stageiii= colData(data)$tumor_stage=="stage iii"
stageiiia= colData(data)$tumor_stage=="stage iiia"
stageiiib= colData(data)$tumor_stage=="stage iiib"
stageiiic= colData(data)$tumor_stage=="stage iiic"
stageall[stageiii|stageiiia|stageiiib|stageiiic]=3

stageiv=colData(data)$tumor_stage=="stage iv"
stageiva= colData(data)$tumor_stage=="stage iva"
stageivb= colData(data)$tumor_stage=="stage ivb"
stageall[stageiv|stageiva|stageivb]=4

rm(list=c("isna","notrep","stagei","stageia","stageii","stageiia","stageiib","stageiic","stageiii","stageiiia","stageiiib","stageiiic","stageiv","stageiva","stageivb"))

#Survival
colData(data)$status=as.numeric(colData(data)$vital_status=="dead")
colData(data)$time=NA
colData(data)$time[colData(data)$status==1 & !is.na(colData(data)$status)]= colData(data)$days_to_death[colData(data)$status==1 & !is.na(colData(data)$status)]
colData(data)$time[colData(data)$status==0 & !is.na(colData(data)$status)]= colData(data)$days_to_last_follow_up[colData(data)$status==0 & !is.na(colData(data)$status)]


ct=censor.time(colData(data)$time, colData(data)$status, time.cens = 1825 )
colData(data)$time=ct$surv.time.cens
colData(data)$status=ct$surv.event.cens

age= as.numeric(colData(data)$age_at_diagnosis/365)
esetClinic= data.frame(primarysite=primarysite,age, days_to_recurrence=NA,rfs_status=NA,days_to_os=colData(data)$time,os_status=as.numeric(colData(data)$status), stageall=as.numeric(stageall),batch="tcga",row.names=colData(data)$patient)
colData(data)=as(esetClinic, "DataFrame")

#Remove stage iv and uncknown stage

Out= is.na(colData(data)$stageall) | colData(data)$stageall == 4
data=data[,!Out]
#Voom normalization
library(limma)
library(edgeR)
Counts=assays(data)$HTSeq

#Filter low count genes
keep= filterByExpr(Counts)
data=data[keep,]
filteredCounts=Counts[keep,]
VoomCounts=voom(filteredCounts,plot=TRUE)

Expression= apply(VoomCounts,1,genefu::rescale,na.rm=TRUE,q=0.05)
Expression=t(Expression)
assays(data)$HTSeq=Expression

#Duplicated genes

Dup=rowData(data)[which(duplicated(rowData(data)$external_gene_name)),"external_gene_name"]
 Var= apply(VoomCounts,1,var)
drop=NULL
for(j in Dup){
  pos=which(rowData(data)$external_gene_name==j)
  drop= c(drop,pos[-which.max(Var[pos])])
  
}
data=data[-drop,]

rm(list=c("Counts","filteredCounts","VoomCounts"))


require(dplyr)

colData(data)$sample=rownames(colData(data))
P=inner_join(as.data.frame(colData(data)), coadCLASS)
rownames(P)=P$sample
P$sample=NULL
P$dataset=NULL
colnames(P)[which(colnames(P)=="CMS_final_network_plus_RFclassifier_in_nonconsensus_samples")]="CMS_Final"
P$time= P$days_to_os
P$status= P$os_status

#identical(rownames(P),rownames(colData(data))[rownames(colData(data)) %in% rownames(P)])
data=data[,rownames(colData(data)) %in% rownames(P)]

phenoData <- new("AnnotatedDataFrame",data=P)

#Reformulate fData
library(hgu133plus2.db)

ID     <- rowData(data)$ensembl_gene_id
Symbol <- rowData(data)$external_gene_name
Entrez <- mapIds(hgu133plus2.db, keys=Symbol, c("ENTREZID"), keytype="SYMBOL")



featureData= new("AnnotatedDataFrame",data=data.frame(ID=ID,EntrezGene.ID=Entrez,gene=Symbol,row.names=ID))

eset6= ExpressionSet(assayData=assays(data)$HTSeq,phenoData=phenoData,featureData=featureData)

eset6=eset6[!is.na(fData(eset6)$EntrezGene.ID),]
featureNames(eset6)=fData(eset6)$EntrezGene.ID

eset6=eset6[,!is.na(Surv(eset6$time,eset6$status))]

esets= append(esets,eset6)
downl=c(downl,"TCGA")
names(esets)=downl

inALL=Reduce(intersect, lapply(esets,featureNames))
inALLpData=Reduce(intersect, lapply(esets,function(x) colnames(pData(x))))
for(i in downl){
  esets[[i]]=esets[[i]][inALL,]
  pData(esets[[i]])=pData(esets[[i]])[,inALLpData]
}


#Check if dimensions match
dim(exprs(esets[[1]]))[1]==dim(exprs(esets[[6]]))[1]
#Check if rownames match
identical(featureNames(esets[[1]]), featureNames(esets[[6]]))

save(esets, file="./Data/RNA_CRC.rda")
