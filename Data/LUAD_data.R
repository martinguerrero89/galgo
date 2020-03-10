

#Centroids from Wilkerson et al.
#Download from "http://cancer.unc.edu/nhayes/publications/adenocarcinoma.2012/wilkerson.2012.LAD.predictor.centroids.csv.zip"

download.file("http://cancer.unc.edu/nhayes/publications/adenocarcinoma.2012/wilkerson.2012.LAD.predictor.centroids.csv.zip", "./Data/wilkerson.2012.LAD.predictor.centroids.csv.zip")
unzip("./Data/wilkerson.2012.LAD.predictor.centroids.csv.zip",exdir="./Data")
WilkCentroids= read.csv("./Data/wilkerson.2012.LAD.predictor.centroids.csv",row.names=1)

destdir = "./Data"

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50081
gds <- getGEO("GSE50081", destdir = destdir)
gds <- gds[[1]]

gds<- getGEO(filename=paste(destdir,"GSE50081_series_matrix.txt.gz",sep="/"))
eset1 <- gds

#Rename "adenocarcinoma" value

pData(eset1)$"histology:ch1"[pData(eset1)$"histology:ch1"=="adenocarcinoma"]="ADC"
histology= pData(eset1)$"histology:ch1"

RFS_time= as.numeric(pData(eset1)$"disease-free survival time:ch1")*365
RFS_status= pData(eset1)$"recurrence:ch1"
RFS_status[pData(eset1)$"recurrence:ch1"=="Y"]= 1
RFS_status[pData(eset1)$"recurrence:ch1"=="N"]= 0
RFS_status[pData(eset1)$"recurrence:ch1"=="U"]= NA
RFS_status=as.numeric(RFS_status)

OS_time= as.numeric(pData(eset1)$"survival time:ch1")*365
OS_status= pData(eset1)$"status:ch1"
OS_status[pData(eset1)$"status:ch1"=="dead"]= 1
OS_status[pData(eset1)$"status:ch1"=="alive"]= 0
OS_status=as.numeric(OS_status)

esetClinic= data.frame(histology=histology,age=as.numeric(pData(eset1)$"age:ch1"),days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=pData(eset1)$"t-stage:ch1",batch="GSE50081",row.names=rownames(pData(eset1)))

pData(eset1)=esetClinic

require(hgu133plus2.db)

ID     <- fData(eset1)$ID
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("ENTREZID"), keytype="PROBEID")
fData(eset1)$EntrezGene.ID=out
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("SYMBOL"), keytype="PROBEID")
fData(eset1)$gene=out

#Drop NA genes
Out=which(is.na(fData(eset1)$EntrezGene.ID))
eset1=eset1[-Out,]

#Drop duplicated genes (keep genes with highest variance)
Dup=unique(fData(eset1)[which(duplicated(fData(eset1)$EntrezGene.ID)),"EntrezGene.ID"])
Var= apply(exprs(eset1),1,var)
drop=NULL
for(j in Dup){
  pos=which(fData(eset1)$EntrezGene.ID==j)
  drop= c(drop,pos[-which.max(Var[pos])])
}
eset1=eset1[-drop,]
featureNames(eset1) <- fData(eset1)$EntrezGene.ID


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



eset1=expandProbesets(eset1)
fData(eset1)=fData(eset1)[,c("ID","EntrezGene.ID","gene")]


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210
gds <- getGEO("GSE31210", destdir = destdir)
gds <- gds[[1]]
gds<- getGEO(filename=paste(destdir,"GSE31210_series_matrix.txt.gz",sep="/"))
eset2 <- gds

pData(eset2)$"description"=as.character(pData(eset2)$"description")
pData(eset2)$"description"[pData(eset2)$"description"=="Gene expression data from primary lung ADC"]="ADC"
histology= pData(eset2)$"description"

RFS_time= as.numeric(pData(eset2)$"days before relapse/censor:ch1")
RFS_status= pData(eset2)$"relapse:ch1"
RFS_status[pData(eset2)$"relapse:ch1"=="relapsed"]= 1
RFS_status[pData(eset2)$"relapse:ch1"=="not relapsed"]= 0
RFS_status=as.numeric(RFS_status)

OS_time= as.numeric(pData(eset2)$"days before death/censor:ch1")
OS_status= pData(eset2)$"death:ch1"
OS_status[pData(eset2)$"death:ch1"=="dead"]= 1
OS_status[pData(eset2)$"death:ch1"=="alive"]= 0
OS_status=as.numeric(OS_status)

stageall=pData(eset2)$"pathological stage:ch1"
stageall[pData(eset2)$"pathological stage:ch1"=="IA"]=1
stageall[pData(eset2)$"pathological stage:ch1"=="IB"]=1
stageall[pData(eset2)$"pathological stage:ch1"=="II"]=2
stageall=as.numeric(stageall)

esetClinic= data.frame(histology=histology,age=as.numeric(pData(eset2)$"age (years):ch1"),days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=stageall,batch="GSE31210",row.names=rownames(pData(eset2)))

pData(eset2)=esetClinic

require(hgu133plus2.db)

ID     <- fData(eset2)$ID
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("ENTREZID"), keytype="PROBEID")
fData(eset2)$EntrezGene.ID=out
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("SYMBOL"), keytype="PROBEID")
fData(eset2)$gene=out

#Drop NA genes
Out=which(is.na(fData(eset2)$EntrezGene.ID))
eset2=eset2[-Out,]

#Drop duplicated genes (keep genes with highest variance)
Dup=unique(fData(eset2)[which(duplicated(fData(eset2)$EntrezGene.ID)),"EntrezGene.ID"])
Var= apply(exprs(eset2),1,var)
drop=NULL
for(j in Dup){
  pos=which(fData(eset2)$EntrezGene.ID==j)
  drop= c(drop,pos[-which.max(Var[pos])])
}
eset2=eset2[-drop,]
featureNames(eset2) <- fData(eset2)$EntrezGene.ID
eset2=expandProbesets(eset2)
fData(eset2)=fData(eset2)[,c("ID","EntrezGene.ID","gene")]


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188
gds <- getGEO("GSE19188", destdir = destdir)
gds <- gds[[1]]
gds<- getGEO(filename=paste(destdir,"GSE19188_series_matrix.txt.gz",sep="/"))
eset3 <- gds

histology= pData(eset3)$"cell type:ch1"

RFS_time= NA
RFS_status= NA

OS_time= as.numeric(pData(eset3)$"overall survival:ch1")*30
OS_status= pData(eset3)$"status:ch1"
OS_status[pData(eset3)$"status:ch1"=="deceased"]= 1
OS_status[pData(eset3)$"status:ch1"=="alive"]= 0
OS_status=as.numeric(OS_status)

#Tumor stage not available but there are very few stage IV according to the original paper
stageall=0

esetClinic= data.frame(histology=histology,age=NA,days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=stageall,batch="GSEGSE19188",row.names=rownames(pData(eset3)))

pData(eset3)=esetClinic

ID     <- fData(eset3)$ID
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("ENTREZID"), keytype="PROBEID")
fData(eset3)$EntrezGene.ID=out
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("SYMBOL"), keytype="PROBEID")
fData(eset3)$gene=out

#Drop NA genes
Out=which(is.na(fData(eset3)$EntrezGene.ID))
eset3=eset3[-Out,]

#Drop duplicated genes (keep genes with highest variance)
Dup=unique(fData(eset3)[which(duplicated(fData(eset3)$EntrezGene.ID)),"EntrezGene.ID"])
Var= apply(exprs(eset3),1,var)
drop=NULL
for(j in Dup){
  pos=which(fData(eset3)$EntrezGene.ID==j)
  drop= c(drop,pos[-which.max(Var[pos])])
}
eset3=eset3[-drop,]
featureNames(eset3) <- fData(eset3)$EntrezGene.ID
eset3=expandProbesets(eset3)
fData(eset3)=fData(eset3)[,c("ID","EntrezGene.ID","gene")]


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37745
gds <- getGEO("GSE37745", destdir = destdir)
gds <- gds[[1]]
gds<- getGEO(filename=paste(destdir,"GSE37745_series_matrix.txt.gz",sep="/"))
eset4 <- gds

pData(eset4)$"histology:ch1"[pData(eset4)$"histology:ch1"=="adeno"]="ADC"
histology= pData(eset4)$"histology:ch1"

RFS_time= as.numeric(pData(eset4)$"days to recurrence / to last visit:ch1")
RFS_status= pData(eset4)$"recurrence:ch1"
RFS_status[pData(eset4)$"recurrence:ch1"=="yes"]= 1
RFS_status[pData(eset4)$"recurrence:ch1"=="no"]= 0
RFS_status=as.numeric(RFS_status)

OS_time= as.numeric(pData(eset4)$"days to determined death status:ch1")
OS_status= pData(eset4)$"dead:ch1"
OS_status[pData(eset4)$"dead:ch1"=="yes"]= 1
OS_status[pData(eset4)$"dead:ch1"=="no"]= 0
OS_status=as.numeric(OS_status)

stageall=pData(eset4)$"tumor stage:ch1"
stageall[pData(eset4)$"tumor stage:ch1"=="1a"]=1
stageall[pData(eset4)$"tumor stage:ch1"=="1b"]=1
stageall[pData(eset4)$"tumor stage:ch1"=="2a"]=2
stageall[pData(eset4)$"tumor stage:ch1"=="2b"]=2
stageall[pData(eset4)$"tumor stage:ch1"=="3a"]=3
stageall[pData(eset4)$"tumor stage:ch1"=="3b"]=3
stageall=as.numeric(stageall)

esetClinic= data.frame(histology=histology,age=as.numeric(pData(eset4)$"age:ch1"),days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=stageall,batch="GSE37745",row.names=rownames(pData(eset4)))

pData(eset4)=esetClinic

ID     <- fData(eset4)$ID
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("ENTREZID"), keytype="PROBEID")
fData(eset4)$EntrezGene.ID=out
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("SYMBOL"), keytype="PROBEID")
fData(eset4)$gene=out

#Drop NA genes
Out=which(is.na(fData(eset4)$EntrezGene.ID))
eset4=eset4[-Out,]

#Drop duplicated genes (keep genes with highest variance)
Dup=unique(fData(eset4)[which(duplicated(fData(eset4)$EntrezGene.ID)),"EntrezGene.ID"])
Var= apply(exprs(eset4),1,var)
drop=NULL
for(j in Dup){
  pos=which(fData(eset4)$EntrezGene.ID==j)
  drop= c(drop,pos[-which.max(Var[pos])])
}
eset4=eset4[-drop,]
featureNames(eset4) <- fData(eset4)$EntrezGene.ID
eset4=expandProbesets(eset4)
fData(eset4)=fData(eset4)[,c("ID","EntrezGene.ID","gene")]


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14814
gds <- getGEO("GSE14814", destdir = destdir)
gds <- gds[[1]]
gds<- getGEO(filename=paste(destdir,"GSE14814_series_matrix.txt.gz",sep="/"))
eset5 <- gds

histology= pData(eset5)$"Histology type:ch1"

RFS_time= as.numeric(pData(eset5)$"DSS time:ch1")*365
RFS_status= pData(eset5)$"DSS status:ch1"
RFS_status[pData(eset5)$"DSS status:ch1"=="dead"]= 1
RFS_status[pData(eset5)$"DSS status:ch1"=="Dead"]= 1
RFS_status[pData(eset5)$"DSS status:ch1"=="Alive"]= 0
RFS_status=as.numeric(RFS_status)

OS_time= as.numeric(pData(eset5)$"OS time:ch1")*365
OS_status= pData(eset5)$"OS status:ch1"
OS_status[pData(eset5)$"OS status:ch1"=="Dead"]= 1
OS_status[pData(eset5)$"OS status:ch1"=="Alive"]= 0
OS_status=as.numeric(OS_status)

stageall=pData(eset5)$"Stage:ch1"
stageall[pData(eset5)$"Stage:ch1"=="I"]=1
stageall[pData(eset5)$"Stage:ch1"=="1B"]=1
stageall[pData(eset5)$"Stage:ch1"=="II"]=2
stageall[pData(eset5)$"Stage:ch1"=="2A"]=2
stageall[pData(eset5)$"Stage:ch1"=="2B"]=2
stageall=as.numeric(stageall)

esetClinic= data.frame(histology=histology,age=as.numeric(pData(eset5)$"age:ch1"),days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=stageall,batch="GSE14814",row.names=rownames(pData(eset5)))

pData(eset5)=esetClinic

require(hgu133a.db)

ID     <- fData(eset5)$ID
out <- mapIds(hgu133a.db, keys=as.character(ID), c("ENTREZID"), keytype="PROBEID")
fData(eset5)$EntrezGene.ID=out
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("SYMBOL"), keytype="PROBEID")
fData(eset5)$gene=out

#Drop NA genes
Out=which(is.na(fData(eset5)$EntrezGene.ID))
eset5=eset5[-Out,]

#Drop duplicated genes (keep genes with highest variance)
Dup=unique(fData(eset5)[which(duplicated(fData(eset5)$EntrezGene.ID)),"EntrezGene.ID"])
Var= apply(exprs(eset5),1,var)
drop=NULL
for(j in Dup){
  pos=which(fData(eset5)$EntrezGene.ID==j)
  drop= c(drop,pos[-which.max(Var[pos])])
}
eset5=eset5[-drop,]
featureNames(eset5) <- fData(eset5)$EntrezGene.ID
eset5=expandProbesets(eset5)
fData(eset5)=fData(eset5)[,c("ID","EntrezGene.ID","gene")]


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465
gds <- getGEO("GSE68465", destdir = destdir)
gds <- gds[[1]]
gds<- getGEO(filename=paste(destdir,"GSE68465_series_matrix.txt.gz",sep="/"))
eset6 <- gds

histology= pData(eset6)$"disease_state:ch1"
histology[histology=="Lung Adenocarcinoma"]="ADC"


RFS_time= as.numeric(pData(eset6)$"months_to_first_progression:ch1")*30
RFS_time[is.na(RFS_time)]=(as.numeric(pData(eset6)$"months_to_last_contact_or_death:ch1")*30)[is.na(RFS_time)]
RFS_status= pData(eset6)$"first_progression_or_relapse:ch1"
RFS_status[pData(eset6)$"first_progression_or_relapse:ch1"=="Yes"]= 1
RFS_status[pData(eset6)$"first_progression_or_relapse:ch1"=="No"]= 0
RFS_status=as.numeric(RFS_status)

OS_time= as.numeric(pData(eset6)$"months_to_last_contact_or_death:ch1")*30
OS_status= pData(eset6)$"vital_status:ch1"
OS_status[pData(eset6)$"vital_status:ch1"=="Dead"]= 1
OS_status[pData(eset6)$"vital_status:ch1"=="Alive"]= 0
OS_status=as.numeric(OS_status)

stageall=pData(eset6)$"disease_stage:ch1"
stageall[pData(eset6)$"disease_stage:ch1"=="pN0pT1"]=1
stageall[pData(eset6)$"disease_stage:ch1"=="pN0pT2"]=1 
stageall[pData(eset6)$"disease_stage:ch1"=="pN0pT3"]=2
stageall[pData(eset6)$"disease_stage:ch1"=="pN0pT4"]=3
stageall[pData(eset6)$"disease_stage:ch1"=="pN1pT1"]=2
stageall[pData(eset6)$"disease_stage:ch1"=="pN1pT2"]=2 
stageall[pData(eset6)$"disease_stage:ch1"=="pN1pT3"]=3
stageall[pData(eset6)$"disease_stage:ch1"=="pN1pT4"]=3
stageall[pData(eset6)$"disease_stage:ch1"=="pN2pT1"]=3
stageall[pData(eset6)$"disease_stage:ch1"=="pN2pT2"]=3
stageall[pData(eset6)$"disease_stage:ch1"=="pN2pT3"]=3
stageall[pData(eset6)$"disease_stage:ch1"=="pN2pT4"]=3
stageall[pData(eset6)$"disease_stage:ch1"=="pNXpT1"]=NA
stageall[pData(eset6)$"disease_stage:ch1"=="pp"]=NA
stageall=as.numeric(stageall)

esetClinic= data.frame(histology=histology,age=as.numeric(pData(eset6)$"age:ch1"),days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=stageall,batch="GSE68465",row.names=rownames(pData(eset6)))

pData(eset6)=esetClinic

require(hgu133a.db)

ID     <- fData(eset6)$ID
out <- mapIds(hgu133a.db, keys=as.character(ID), c("ENTREZID"), keytype="PROBEID")
fData(eset6)$EntrezGene.ID=out
out <- mapIds(hgu133plus2.db, keys=as.character(ID), c("SYMBOL"), keytype="PROBEID")
fData(eset6)$gene=out

#Drop NA genes
Out=which(is.na(fData(eset6)$EntrezGene.ID))
eset6=eset6[-Out,]

#Drop duplicated genes (keep genes with highest variance)
Dup=unique(fData(eset6)[which(duplicated(fData(eset6)$EntrezGene.ID)),"EntrezGene.ID"])
Var= apply(exprs(eset6),1,var)
drop=NULL
for(j in Dup){
  pos=which(fData(eset6)$EntrezGene.ID==j)
  drop= c(drop,pos[-which.max(Var[pos])])
}
eset6=eset6[-drop,]
featureNames(eset6) <- fData(eset6)$EntrezGene.ID
eset6=expandProbesets(eset6)
fData(eset6)=fData(eset6)[,c("ID","EntrezGene.ID","gene")]



#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")
library(TCGAbiolinks)
setwd(destdir)


query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  sample.type=c("Primary solid Tumor"),
                  workflow.type = "HTSeq - Counts",
                  legacy=FALSE)


GDCdownload(query)

data <- GDCprepare(query, save = TRUE, 
                   save.filename = "RNA_LUAD.rda",
                   remove.files.prepared = F)

#Open expression matrix
setwd('..')
load("./Data/RNA_LUAD.rda")
library(SummarizedExperiment)

data=data[,!duplicated(colData(data)$patient)]
data=data[!duplicated(rowData(data)$ensembl_gene_id),]


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


histology= "ADC"

RFS_time= NA
RFS_status= NA


OS_status=colData(data)$"vital_status"
OS_status[colData(data)$"vital_status"=="dead"]= 1
OS_status[colData(data)$"vital_status"=="alive"]= 0
OS_status=as.numeric(OS_status)

OS_time= as.numeric(colData(data)$"days_to_death")
OS_time[is.na(OS_time)]=as.numeric(colData(data)$"days_to_last_follow_up")[is.na(OS_time)]

stageall=colData(data)$"tumor_stage"
stageall[colData(data)$"tumor_stage"=="stage i"]=1
stageall[colData(data)$"tumor_stage"=="stage ia"]=1
stageall[colData(data)$"tumor_stage"=="stage ib"]=1
stageall[colData(data)$"tumor_stage"=="stage ii"]=2
stageall[colData(data)$"tumor_stage"=="stage ii"]=2
stageall[colData(data)$"tumor_stage"=="stage iia"]=2
stageall[colData(data)$"tumor_stage"=="stage iib"]=2
stageall[colData(data)$"tumor_stage"=="stage iiia"]=3
stageall[colData(data)$"tumor_stage"=="stage iiib"]=3
stageall[colData(data)$"tumor_stage"=="stage iv"]=4
stageall=as.numeric(stageall)

esetClinic= data.frame(histology=histology,age=as.numeric(colData(data)$"age_at_diagnosis")/365,days_to_recurrence=RFS_time,rfs_status=RFS_status,days_to_os=OS_time,os_status=OS_status, stageall=stageall,batch="TCGA",row.names=colData(data)$patient)

phenoData <- new("AnnotatedDataFrame",data=esetClinic)

#Reformulate fData
library(hgu133plus2.db)

ID     <- rowData(data)$ensembl_gene_id
Symbol <- rowData(data)$external_gene_name
Entrez <- mapIds(hgu133plus2.db, keys=Symbol, c("ENTREZID"), keytype="SYMBOL")


featureData= new("AnnotatedDataFrame",data=data.frame(ID=ID,EntrezGene.ID=Entrez,gene=Symbol,row.names=ID))
eset7= ExpressionSet(assayData=VoomCounts$E,phenoData=phenoData,featureData=featureData)


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



downl= c("GSE50081","GSE31210","GSE19188","GSE37745","GSE14814","GSE68465","TCGA")
esets=list(eset1,eset2,eset3,eset4,eset5,eset6,eset7)
names(esets)=downl


for(i in downl){

require(dplyr)

pData(esets[[i]])$sample=rownames(pData(esets[[i]]))

 #Apply robust linear scaling
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3283537/#bib61
Expression= apply(exprs(esets[[i]]),1,genefu::rescale,na.rm=TRUE,q=0.05)
Expression=t(Expression)
if(anyNA(Expression)){
NAs=apply(Expression,1,anyNA)
Expression[NAs,]=0}
exprs(esets[[i]])=Expression


ct=censor.time(pData(esets[[i]])$days_to_os, pData(esets[[i]])$os_status, time.cens = 1825 )
pData(esets[[i]])$time=ct$surv.time.cens
pData(esets[[i]])$status=ct$surv.event.cens

#Keep only ADC
esets[[i]]=esets[[i]][,pData(esets[[i]])$histology=="ADC"]

#Remove NAs in survival and stage 4 samples
esets[[i]]=esets[[i]][,!c(is.na(Surv(pData(esets[[i]])$time, pData(esets[[i]])$status)) | c(pData(esets[[i]])$stageall==4)| is.na(pData(esets[[i]])$stageall)) ]


#Classify into Wilkerson et al. subtypes
featureNames(esets[[i]])= fData(esets[[i]])$EntrezGene.ID


wilkGenes= rownames(WilkCentroids)
SubCentroid= as.matrix(WilkCentroids[wilkGenes %in% fData(esets[[i]])$gene,])
esets[[i]]$Wilk.Subtype= classify(exprs(esets[[i]])[match(rownames(SubCentroid),fData(esets[[i]])$gene),],SubCentroid)
}


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

save(esets, file="/home/mguerrero/Genetic_alg/App_FINAL/Data/RNA_LUAD_all.rda")
                                    
