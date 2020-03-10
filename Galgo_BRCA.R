
load("./Data/RNA_BRCA.rda") #Directory where RNA_BRCA.rda was saved
downl=names(esets)

##Galgo Hyperparameters

population= 300                   #Number of individuals to evaluate
generations=500                 #Number of generations
nCV=5                             #Number of crossvalidations for function "crossvalidation"
GPU= TRUE                         # to use gpuR
distancetype="pearson"          #Options are: "pearson","uncentered","spearman","euclidean"
TournamentSize=2
period=5475



trainSet=esets[["metabric"]]
prob_matrix= exprs(trainSet)
clinical=pData(trainSet)
OS=Surv(time=clinical$time,event=clinical$status)
chrom_length= nrow(prob_matrix)   #length of chromosome



h="BRCA"
dir.create(paste0("./Results/Results_",h),recursive=TRUE)
resultdir=paste0("./Results/Results_",h,"/")


output= galgoR::galgo(generations = generations, population = population,prob_matrix = prob_matrix, OS=OS,usegpu=GPU,nCV= nCV,
                   distancetype=distancetype, TournamentSize=TournamentSize, period=period,res_dir=resultdir)


#file= paste0("results",generations,".rda")
#load(paste0(resultdir,file))
X1=output$Solutions
PARETO=output$ParetoFront
dots=do.call(rbind,PARETO)

rownames(X1)=paste("result",1:nrow(X1),sep=".")

#Plot pareto front
plot(dots)
lines(X1[names(which(X1[order(X1[,(chrom_length+2)]),"rnkIndex"]==1)),(chrom_length+2):(chrom_length+3)],col="red")


RES=data.frame(solution=as.character(),k=as.numeric(),ngenes=as.numeric(),trainSil=as.numeric(),trainC=as.numeric(),stringsAsFactors = FALSE)

if(GPU==TRUE){
  mydist= galgoR::calculate_distance_pearson_gpu
}else{
    mydist= galgoR::calculate_distance_pearson_cpu}


for(i in 1:sum(X1[,"rnkIndex"]==1)){

  RESULT=names(which(X1[order(X1[,(chrom_length+2)]),"rnkIndex"]==1))[i]
  R=as.logical(X1[RESULT,1:chrom_length])
  k=X1[RESULT,"k"]
  D=mydist(prob_matrix[R,])
  
  
  hsp_class=galgoR::cluster_algorithm(D,k)
  C=galgoR::kcentroid(prob_matrix[R,],hsp_class$cluster)
  hsp_class= as.factor(galgoR::classify(prob_matrix[R,],C,method=distancetype))
  hsp_classdf= as.data.frame(hsp_class)
  
  mysurv <- OS
  tumortotal <- survfit(mysurv~ hsp_classdf$hsp_class)
  totalsdf <- survdiff(mysurv~ hsp_classdf$hsp_class)
  tumortotalpval <- 1 - pchisq(totalsdf$chisq, length(totalsdf$n) - 1)
  tumortotalpval <- format(tumortotalpval, digits=4)
    
  coxsimple=coxph(mysurv~hsp_class,data=hsp_classdf)
  a=concordance.index(predict(coxsimple),surv.time=OS[,1],surv.event=OS[,2],outx=FALSE)$c.index
    
  c=mean(silhouette(as.numeric(hsp_class),D)[,3])
  
  row= c(RESULT,k,sum(R),c,a)
  RES[nrow(RES)+1,]=row
  print(row)
  
  }
  

  testSol=NULL
  for(i in unique(RES$k)){
  testSol=c(testSol,RES[RES$k==i,"solution"][which.max(RES[RES$k==i,"trainC"])])
    }

  
  galgo=NULL
  CentroidsList=list()
  for(j in testSol){
  RESULT=j
  R=as.logical(X1[RESULT,1:chrom_length])
  k=X1[RESULT,"k"]
  D=mydist(prob_matrix[R,])
  hsp_class=galgoR::cluster_algorithm(D,k)$cluster
  hsp_classdf= data.frame(hsp_class=as.factor(hsp_class))
  coxsimple=coxph(mysurv~hsp_class,data=hsp_classdf)
  C=galgoR::kcentroid(prob_matrix[R,],hsp_class)
  labels=rank(c(0,coxsimple$coefficients))
  colnames(C)= labels
  CentroidsList[[j]]=C
  name= paste(k,j,sep="_")
  galgo=c(galgo,name)
     for(i in downl){
    hsp_class= as.factor(galgoR::classify(exprs(esets[[i]])[R,],C,method=distancetype))
  hsp_classdf= data.frame(hsp_class=hsp_class)
  pred= predict(coxsimple,hsp_classdf)
  levels(hsp_class)=labels
  hsp_class = factor(hsp_class,levels(hsp_class)[order(labels)])
  pData(esets[[i]])[,name]=hsp_class
  pData(esets[[i]])[,paste0(name,".Pred")]=pred
  
    }
  }
  

SSP= c("pam50","scm","intclust","AIMS")
  for(j in SSP){
  class=clinical[,j]
  classdf= data.frame(class=as.factor(class))
  coxsimple=coxph(OS~class,data=classdf)
     for(i in downl){
    class= pData(esets[[i]])[,j]
  classdf= data.frame(class=class)
  pred= predict(coxsimple,classdf)
  pData(esets[[i]])[,paste0(j,".Pred")]=pred
      }
  }
  
  
  
riskP=c("oncotypedx","endopredict","Mammaprint",paste0(c(galgo,SSP),".Pred")) 
resMatrix <- as.list(NULL)
for(i in downl){
res=pData(esets[[i]])
SURV=Surv(res[,"time"],res[,"status"])

for (Dat in riskP){
x=t(res[,Dat])
cindex <- t(apply(X=x, MARGIN=1, function(x, y) {
tt <- concordance.index(x=x, surv.time=y[,1], surv.event=y[,2], method="noether",outx=FALSE, na.rm=TRUE);
return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)) },
y=SURV))
resMatrix[[Dat]] <- rbind(resMatrix[[Dat]], cindex)
}
}

lapply(resMatrix,colMeans)

#Comparison of signatures
 ccmData <- tt <- rr <- NULL

for(i in downl){
CIlist <- t(apply(X=pData(esets[[i]])[,riskP], MARGIN=2, function(x, y, z) {
     tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", outx=FALSE, na.rm=TRUE);
     return(tt); }, y=pData(esets[[i]])$time, z=pData(esets[[i]])$status))
assign(i,CIlist)
}

for(i in 1:length(riskP)){
 tt <- NULL
 listOne <- list()
 for(j in downl){
 listOne[[j]]=get(j)[[i]]
 }
 for(z in 1:length(riskP)){
listTwo <- list()
 for(j in downl){
 listTwo[[j]]=get(j)[[z]]
 }
 rr <- cindex.comp.meta(list.cindex1=listOne, list.cindex2=listTwo)
 tt <- cbind(tt, rr$p.value) ##list(round(rr$p.value,5)))
}
ccmData <- rbind(ccmData, tt)
 }

ccmData <- as.data.frame(ccmData)
colnames(ccmData) <- riskP
rownames(ccmData) <- riskP


#CINDEX FOREST SUPERPLOT

resMatrixForest=resMatrix
for( i in riskP){
tt=combine.est(as.numeric(resMatrixForest[[i]][,"cindex"]),as.numeric(resMatrixForest[[i]][,"cindex.se"]))
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
resMatrixForest[[i]]=rbind(resMatrixForest[[i]],unlist(tt))

resMatrixForest[[i]]=cbind(resMatrixForest[[i]],dataset=c(downl,"Combined"))
resMatrixForest[[i]]=cbind(resMatrixForest[[i]],signature=i)
}

RR_data= as.data.frame(do.call(rbind,resMatrixForest),stringsAsFactors=FALSE)
RR_data$cindex=as.numeric(RR_data$cindex)
RR_data$lower=as.numeric(RR_data$lower)
RR_data$upper=as.numeric(RR_data$upper)
RR_data$dataset <- factor(RR_data$dataset, levels = c(downl,"Combined"))


p = ggplot(data=RR_data,
    aes(x = signature,y = cindex, ymin = lower, ymax = upper ))+
    geom_pointrange(aes(col=signature))+
    geom_hline(aes(fill=signature),yintercept =0.5, linetype=2)+
    xlab('Signature for each dataset')+ ylab("C.Index (95% Confidence Interval)")+
    geom_errorbar(aes(ymin=lower, ymax=upper,col=signature),width=0.5,cex=1)+ 
    facet_wrap(~dataset,strip.position="left",nrow=9,scales = "free_y") +
    theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
    coord_flip()
 p



Comb=combineTwoExpressionSet(esets[[1]],esets[[2]])
for(i in 3:length(downl[-which(downl=="metabric")])){
  Comb=combineTwoExpressionSet(Comb,esets[[i]])
}


finalSig= RR_data[RR_data$dataset=="Combined" & RR_data$signature %in% paste0(galgo,".Pred")[which.max(RR_data[RR_data$dataset=="Combined" & RR_data$signature %in% paste0(galgo,".Pred"),"cindex"])],"signature"]
finalSig= substr(finalSig, 1, nchar(finalSig)-5)


  tumortotal1 <- survfit(Surv(Comb$time,Comb$status)~ Comb$pam50)
  tumortotal2 <- survfit(Surv(Comb$time,Comb$status)~ Comb[[finalSig]])
  tumortotal1diff <- survdiff(Surv(Comb$time,Comb$status)~ Comb$pam50)
  tumortotal2diff <- survdiff(Surv(Comb$time,Comb$status)~ Comb[[finalSig]])
  tumortotal1pval<- pchisq(tumortotal1diff$chisq, length(tumortotal1diff$n) - 1, lower.tail = FALSE) 
  tumortotal2pval<- pchisq(tumortotal2diff$chisq, length(tumortotal2diff$n) - 1, lower.tail = FALSE) 
  
  
  SURV=list(PAM50=tumortotal1, GALGO=tumortotal2)
  COLS=c(1:8,10)
  par(cex=1.35, mar=c(3.8, 3.8, 2.5, 2.5) + 0.1)
  p=ggsurvplot(SURV,combine=TRUE,data=Comb,risk.table=TRUE,pval=TRUE,palette="dark2", title="Combined Set \n Breast survival comparison", surv.scale="percent", conf.int=FALSE, xlab="time (days)", ylab="survival(%)", xlim=c(0,period),break.time.by = 365, ggtheme = theme_minimal(), risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,censor=FALSE)
  print(p)
  
  tumortotal1 <- survfit(Surv(esets[["metabric"]]$time,esets[["metabric"]]$status)~ esets[["metabric"]]$pam50)
  tumortotal2 <- survfit(Surv(esets[["metabric"]]$time,esets[["metabric"]]$status)~ esets[["metabric"]][[finalSig]])
  tumortotal1diff <- survdiff(Surv(esets[["metabric"]]$time,esets[["metabric"]]$status)~ esets[["metabric"]]$pam50)
  tumortotal2diff <- survdiff(Surv(esets[["metabric"]]$time,esets[["metabric"]]$status)~ esets[["metabric"]][[finalSig]])
  tumortotal1pval<- pchisq(tumortotal1diff$chisq, length(tumortotal1diff$n) - 1, lower.tail = FALSE) 
  tumortotal2pval<- pchisq(tumortotal2diff$chisq, length(tumortotal2diff$n) - 1, lower.tail = FALSE) 
  
  
  
  SURV=list(PAM50=tumortotal1, GALGO=tumortotal2)
  par(cex=1.35, mar=c(3.8, 3.8, 2.5, 2.5) + 0.1)
  p=ggsurvplot(SURV,combine=TRUE,data=esets[["metabric"]],risk.table=TRUE,pval=TRUE,palette="dark2", title="Train Set \n Breast survival comparison", surv.scale="percent", conf.int=FALSE, xlab="time (days)", ylab="survival(%)", xlim=c(0,5500),break.time.by = 365, ggtheme = theme_minimal(), risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,censor=FALSE)
  print(p)


V=unlist(strsplit(finalSig,"_"))
R=as.logical(X1[V[2],1:chrom_length])
  k=as.numeric(V[1])
  

trainset="metabric"
TrainExprs=exprs(esets[[trainset]])
rownames(TrainExprs)= fData(esets[[trainset]])$gene  
TrainClass= pData(esets[[trainset]])[,riskP]


colfuncR <- colorRampPalette(rev(brewer.pal(11 , "Spectral" )))

OncotypeCol= colfuncR(100)[as.numeric(cut(as.numeric(TrainClass[,"oncotypedx"]),breaks=100))]
endoCol= colfuncR(100)[as.numeric(cut(as.numeric(TrainClass[,"endopredict"]),breaks=100))]
MammaCol= colfuncR(100)[as.numeric(cut(as.numeric(TrainClass[,"Mammaprint"]),breaks=100))]

IntClustCol=colfuncR(10)[as.numeric(as.factor(TrainClass[,"intclust.Pred"]))]
scmCol=colfuncR(4)[as.numeric(as.factor(TrainClass[,"scm.Pred"]))]
pam50Col=colfuncR(5)[as.numeric(as.factor(TrainClass[,"pam50.Pred"]))]
AIMSCol=colfuncR(5)[as.numeric(as.factor(TrainClass[,"AIMS.Pred"]))]

for(i in 1:length(galgo)){
 V= unlist(strsplit(galgo[i],"_"))
 n= as.numeric(V[1]) 
  assign(paste0("galgo",n), colfuncR(n)[as.numeric(as.factor(TrainClass[,paste(galgo[i],"Pred",sep=".")]))])
}

GalgoCols=do.call(cbind,mget(paste0("galgo",sub("\\_.*", "", galgo))))
SignatureCols= data.frame(AIMSCol,pam50Col,scmCol,IntClustCol,endoCol,MammaCol,OncotypeCol)
COLS=cbind(GalgoCols,SignatureCols)


CLASS=TrainClass[,paste(finalSig,"Pred",sep=".")]
Ord= order(CLASS)

source("./Functions/heatmap_3.R")

colfunc <- colorRampPalette(c("blue", "white", "red"))
breaks=seq(-2, 2, by=0.1) 
#now add outliers
breaks=append(breaks, 5)
breaks=append(breaks, -5, 0)

mydist2= function(x) mydist(t(x))

HM1=heatmap.3(TrainExprs[R,Ord],Colv= FALSE,hclustfun=hclust,dendrogram="none", distfun=mydist2,trace="none",scale="row",col=colfunc,breaks=breaks,ColSideColors =as.matrix(COLS[Ord,]),ColSideColorsSize=3)


#Testset (Comb)


TestExprs=ComBat(exprs(Comb),pData(Comb)$batch)
rownames(TestExprs)= fData(Comb)$gene  
TestClass= pData(Comb)[,riskP]

colfuncR <- colorRampPalette(rev(brewer.pal(11 , "Spectral" )))
OncotypeCol= colfuncR(100)[as.numeric(cut(as.numeric(TestClass[,"oncotypedx"]),breaks=100))]
endoCol= colfuncR(100)[as.numeric(cut(as.numeric(TestClass[,"endopredict"]),breaks=100))]
MammaCol= colfuncR(100)[as.numeric(cut(as.numeric(TestClass[,"Mammaprint"]),breaks=100))]

IntClustCol=colfuncR(10)[as.numeric(as.factor(TestClass[,"intclust.Pred"]))]
scmCol=colfuncR(4)[as.numeric(as.factor(TestClass[,"scm.Pred"]))]
pam50Col=colfuncR(5)[as.numeric(as.factor(TestClass[,"pam50.Pred"]))]
AIMSCol=colfuncR(5)[as.numeric(as.factor(TestClass[,"AIMS.Pred"]))]

for(i in 1:length(galgo)){
 V= unlist(strsplit(galgo[i],"_"))
 n= as.numeric(V[1]) 
  assign(paste0("galgo",n), colfuncR(n)[as.numeric(as.factor(TestClass[,paste(galgo[i],"Pred",sep=".")]))])
}

GalgoCols=do.call(cbind,mget(paste0("galgo",substr(galgo,1,1))))
SignatureCols= data.frame(AIMSCol,pam50Col,scmCol,IntClustCol,endoCol,MammaCol,OncotypeCol)
COLS=cbind(GalgoCols,SignatureCols)


CLASS=TestClass[,paste(finalSig,"Pred",sep=".")]
Ord= order(CLASS)

HM2=heatmap.3(TestExprs[R,Ord][rev(HM1$rowInd),],Colv= FALSE,Rowv=FALSE,hclustfun=hclust,dendrogram="none", distfun=mydist2,trace="none",scale="row",col=colfunc,breaks=breaks,ColSideColors =as.matrix(COLS[Ord,]),ColSideColorsSize=3)

#Save centroids and metanalisis

for(i in names(CentroidsList)){
write.xlsx(CentroidsList[[i]], file=paste(resultdir,"Breast_Galgo_Centroids.xlsx",sep="/"), sheetName=i,append=TRUE, row.names=TRUE)
}


write.xlsx(RR_data[,c(6,5,1:4)], file=paste(resultdir,"Breast_CI_Meta.xlsx",sep="/"), sheetName="CI",append=TRUE, row.names=FALSE)
write.xlsx(ccmData, file=paste(resultdir,"Breast_CI_Meta.xlsx",sep="/"), sheetName="Comparison",append=TRUE, row.names=TRUE)


#Gage pathway analysis


#data(kegg.gs)
data(go.sets.hs)
kegg.gs=go.sets.hs

CLASS=as.factor(TrainClass[,paste(finalSig,"Pred",sep=".")])
TrainExprs= exprs(esets[[trainset]])


for(i in 1:length(levels(CLASS))){
   ind= which(CLASS==levels(CLASS)[i])
   clust.kegg.p = gage(TrainExprs,ref= c(1:ncol(TrainExprs))[-ind], gsets=kegg.gs,compare="as.group")
   assign(paste0("clust.kegg.p_",i), clust.kegg.p)
  
##Sort and count signficant gene sets based on q- or p-value cutoffs:
   clust.kegg.p.sig<- sigGeneSet(clust.kegg.p, outname = paste0("clust.kegg.",i), heatmap = FALSE)
clust.kegg.p.esg.up <- esset.grp(clust.kegg.p$greater,
                                  TrainExprs,ref=c(1:ncol(TrainExprs))[-ind], gsets = kegg.gs,
                                  outname = paste0("clust.kegg.UP_",i),
                                  test4up = TRUE, output = TRUE,compare="as.group")
   assign(paste0("clust.up_",i), clust.kegg.p.esg.up)
   
   clust.kegg.p.esg.dn <- esset.grp(clust.kegg.p$less,
                                  TrainExprs,ref=c(1:ncol(TrainExprs))[-ind], gsets = kegg.gs,
                                  outname = paste0("clust.kegg.DN_",i),
                                  test4up = TRUE, output = TRUE,compare="as.group")
   assign(paste0("clust.dn_",i), clust.kegg.p.esg.dn)

   assign(paste0("g",i), clust.kegg.p.esg.up$essentialSets[1:3])
   assign(paste0("l",i), clust.kegg.p.esg.dn$essentialSets[1:3])
   
 }


#Ploting Gage

Paths= c(unlist(mget(paste0("g", 1:length(levels(CLASS))))), unlist(mget(paste0("l", 1:length(levels(CLASS))))))
Paths= unique(Paths)
Path_name<- substring(Paths,12)

BAR=list()
for( i in 1:length(levels(CLASS))){
bar= data.frame(stat=get(paste0("clust.kegg.p_",i))$stats[Paths,1],path=Path_name,subtype=i)
bar= bar[order(bar$stat),]
BAR[[i]]=bar
}

BAR=do.call(rbind,BAR)
BAR$path <- factor(BAR$path, levels=unique(BAR$path) )


# Basic barplot
p<-ggplot(data=BAR, aes(x=path, y=stat,fill=as.factor(subtype))) +
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values=c(1:length(levels(CLASS))))
p

# Horizontal bar plot
p + coord_flip()

## heatmap data prep ##


split_bar<- split(BAR, BAR$subtype)

heatmap_table= list()
for(i in 1:length(split_bar)){
 
st<- as.data.frame(split_bar[[i]])
st<- arrange(st, path)
 heatmap_table[[i]]= st$stat
  }
heatmap_table= do.call(cbind, heatmap_table)
colnames(heatmap_table) = paste0("subtype.",1:length(split_bar))
rownames(heatmap_table)<- as.character(st$path)

  
## heatmap ##
colfuncR <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))
col<-colfuncR(100)
heatmap.3(as.matrix(heatmap_table), col= col, trace="none",
           sepcolor="white",
           colsep=1:ncol(heatmap_table),
           rowsep=1:nrow(heatmap_table),
           dendrogram = "none",margin=c(6,24),cexRow=1,cexCol = 1)

##
###Comparison with random signatures
## Montecarlo randomization pval

CdistBreast=list()
CresBreast=list()
for(j in testSol){
iter=1000


k=RES[RES$solution==j,"k"]
ngenes=RES[RES$solution==j,"ngenes"]
W=matrix(0,nrow=iter,ncol=chrom_length)
for(it in 1:iter){
  ind=sample(1:chrom_length,ngenes,replace=F)
  W[it,ind]=1
}

Cdist=NULL
for(it in 1:iter){
R=as.logical(W[it,])
D=mydist(prob_matrix[R,])
  hsp_class=as.factor(galgoR::cluster_algorithm(D,k)$cluster)
  coxsimple=coxph(OS~hsp_class)
  pred= predict(coxsimple)
  cind= concordance.index(x=pred, surv.time=OS[,1], surv.event=OS[,2], method="noether",outx=FALSE, na.rm=TRUE)$c.index
  Cdist=c(Cdist,cind)
}

Name=paste("Breast",k,sep=".")
CdistBreast[[Name]]=Cdist
CresBreast[[Name]]= as.numeric(RES[RES$solution==j,"trainC"])

}

