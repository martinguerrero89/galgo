load("./Data/RNA_OV_all.rda")
downl=names(esets)


#Run galgo in ovarian

population= 300                   #Number of individuals to evaluate
generations=500                 #Number of generations
nCV=5                             #Number of crossvalidations for function "crossvalidation"
GPU= TRUE                         # to use gpuR
TournamentSize=2
period=1825

trainSet=esets[["TCGA.RNAseq"]] ;h="HGSOC"
prob_matrix= exprs(trainSet)
clinical=pData(trainSet)
OS=Surv(time=clinical$time,event=clinical$status)
chrom_length= nrow(prob_matrix)   #length of chromosome


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


SSP= c("Verhaak.subtypes","Bentink.subtypes","Helland.subtypes","Konecny.subtypes","Consensus.subtypes")
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
  
  
  
riskP=c(paste0(c(galgo,SSP),".Pred")) 
resMatrix <- as.list(NULL)
for(i in downl){
res=pData(esets[[i]])
SURV=Surv(res[,"time"],res[,"status"])
#In=!is.na(SURV)

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
for(i in 3:length(downl[-which(downl=="TCGA.RNAseq")])){
  Comb=combineTwoExpressionSet(Comb,esets[[i]])
}



library(survminer)
  tumortotal1 <- survfit(Surv(Comb$time,Comb$status)~ Comb$Consensus.subtypes)
  tumortotal2 <- survfit(Surv(Comb$time,Comb$status)~ Comb$"4_result.9")
  tumortotal1diff <- survdiff(Surv(Comb$time,Comb$status)~ Comb$Consensus.subtypes)
  tumortotal2diff <- survdiff(Surv(Comb$time,Comb$status)~ Comb$"4_result.9")
  tumortotal1pval<- pchisq(tumortotal1diff$chisq, length(tumortotal1diff$n) - 1, lower.tail = FALSE) 
  tumortotal2pval<- pchisq(tumortotal2diff$chisq, length(tumortotal2diff$n) - 1, lower.tail = FALSE) 
  
  
  SURV=list(Cons=tumortotal1, GALGO=tumortotal2)
  COLS=c(1:8,10)
  par(cex=1.35, mar=c(3.8, 3.8, 2.5, 2.5) + 0.1)
  p=ggsurvplot(SURV,combine=TRUE,data=Comb,risk.table=TRUE,pval=TRUE,palette="dark2", title="Combined Set \n Ovarian survival comparison", surv.scale="percent", conf.int=FALSE, xlab="time (days)", ylab="survival(%)", xlim=c(0,period),break.time.by = 365, ggtheme = theme_minimal(), risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,censor=FALSE)
  print(p)


  tumortotal1 <- survfit(Surv(esets[["TCGA.RNAseq"]]$time,esets[["TCGA.RNAseq"]]$status)~ esets[["TCGA.RNAseq"]]$Consensus.subtypes)
  tumortotal2 <- survfit(Surv(esets[["TCGA.RNAseq"]]$time,esets[["TCGA.RNAseq"]]$status)~ esets[["TCGA.RNAseq"]]$"4_result.9")
  tumortotal1diff <- survdiff(Surv(esets[["TCGA.RNAseq"]]$time,esets[["TCGA.RNAseq"]]$status)~ esets[["TCGA.RNAseq"]]$Consensus.subtypes)
  tumortotal2diff <- survdiff(Surv(esets[["TCGA.RNAseq"]]$time,esets[["TCGA.RNAseq"]]$status)~ esets[["TCGA.RNAseq"]]$"4_result.9")
  tumortotal1pval<- pchisq(tumortotal1diff$chisq, length(tumortotal1diff$n) - 1, lower.tail = FALSE) 
  tumortotal2pval<- pchisq(tumortotal2diff$chisq, length(tumortotal2diff$n) - 1, lower.tail = FALSE) 
  
  
  
  SURV=list(Cons=tumortotal1, GALGO=tumortotal2)
  par(cex=1.35, mar=c(3.8, 3.8, 2.5, 2.5) + 0.1)
  p=ggsurvplot(SURV,combine=TRUE,data=Comb,risk.table=TRUE,pval=TRUE,palette="dark2", title="Train Set \n Ovarian survival comparison", surv.scale="percent", conf.int=FALSE, xlab="time (days)", ylab="survival(%)", xlim=c(0,period),break.time.by = 365, ggtheme = theme_minimal(), risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,censor=FALSE)
  print(p)


finalSig="4_result.9"
V=unlist(strsplit(finalSig,"_"))
R=as.logical(X1[V[2],1:chrom_length])
  k=as.numeric(V[1])
  
trainset="TCGA.RNAseq"
TrainExprs=exprs(esets[[trainset]])
rownames(TrainExprs)= fData(esets[[trainset]])$gene  
TrainClass= pData(esets[[trainset]])[,riskP]

library(RColorBrewer)
colfuncR <- colorRampPalette(rev(brewer.pal(11 , "Spectral" )))

Consensus= colfuncR(4)[as.numeric(as.factor(TrainClass[,riskP[11]]))]
Konecny= colfuncR(4)[as.numeric(as.factor(TrainClass[,riskP[10]]))]
Helland= colfuncR(4)[as.numeric(as.factor(TrainClass[,riskP[9]]))]
Bentink= colfuncR(2)[as.numeric(as.factor(TrainClass[,riskP[8]]))]
Verhaak= colfuncR(4)[as.numeric(as.factor(TrainClass[,riskP[7]]))]
galgo2=colfuncR(2)[as.numeric(as.factor(TrainClass[,riskP[6]]))]
galgo3=colfuncR(3)[as.numeric(as.factor(TrainClass[,riskP[5]]))]
galgo4=colfuncR(4)[as.numeric(as.factor(TrainClass[,riskP[4]]))]
galgo6=colfuncR(6)[as.numeric(as.factor(TrainClass[,riskP[3]]))]
galgo5=colfuncR(5)[as.numeric(as.factor(TrainClass[,riskP[2]]))]
galgo10=colfuncR(10)[as.numeric(as.factor(TrainClass[,riskP[1]]))]


COLS=data.frame(galgo4, galgo2,galgo3,galgo5,galgo6, galgo10,Verhaak, Bentink, Helland,Konecny,Consensus)

CLASS=TrainClass[,paste(finalSig,"Pred",sep=".")]
Ord= order(CLASS)
library(gplots)

colfunc <- colorRampPalette(c("blue", "white", "red"))
breaks=seq(-2, 2, by=0.1) 
#now add outliers
breaks=append(breaks, 5)
breaks=append(breaks, -5, 0)

HM1=heatmap.3(TrainExprs[R,Ord],Colv= FALSE,hclustfun=hclust,dendrogram="none", distfun=mydist2,trace="none",scale="row",col=colfunc,breaks=breaks,ColSideColors =as.matrix(COLS[Ord,]),ColSideColorsSize=3)


#Testset (Comb)

TestExprs=ComBat(exprs(Comb),pData(Comb)$batch)
rownames(TestExprs)= fData(Comb)$gene  
TestClass= pData(Comb)[,riskP]


Consensus= colfuncR(4)[as.numeric(as.factor(TestClass[,riskP[11]]))]
Konecny= colfuncR(4)[as.numeric(as.factor(TestClass[,riskP[10]]))]
Helland= colfuncR(4)[as.numeric(as.factor(TestClass[,riskP[9]]))]
Bentink= colfuncR(2)[as.numeric(as.factor(TestClass[,riskP[8]]))]
Verhaak= colfuncR(4)[as.numeric(as.factor(TestClass[,riskP[7]]))]
galgo2=colfuncR(2)[as.numeric(as.factor(TestClass[,riskP[6]]))]
galgo3=colfuncR(3)[as.numeric(as.factor(TestClass[,riskP[5]]))]
galgo4=colfuncR(4)[as.numeric(as.factor(TestClass[,riskP[4]]))]
galgo6=colfuncR(6)[as.numeric(as.factor(TestClass[,riskP[3]]))]
galgo5=colfuncR(5)[as.numeric(as.factor(TestClass[,riskP[2]]))]
galgo10=colfuncR(10)[as.numeric(as.factor(TestClass[,riskP[1]]))]


COLS=data.frame(galgo4, galgo2,galgo3,galgo5,galgo6, galgo10,Verhaak, Bentink, Helland,Konecny,Consensus)

CLASS=TestClass[,paste(finalSig,"Pred",sep=".")]
Ord= order(CLASS)

HM2=heatmap.3(TestExprs[R,Ord][rev(HM1$rowInd),],Colv= FALSE,Rowv=FALSE,hclustfun=hclust,dendrogram="none", distfun=mydist2,trace="none",scale="row",col=colfunc,breaks=breaks,ColSideColors =as.matrix(COLS[Ord,]),ColSideColorsSize=3)

library(gage)
library(gageData)
data(go.sets.hs)
kegg.gs=go.sets.hs

CLASS=as.factor(TrainClass[,paste(finalSig,"Pred",sep=".")])
TrainExprs= exprs(esets[[trainset]])


ind1<- which( CLASS == levels(CLASS)[1])
ind2<- which( CLASS == levels(CLASS)[2])
ind3<- which( CLASS == levels(CLASS)[3])
ind4<- which( CLASS == levels(CLASS)[4])


clust1.kegg.p <- gage(TrainExprs,ref=c(ind2,ind3,ind4), gsets = kegg.gs,compare="as.group")
clust2.kegg.p <- gage(TrainExprs,ref=c(ind1,ind3,ind4), gsets = kegg.gs,compare="as.group")
clust3.kegg.p <- gage(TrainExprs,ref=c(ind1,ind2,ind4), gsets = kegg.gs,compare="as.group")
clust4.kegg.p <- gage(TrainExprs,ref=c(ind1,ind2,ind3), gsets = kegg.gs,compare="as.group")

##Sort and count signficant gene sets based on q- or p-value cutoffs:
##plot heatmap by clust

clust1.kegg.p.sig<- sigGeneSet(clust1.kegg.p, outname = "clust1.kegg", heatmap = FALSE)
clust2.kegg.p.sig<- sigGeneSet(clust2.kegg.p, outname = "clust2.kegg", heatmap = FALSE)
clust3.kegg.p.sig<- sigGeneSet(clust3.kegg.p, outname = "clust3.kegg", heatmap = FALSE)
clust4.kegg.p.sig<- sigGeneSet(clust4.kegg.p, outname = "clust4.kegg", heatmap = FALSE)

## extract a non-redundant signcant gene set list

# subtype I #
clust1.kegg.p.esg.up <- esset.grp(clust1.kegg.p$greater,
                                  TrainExprs,ref=c(ind2,ind3,ind4), gsets = kegg.gs,
                                  outname = "clust1.kegg.up",
                                  test4up = TRUE, output = TRUE,compare="as.group")

clust1.kegg.p.esg.dn <- esset.grp(clust1.kegg.p$less,
                                  TrainExprs,ref=c(ind2,ind3,ind4), gsets = kegg.gs, 
                                  test4up = FALSE, output = TRUE, 
                                  outname = "clust1.kegg.dn",compare="as.group")

g1<- clust1.kegg.p.esg.up$essentialSets[1:3]
l1<- clust1.kegg.p.esg.dn$essentialSets[1:3]

# subtype II #

clust2.kegg.p.esg.up <- esset.grp(clust2.kegg.p$greater,
                                 TrainExprs,ref=c(ind1,ind3,ind4), gsets = kegg.gs,
                                  outname = "clust2.kegg.up",
                                  test4up = TRUE, output = TRUE,compare="as.group")

clust2.kegg.p.esg.dn <- esset.grp(clust2.kegg.p$less,
                                  TrainExprs,ref=c(ind1,ind3,ind4), gsets = kegg.gs, 
                                  test4up = FALSE, output = TRUE, 
                                  outname = "clust2.kegg.dn",compare="as.group")

g2<- clust2.kegg.p.esg.up$essentialSets[1:3]
l2<- clust2.kegg.p.esg.dn$essentialSets[1:3]

# subtype III #

clust3.kegg.p.esg.up <- esset.grp(clust3.kegg.p$greater,
                                 TrainExprs,ref=c(ind1,ind2,ind4), gsets = kegg.gs,
                                  outname = "clust3.kegg.up",
                                  test4up = TRUE, output = TRUE,compare="as.group")

clust3.kegg.p.esg.dn <- esset.grp(clust3.kegg.p$less,
                                  TrainExprs,ref=c(ind1,ind2,ind4), gsets = kegg.gs, 
                                  test4up = FALSE, output = TRUE, 
                                  outname = "clust3.kegg.dn",compare="as.group")

g3<- clust3.kegg.p.esg.up$essentialSets[1:3]
l3<- clust3.kegg.p.esg.dn$essentialSets[1:3]


# subtype IV #

clust4.kegg.p.esg.up <- esset.grp(clust4.kegg.p$greater,
                                 TrainExprs,ref=c(ind1,ind2,ind3), gsets = kegg.gs,
                                  outname = "clust4.kegg.up",
                                  test4up = TRUE, output = TRUE,compare="as.group")

clust4.kegg.p.esg.dn <- esset.grp(clust4.kegg.p$less,
                                  TrainExprs,ref=c(ind1,ind2,ind3), gsets = kegg.gs, 
                                  test4up = FALSE, output = TRUE, 
                                  outname = "clust4.kegg.dn",compare="as.group")

g4<- clust4.kegg.p.esg.up$essentialSets[1:3]
l4<- clust4.kegg.p.esg.dn$essentialSets[1:3]



Paths<- c(g1,g2,g3,g4,l1,l2,l3,l4)
Paths<- unique(c(g1,g2,g3,g4,l1,l2,l3,l4))
Path_name<- substring(Paths,12)

bar1=data.frame(stat=clust1.kegg.p$stats[Paths,1],path=Path_name,subtype=1)
bar1= bar1[order(bar1$stat),]
bar2=data.frame(stat=clust2.kegg.p$stats[Paths,1],path=Path_name,subtype=2)
bar2= bar2[order(bar2$stat),]
bar3=data.frame(stat=clust3.kegg.p$stats[Paths,1],path=Path_name,subtype=3)
bar3= bar3[order(bar3$stat),]
bar4=data.frame(stat=clust4.kegg.p$stats[Paths,1],path=Path_name,subtype=4)
bar4= bar4[order(bar4$stat),]


bar=rbind(bar1,bar2,bar3,bar4)

bar$path <- factor(bar$path, levels=unique(bar$path) )


## heatmap data prep ##
library(dplyr)

split_bar<- split(bar, bar$subtype)

st_I<- as.data.frame(split_bar[1])
st_I<- arrange(st_I, X1.path)

st_II<- as.data.frame(split_bar[2])
st_II<- arrange(st_II, X2.path)


st_III<- as.data.frame(split_bar[3])
st_III<- arrange(st_III, X3.path)


st_IV<- as.data.frame(split_bar[4])
st_IV<- arrange(st_IV, X4.path)


heatmap_table<- as.data.frame(cbind(st_I$X1.stat, st_II$X2.stat,st_III$X3.stat,st_IV$X4.stat))
colnames(heatmap_table)<- c("Subtype I", "Subtype II","Subtype III","Subtype IV")
rownames(heatmap_table)<- as.character(st_I$X1.path)

## heatmap ##
library(gplots)
colfuncR <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))
col<-colfuncR(100)


 heatmap.3(as.matrix(heatmap_table), col= col, trace="none",
           sepcolor="white",
           colsep=1:ncol(heatmap_table),
           rowsep=1:nrow(heatmap_table),
           dendrogram = "none",margin=c(6,28),cexRow=1,cexCol = 1)





library(xlsx)
for(i in names(CentroidsList)){
write.xlsx(CentroidsList[[i]], file="Ovarian_Galgo_Centroids.xlsx", sheetName=i,append=TRUE, row.names=TRUE)
}


write.xlsx(RR_data[,c(6,5,1:4)], file="Ovarian_CI_Meta.xlsx", sheetName="CI",append=TRUE, row.names=FALSE)
write.xlsx(ccmData, file="Ovarian_CI_Meta.xlsx", sheetName="Comparison",append=TRUE, row.names=TRUE)


#Cindex random distribution

CdistOvarian=list()
CresOvarian=list()
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


Name=paste("Ovarian",j,sep=".")
CdistOvarian[[Name]]=Cdist
CresOvarian[[Name]]= as.numeric(RES[RES$solution==j,"trainC"])

}

