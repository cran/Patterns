## ----setup, include = FALSE----------------------------------------------
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----loadPatterns, echo=FALSE, cache=FALSE, eval = LOCAL-----------------
require(Patterns)

## ----mergeselection, warning=FALSE, fig.keep='first', eval = LOCAL-------
selection<-Patterns::unionMicro(list(selection1,selection2,selection3,selection4))
summary(selection)

## ----addgroupselection, warning=FALSE, fig.keep='last', eval = LOCAL-----
selection@group <- rep(NA, length(selection@name))
names(selection@group) <- selection@name
selection@group[selection@name %in% selection4@name] <- 4
selection@group[selection@name %in% selection3@name] <- 3
selection@group[selection@name %in% selection2@name] <- 2
selection@group[selection@name %in% selection1@name] <- 1
plot(selection)

## ----plotF, eval = LOCAL-------------------------------------------------
plotF(network@F, choice='F')

## ----saveinference, eval=FALSE-------------------------------------------
#  save(list=c("selection"),file="selection.RData")
#  save(list=c("infos"),file="infos.RData")

## ----plotTFinsel, warning=FALSE, eval = LOCAL----------------------------
matplot(t(selection@microarray[tfs,]),type="l",lty=1)

## ----plotTFinsel2, warning=FALSE, eval = LOCAL---------------------------
kk<-kmeans((selection@microarray[tfs,]),10)
matplot(t(kk$centers),type="l",lty=1)

## ----TODO, warning=FALSE, echo=FALSE, eval=FALSE-------------------------
#  #TO DO
#  #Focus on TF that were not selected.
#  
#  indice<-which(CLL[,2] %in% TF[tfs<-which(! TF %in% selection@gene_ID)])
#  a<-1:200
#  matplot(log(t(agg_S[indice[a],]/agg_US[indice[a],])),lty=1,type="l")
#  kkk<-kmeans(log((agg_S[indice,]/agg_US[indice,])),10)
#  matplot(t(kkk$centers),type="l",lty=1)
#  
#  poi<-indice[which(kkk$cluster==2 )]
#  matmat<-log((agg_S[poi,]/agg_US[poi,]))
#  
#  addna<-function(mat,t,p){
#  
#  
#    mat2<-mat[,1:t]
#    for(i in 2:p){
#      print(1:t+(i-1)*t)
#      mat2<-cbind(mat2,rep(NA,nrow(mat2)),mat[,1:t+(i-1)*t])
#    }
#    return(mat2)
#  }
#  
#  pdf("forgotten_TF.pdf",width=15,height=5)
#  for(i in 1:15){
#  poi<-indice[which(kkk$cluster==i )]
#  if(length(poi)>2){
#  matmat<-log((agg_S[poi,]/agg_US[poi,]))
#  #matplot(t(matmat),lty=1,type="l")
#  matplot(t(addna(matmat,4,5)),lty=1,type="l")}
#  }
#  dev.off()
#  abline(v=c(2,6,10,14,18))
#  
#  
#  poi<-indice[which(kkk$cluster==1 )]
#  matplot(log(t(agg_S[poi,]/agg_US[poi,])),lty=1,type="l")
#  TFi<-function(x) length(which(TF %in% x))
#  
#  
#  
#  n<-40
#  kre<-kmeans(selection@microarray,n)
#  kre
#  lll<-split(selection@gene_ID,kre$cluster)
#  
#  require(DCGL)
#  require("clusterProfiler")
#  require("AnnotationFuncs")
#  require(org.Hs.eg.db)
#  
#  pp<-list()
#  
#  for(k in 1:2){
#    print(k)
#  pp[[k]]<-translate(lll[[k]],from=org.Hs.egSYMBOL2EG,simplify=TRUE)
#  # GOs[[k]]<-enrichGO(pp, organism = "human", ont = "MF", pvalueCutoff = 0.05,
#  #          pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 5,
#  #          readable = FALSE)
#  
#  }
#  
#  names(pp)<-paste("X",1:2,sep="")
#  test<-compareCluster(pp,fun="enrichGO", organism="human", pvalueCutoff=0.05)
#  
#  
#  plot(test)
#  
#  translate(lll[[k]],from=org.Hs.egSYMBOL2EG,simplify=TRUE)
#  
#  
#  TFu<-(unlist(lapply(pp,TFi)))
#  TFy<-unlist(lapply(pp,length))
#  
#  plot(TFu/TFy)
#  plot(TFu)
#  sum(TFu)
#  
#  entrez<-translate(selection@gene_ID,from=org.Hs.egSYMBOL2EG,simplify=TRUE)
#  
#  geneName<-translate(entrez[which(TF %in% entrez)],from=org.Hs.egSYMBOL,simplify=TRUE)
#  
#  which(selection@gene_ID %in% "EGR1")

