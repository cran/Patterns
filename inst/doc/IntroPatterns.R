## ----setup, include = FALSE----------------------------------------------
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5
)

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("Patterns")

## ---- eval = FALSE-------------------------------------------------------
#  devtools::install_github("fbertran/Patterns")

## ----microarrayclass, message=FALSE, warning=FALSE, eval = LOCAL---------
library(Patterns)
if(!require(CascadeData)){install.packages("CascadeData")}
data(micro_US)
micro_US<-as.micro_array(micro_US[1:100,],time=c(60,90,210,390),subject=6)
str(micro_US)

## ----plotmicroarrayclass, fig.keep='all', eval = LOCAL-------------------
summary(micro_US)
plot(micro_US)

## ----createF, eval = LOCAL-----------------------------------------------
Ti<-4
Ngrp<-4

Fmat=array(0,dim=c(Ti,Ti,Ngrp^2))

for(i in 1:(Ti^2)){
  if(((i-1) %% Ti) > (i-1) %/% Ti){
    Fmat[,,i][outer(1:Ti,1:Ti,function(x,y){0<(x-y) & (x-y)<2})]<-1
    }
}

## ----CascadeInit, eval = LOCAL-------------------------------------------
Fbis=Patterns::CascadeFinit(Ti,Ngrp,low.trig=FALSE)
str(Fbis)

## ----CascadeInitCheck, eval = LOCAL--------------------------------------
print(all(Fmat==Fbis))

## ----CascadeFinalize, eval = LOCAL---------------------------------------
Fmat[,,3]<-Fmat[,,3]*0.2
Fmat[3,1,3]<-1
Fmat[4,2,3]<-1
Fmat[,,4]<-Fmat[,,3]*0.3
Fmat[4,1,4]<-1
Fmat[,,8]<-Fmat[,,3]

## ----randomN, eval = LOCAL-----------------------------------------------
set.seed(1)
Net<-Patterns::network_random(
  nb=100,
  time_label=rep(1:4,each=25),
  exp=1,
  init=1,
  regul=round(rexp(100,1))+1,
  min_expr=0.1,
  max_expr=2,
  casc.level=0.4
)
Net@F<-Fmat
str(Net)

## ----plotnet1, eval = LOCAL----------------------------------------------
Patterns::plot(Net, choice="network")

## ----plotnet2, eval = LOCAL----------------------------------------------
plot(Net, choice="network", gr=rep(1:4,each=25))

## ----plotF, eval = LOCAL-------------------------------------------------
plot(Net, choice="F")

## ----plotFpixmap, eval = LOCAL-------------------------------------------
plot(Net, choice="Fpixmap")

## ----genesimul, message=FALSE, warning=FALSE, eval = LOCAL---------------
set.seed(1)
M <- Patterns::gene_expr_simulation(
  network=Net,
  time_label=rep(1:4,each=25),
  subject=5,
  peak_level=200,
  act_time_group=1:4)
str(M)

## ----Fresults, eval = LOCAL----------------------------------------------
plot(Net_inf_P, choice="F")

## ----heatresults, eval = LOCAL-------------------------------------------
stats::heatmap(Net_inf_P@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----Finitshow, eval = LOCAL---------------------------------------------
Ti<-4;
ngrp<-4
nF<-ngrp^2
Finit<-array(0,c(Ti,Ti,nF))	
              for(ii in 1:nF){    
                if((ii%%(ngrp+1))==1){
                  Finit[,,ii]<-0
                } else {
                  Finit[,,ii]<-cbind(rbind(rep(0,Ti-1),diag(1,Ti-1)),rep(0,Ti))+rbind(cbind(rep(0,Ti-1),diag(1,Ti-1)),rep(0,Ti))
                }
              }

## ----Fshapeshow, eval = LOCAL--------------------------------------------
Fshape<-array("0",c(Ti,Ti,nF)) 
for(ii in 1:nF){  
  if((ii%%(ngrp+1))==1){
    Fshape[,,ii]<-"0"
  } else {
    lchars <- paste("a",1:(2*Ti-1),sep="")
    tempFshape<-matrix("0",Ti,Ti)
    for(bb in (-Ti+1):(Ti-1)){
      tempFshape<-replaceUp(tempFshape,matrix(lchars[bb+Ti],Ti,Ti),-bb)
    }
    tempFshape <- replaceBand(tempFshape,matrix("0",Ti,Ti),0)
    Fshape[,,ii]<-tempFshape
  }
}

## ----Fshapeothershow, eval = LOCAL---------------------------------------
TestIndic=matrix(!((1:(Ti^2))%%(ngrp+1)==1),byrow=TRUE,ngrp,ngrp)
TestIndic

## ----Fshapeothershow2, eval = LOCAL--------------------------------------
IndicFinit(Ti,ngrp,TestIndic)
IndicFshape(Ti,ngrp,TestIndic)

## ----plotfshape1, eval = LOCAL-------------------------------------------
plotF(Fshape,choice="Fshape")

## ----plotfshape2, fig.keep="none", eval = LOCAL--------------------------
plotF(CascadeFshape(4,4),choice="Fshape")

## ----plotfshape3, fig.keep="none", eval = LOCAL--------------------------
plotF(IndicFshape(Ti,ngrp,TestIndic),choice="Fshape")

## ----FresultsLC, eval = LOCAL--------------------------------------------
plot(Net_inf_P_S, choice="F")

## ----heatresultsLC, fig.keep="none", eval = LOCAL------------------------
stats::heatmap(Net_inf_P_S@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----Fresultslasso2, fig.keep="none", eval = LOCAL-----------------------
plot(Net_inf_P_Lasso2, choice="F")

## ----heatresultslasso2, fig.keep="none", eval = LOCAL--------------------
stats::heatmap(Net_inf_P_Lasso2@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----netinfPriors, eval = LOCAL------------------------------------------
Weights_Net=slot(Net,"network")
Weights_Net[Net@network!=0]=.1        
Weights_Net[Net@network==0]=1000

## ----Fresultslasso2Weighted, fig.keep="none", eval = LOCAL---------------
plot(Net_inf_P_Lasso2_Weighted, choice="F")

## ----heatresultslasso2Weighted, fig.keep="none", eval = LOCAL------------
stats::heatmap(Net_inf_P_Lasso2_Weighted@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----FresultsSPLS, fig.keep="none", eval = LOCAL-------------------------
plot(Net_inf_P_SPLS, choice="F")

## ----heatresultsSPLS, fig.keep="none", eval = LOCAL----------------------
stats::heatmap(Net_inf_P_SPLS@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----FresultsEN, fig.keep="none", eval = LOCAL---------------------------
plot(Net_inf_P_ELASTICNET, choice="F")

## ----heatresultsEN, fig.keep="none", eval = LOCAL------------------------
stats::heatmap(Net_inf_P_ELASTICNET@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----FresultsStab, fig.keep="none", eval = LOCAL-------------------------
plot(Net_inf_P_stability, choice="F")

## ----heatresultsStab, fig.keep="none", eval = LOCAL----------------------
stats::heatmap(Net_inf_P_stability@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----FresultsStabWeight, fig.keep="none", eval = LOCAL-------------------
plot(Net_inf_P_StabWeight, choice="F")

## ----heatresultsStabWeight, fig.keep="none", eval = LOCAL----------------
stats::heatmap(Net_inf_P_StabWeight@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----FresultsRobust, fig.keep="none", eval = LOCAL-----------------------
plot(Net_inf_P_Robust, choice="F")

## ----heatresultsRobust, fig.keep="none", eval = LOCAL--------------------
stats::heatmap(Net_inf_P_Robust@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----FresultsSBpre, fig.keep="none", eval = LOCAL------------------------
detach("package:Patterns", unload=TRUE)
library(Patterns)

## ----FresultsSB, fig.keep="none", eval = LOCAL---------------------------
plot(Net_inf_P_SelectBoost, choice="F")

## ----heatresultsSB, fig.keep="none", eval = LOCAL------------------------
stats::heatmap(Net_inf_P_SelectBoost@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----FresultsSBWpre, fig.keep="none", eval = LOCAL-----------------------
detach("package:Patterns", unload=TRUE)
library(Patterns)

## ----FresultsSBW, fig.keep="none", eval = LOCAL--------------------------
plot(Net_inf_P_SelectBoostWeighted, choice="F")

## ----heatresultsSBW, fig.keep="none", eval = LOCAL-----------------------
stats::heatmap(Net_inf_P_SelectBoostWeighted@network, Rowv = NA, Colv = NA, scale="none", revC=TRUE)

## ----evolution, warning=FALSE, eval=FALSE--------------------------------
#  data(network)
#  sequence<-seq(0,0.2,length.out=20)
#  evolution(network,sequence,type.ani = "gif")
#  evolution(network,sequence,type.ani = "html")

## ----compare, warning=FALSE, eval = LOCAL--------------------------------
data(simul)

#Comparing true and inferred networks
Crit_values=NULL

#Here are the cutoff level tested
test.seq<-seq(0,max(abs(Net_inf_PL@network*0.9)),length.out=200)
for(u in test.seq){
	Crit_values<-rbind(Crit_values,Patterns::compare(Net,Net_inf_PL,u))
}
matplot(test.seq,Crit_values,type="l",ylab="Criterion value",xlab="Cutoff level",lwd=2)
legend(x="topleft", legend=colnames(Crit_values), lty=1:5,col=1:5,ncol=2,cex=.9)

## ----plotnet, warning=FALSE, eval = LOCAL--------------------------------
data(Selection)
plot(networkCascade,nv=0.133, gr=Selection@group)

## ----microselection, warning=FALSE, cache=FALSE, eval = LOCAL------------
library(Patterns)
library(CascadeData)
data(micro_S)
micro_S<-as.micro_array(micro_S,time=c(60,90,210,390),subject=6,gene_ID=rownames(micro_S))
data(micro_US)
micro_US<-as.micro_array(micro_US,time=c(60,90,210,390),subject=6,gene_ID=rownames(micro_US))

## ----microselection5, warning=FALSE, eval = LOCAL------------------------
Selection<-unionMicro(Selection1,Selection2)
Selection<-unionMicro(Selection,Selection3)
Selection<-unionMicro(Selection,Selection4)
head(Selection)

## ----microselection6, warning=FALSE, fig.keep="none", eval = LOCAL-------
summary(Selection)

## ----microselection7, warning=FALSE, fig.keep="none", eval = LOCAL-------
plot(Selection)

