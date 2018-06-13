
rm(list=ls())

library(dplyr)
library(SuperLearner)
library(tmle)
library(washb)
library(caret)
library(zoo)

source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")

setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned Data")
load("analysis_datasets.Rdata")



##########################
#Unadjusted analysis
##########################

#Examine mean rain exposure in heavy rain and non heavy rain groups

HR.rainmean7<-d %>% group_by(HeavyRain.lag7) %>% summarize(mean.rainfall=mean(rain.ave7.lag7, na.rm=T)) %>% as.data.frame()
HR.rainmean14<-d %>% group_by(HeavyRain.lag14) %>% summarize(mean.rainfall=mean(rain.ave7.lag14, na.rm=T)) %>% as.data.frame()
HR.rainmean21<-d %>% group_by(HeavyRain.lag21) %>% summarize(mean.rainfall=mean(rain.ave7.lag21, na.rm=T)) %>% as.data.frame()

save(HR.rainmean7, 
     HR.rainmean14,
     HR.rainmean21,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/HRmean.Results.Rdata")


Hrain.unadj<-matrix(NA, nrow=3, ncol=11)
Hrain.unadj[1,]<-HeavyRainTMLE(dat=d,i=7)
Hrain.unadj[2,]<-HeavyRainTMLE(dat=d,i=14)
Hrain.unadj[3,]<-HeavyRainTMLE(dat=d,i=21)
colnames(Hrain.unadj)<-c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.unadj



#heavy rain - stratified
Hrain.unadjLT3<-Hrain.unadjLT2<-Hrain.unadjLT1<-matrix(NA, nrow=3, ncol=11)


HR.strat.rainmean14<- d %>% group_by(LT8_T, HeavyRain.lag14) %>% summarize(mean.rainfall=mean(rain.ave7.lag14, na.rm=T), mean.h2s=mean(H2S, na.rm=T)) %>% as.data.frame()
  Hrain.unadjLT1[1,]<-HeavyRainTMLE(dat=d[d$LT8_T=="T1",],i=7)
  Hrain.unadjLT2[1,]<-HeavyRainTMLE(dat=d[d$LT8_T=="T2",],i=7)
  Hrain.unadjLT3[1,]<-HeavyRainTMLE(dat=d[d$LT8_T=="T3",],i=7)
  

HR.strat.rainmean14<- d %>% group_by(LT15_T, HeavyRain.lag14) %>% summarize(mean.rainfall=mean(rain.ave7.lag14, na.rm=T), mean.h2s=mean(H2S, na.rm=T)) %>% as.data.frame()
  Hrain.unadjLT1[2,]<-HeavyRainTMLE(dat=d[d$LT15_T=="T1",],i=14)
  Hrain.unadjLT2[2,]<-HeavyRainTMLE(dat=d[d$LT15_T=="T2",],i=14)
  Hrain.unadjLT3[2,]<-HeavyRainTMLE(dat=d[d$LT15_T=="T3",],i=14)

HR.strat.rainmean21<- d %>% group_by(LT22_T, HeavyRain.lag21) %>% summarize(mean.rainfall=mean(rain.ave7.lag21, na.rm=T), mean.h2s=mean(H2S, na.rm=T)) %>% as.data.frame()
  Hrain.unadjLT1[3,]<-HeavyRainTMLE(dat=d[d$LT22_T=="T1" ,],i=21)
  Hrain.unadjLT2[3,]<-HeavyRainTMLE(dat=d[d$LT22_T=="T2" ,],i=21)
  Hrain.unadjLT3[3,]<-HeavyRainTMLE(dat=d[d$LT22_T=="T3" ,],i=21)

colnames(Hrain.unadjLT3)<-colnames(Hrain.unadjLT2)<-colnames(Hrain.unadjLT1)<-c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.unadjLT1
Hrain.unadjLT2
Hrain.unadjLT3





##############
#Effect of heavy rain on 
#hydrogen sulfide in the 
#stored drinking water
##############

d$h2s <- d$H2S

H2S.Hrain.unadj<-matrix(NA, nrow=3, ncol=11)
H2S.Hrain.unadj[1,]<-H2STMLE(dat=d,i=1)
H2S.Hrain.unadj[2,]<-H2STMLE(dat=d,i=7)
H2S.Hrain.unadj[3,]<-H2STMLE(dat=d,i=14)
colnames(H2S.Hrain.unadj)<-c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
H2S.Hrain.unadj


H2s.unadjLT3<-H2s.unadjLT2<-H2s.unadjLT1<-matrix(NA, nrow=3, ncol=11)

  H2s.unadjLT1[1,]<-H2STMLE(dat=d[d$LT1_T=="T1",],i=1)
  H2s.unadjLT2[1,]<-H2STMLE(dat=d[d$LT1_T=="T2",],i=1)
  H2s.unadjLT3[1,]<-H2STMLE(dat=d[d$LT1_T=="T3",],i=1)
  
  H2s.unadjLT1[2,]<-H2STMLE(dat=d[d$LT8_T=="T1",],i=7)
  H2s.unadjLT2[2,]<-H2STMLE(dat=d[d$LT8_T=="T2",],i=7)
  H2s.unadjLT3[2,]<-H2STMLE(dat=d[d$LT8_T=="T3",],i=7)

  H2s.unadjLT1[3,]<-H2STMLE(dat=d[d$LT15_T=="T1",],i=14)
  H2s.unadjLT2[3,]<-H2STMLE(dat=d[d$LT15_T=="T2",],i=14)
  H2s.unadjLT3[3,]<-H2STMLE(dat=d[d$LT15_T=="T3",],i=14)

colnames(H2s.unadjLT3)<-colnames(H2s.unadjLT2)<-colnames(H2s.unadjLT1)<-c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
H2s.unadjLT1
H2s.unadjLT2
H2s.unadjLT3







#save unadjusted objects
save(
Hrain.unadj,
Hrain.unadjLT1,
Hrain.unadjLT2,
Hrain.unadjLT3,
H2S.Hrain.unadj,
H2s.unadjLT1,
H2s.unadjLT2,
H2s.unadjLT3,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/rain.Results.unadjusted.Rdata")

#Save dataset
save(d, file="C:/Users/andre/Dropbox/Trichy analysis/Results/rain.datasets.Rdata")







####################
#Adjusted TMLE
####################


library=c("SL.mean","SL.glm","SL.glmnet", "SL.gam", "SL.bayesglm")


d$intdate<-as.numeric(d$intdate)


HR.adj<-matrix(NA, nrow<-3, ncol=8)
  set.seed(12345)
  HR.adj[1,]<-HeavyRainTMLE(d=d,Wvars=Wvars,i=7,sl.library=library)
  HR.adj[2,]<-HeavyRainTMLE(d=d,Wvars=Wvars,i=14,sl.library=library)
  HR.adj[3,]<-HeavyRainTMLE(d=d,Wvars=Wvars,i=21,sl.library=library)
colnames(HR.adj)<-c("a","b","c","d","PR","ci.lb","ci.ub","pvalue")
  rownames(HR.adj)<-c("1week.lag", "2week.lag", "3week.lag")
Hrain.adj<-HR.adj
Hrain.adj



H2S.Hrain.adj<-matrix(NA, nrow=3, ncol=8)
H2S.Hrain.adj[1,]<-H2STMLE(dat=d[!is.na(d$h2s),],Wvars=Wvars,i=1, sl.library=library)
H2S.Hrain.adj[2,]<-H2STMLE(dat=d[!is.na(d$h2s),],Wvars=Wvars,i=7, sl.library=library)
H2S.Hrain.adj[3,]<-H2STMLE(dat=d[!is.na(d$h2s),],Wvars=Wvars,i=14, sl.library=library)
colnames(H2S.Hrain.adj)<-c("a","b","c","d","PR", "ci.lb", "ci.ub", "pval")
H2S.Hrain.adj




####################
#Adjusted subgroup analysis
####################
#heavy rain
Hrain.adjLT3<-Hrain.adjLT2<-Hrain.adjLT1<-matrix(NA, nrow=3, ncol=8)
  LT<-d$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  d$LTQ<-2
  d$LTQ[LT < ntiles[1]]<-1
  d$LTQ[LT >= ntiles[2]]<-3
  Hrain.adjLT1[1,]<-HeavyRainTMLE(dat=d[d$LTQ==1,],i=7,Wvars=Wvars,sl.library=library)
  Hrain.adjLT2[1,]<-HeavyRainTMLE(dat=d[d$LTQ==2,],i=7,Wvars=Wvars,sl.library=library)
  Hrain.adjLT3[1,]<-HeavyRainTMLE(dat=d[d$LTQ==3,],i=7,Wvars=Wvars,sl.library=library)
  
  LT<-d$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  d$LTQ<-2
  d$LTQ[LT < ntiles[1]]<-1
  d$LTQ[LT >= ntiles[2]]<-3
  Hrain.adjLT1[2,]<-HeavyRainTMLE(dat=d[d$LTQ==1,],i=14,Wvars=Wvars,sl.library=library)
  Hrain.adjLT2[2,]<-HeavyRainTMLE(dat=d[d$LTQ==2,],i=14,Wvars=Wvars,sl.library=library)
  Hrain.adjLT3[2,]<-HeavyRainTMLE(dat=d[d$LTQ==3,],i=14,Wvars=Wvars,sl.library=library)
  
  LT<-d$LT22
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  d$LTQ<-2
  d$LTQ[LT < ntiles[1]]<-1
  d$LTQ[LT >= ntiles[2]]<-3
  Hrain.adjLT1[3,]<-HeavyRainTMLE(dat=d[d$LTQ==1,],i=21,Wvars=Wvars,sl.library=library)
  Hrain.adjLT2[3,]<-HeavyRainTMLE(dat=d[d$LTQ==2,],i=21,Wvars=Wvars,sl.library=library)
  Hrain.adjLT3[3,]<-HeavyRainTMLE(dat=d[d$LTQ==3,],i=21,Wvars=Wvars,sl.library=library)

colnames(Hrain.adjLT3)<-colnames(Hrain.adjLT2)<-colnames(Hrain.adjLT1)<-c("a","b","c","d","PR", "ci.lb", "ci.ub",   "p")
Hrain.adjLT1
Hrain.adjLT2
Hrain.adjLT3


#h2s
H2S.adjLT3<-H2S.adjLT2<-H2S.adjLT1<-matrix(NA, nrow=3, ncol=8)
  LT<-d$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  d$LTQ<-2
  d$LTQ[LT < ntiles[1]]<-1
  d$LTQ[LT >= ntiles[2]]<-3
  H2S.adjLT1[1,]<-H2STMLE(dat=d[d$LTQ==1,],i=1,Wvars=Wvars,sl.library=library)
  H2S.adjLT2[1,]<-H2STMLE(dat=d[d$LTQ==2,],i=1,Wvars=Wvars,sl.library=library)
  H2S.adjLT3[1,]<-H2STMLE(dat=d[d$LTQ==3,],i=1,Wvars=Wvars,sl.library=library)
  
  LT<-d$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  d$LTQ<-2
  d$LTQ[LT < ntiles[1]]<-1
  d$LTQ[LT >= ntiles[2]]<-3
  H2S.adjLT1[2,]<-H2STMLE(dat=d[d$LTQ==1,],i=7,Wvars=Wvars,sl.library=library)
  H2S.adjLT2[2,]<-H2STMLE(dat=d[d$LTQ==2,],i=7,Wvars=Wvars,sl.library=library)
  H2S.adjLT3[2,]<-H2STMLE(dat=d[d$LTQ==3,],i=7,Wvars=Wvars,sl.library=library)
  
  LT<-d$LT22
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  d$LTQ<-2
  d$LTQ[LT < ntiles[1]]<-1
  d$LTQ[LT >= ntiles[2]]<-3
  H2S.adjLT1[3,]<-H2STMLE(dat=d[d$LTQ==1,],i=14,Wvars=Wvars,sl.library=library)
  H2S.adjLT2[3,]<-H2STMLE(dat=d[d$LTQ==2,],i=14,Wvars=Wvars,sl.library=library)
  H2S.adjLT3[3,]<-H2STMLE(dat=d[d$LTQ==3,],i=14,Wvars=Wvars,sl.library=library)

colnames(H2S.adjLT3)<-colnames(H2S.adjLT2)<-colnames(H2S.adjLT1)<-c("a","b","c","d","PR", "ci.lb", "ci.ub",   "p")
H2S.adjLT1
H2S.adjLT2
H2S.adjLT3


save(
Hrain.adj,
Hrain.adjLT1,
Hrain.adjLT2,
Hrain.adjLT3,
H2S.Hrain.adj,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/rain.Results.adjusted.Rdata")

save(
  H2S.adjLT1,
H2S.adjLT2,
H2S.adjLT3,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/H2S.Results.adjusted.Rdata")








