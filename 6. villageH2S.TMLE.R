
rm(list=ls())

library(tidyverse)
library(SuperLearner)
library(tmle)
library(washb)
library(caret)
library(foreign)


setwd("C:/Users/andre/Dropbox/Trichy analysis/Data")
df<-read.dta("trichy_fu_vilwat.dta")
load("Cleaned Data/LaggedWeather.Rdata")



###################
#Weather data processing
###################
head(weather)

rainvars<-sapply(1:100, function(x) paste0("rain.",x))
rain<-subset(weather, select=c("year","month","day",rainvars))
#Change "trace" to 0.01 in rainfall variables
for(i in 1:100){
  levels(rain[,3+i])[which(levels(rain[,3+i])=="Trace")]<-"0.1"
}

#change rain vars to numeric:
rain<-apply(rain, 2, function(x) as.numeric(as.character(x))) 


ave.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) mean(x))) #NOTE: This still computes averages when there are NA's
    return(windowmean)
}

averain<-matrix(data=NA,nrow=nrow(rain),ncol=21)
ave.names<-rep(NA, 21)
for(i in 1:21){
  averain[,i]<-ave.window(i, 7, "rain", rain)
  ave.names[i]<- paste0("rain.ave7.lag",i)
}
averain<-cbind(rain[,1:3],averain)
colnames(averain)[4:ncol(averain)]<-ave.names
averain<-as.data.frame(averain)

#Create long term rain average
averain$LT1<-ave.window(8, 60, "rain", rain) #set i to 7 days after the associated 7-day lagged variable
averain$LT8<-ave.window(15, 60, "rain", rain)
averain$LT15<-ave.window(22, 60, "rain", rain)
averain$LT22<-ave.window(29, 60, "rain", rain)


#Extract means of long-term quartiles for the stratified tables
    LT<-averain$LT1
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  averain$LTQ<-2
  averain$LTQ[LT < ntiles[1]]<-1
  averain$LTQ[LT >= ntiles[2]]<-3
  averain %>% group_by(LTQ) %>% summarize(meanLT1=mean(LT1, na.rm=T))
  
    LT<-averain$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  averain$LTQ<-2
  averain$LTQ[LT < ntiles[1]]<-1
  averain$LTQ[LT >= ntiles[2]]<-3
  averain %>% group_by(LTQ) %>% summarize(meanLT8=mean(LT8, na.rm=T))
  
    LT<-averain$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  averain$LTQ<-2
  averain$LTQ[LT < ntiles[1]]<-1
  averain$LTQ[LT >= ntiles[2]]<-3
  averain %>% group_by(LTQ) %>% summarize(meanLT15=mean(LT15, na.rm=T))
  
    LT<-averain$LT22
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  averain$LTQ<-2
  averain$LTQ[LT < ntiles[1]]<-1
  averain$LTQ[LT >= ntiles[2]]<-3
  averain %>% group_by(LTQ) %>% summarize(meanLT22=mean(LT22, na.rm=T))
  


#Set date to date format for merge
  
averain$intdate<-as.Date(paste(averain$month,averain$day,averain$year, sep="/"),"%m/%d/%Y")
averain<-subset(averain, select=-c(year,month,day))


#Calculate heavy rainfall events
#90th percentile of rainfall
HeavyRainThres<-as.numeric(quantile(rain[rain[,4]!=0,4],probs = seq(0, 1, 0.1) ,na.rm=T)[9])


HeavyRain<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres,1,0))
table(HeavyRain)

 

HeavyRain.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) ifelse(sum(x)>0,1,0)))
  return(windowmean)
}


#PriorHeavyRain<-matrix(data=NA,nrow=nrow(rain),ncol=15)
PriorHeavyRain<-matrix(data=NA,nrow=nrow(rain),ncol=21)
PriorHeavyRain.names<-rep(NA, 21)
for(i in 1:21){
  PriorHeavyRain[,i]<-HeavyRain.window(i, 7, "rain", HeavyRain)
  PriorHeavyRain.names[i]<- paste0("HeavyRain.lag",i)
}
PriorHeavyRain<-cbind(rain[,1:3],PriorHeavyRain)
colnames(PriorHeavyRain)[4:ncol(PriorHeavyRain)]<-PriorHeavyRain.names
PriorHeavyRain<-as.data.frame(PriorHeavyRain)
head(PriorHeavyRain,30)



###################
#village Survey data processing
###################


#clean and select needed variables
df <- df %>% select(sid, svy, vdatecol, vdatearr, v10) %>% rename(h2s=v10, id=sid) %>% filter(!is.na(h2s))
df$Y <- df$h2s
head(df)

#fix errors in date input
df$vdatecol[31] <- "03/06/08"
df$vdatecol[87] <- "20/08/08"

df$date <- as.Date(df$vdatecol,format='%d/%m/%y')
df$year <- as.numeric(year(df$date))
df$month <- as.numeric(month(df$date))
df$day <- as.numeric(days(df$date))

averain <- averain %>% rename(date=intdate)

#merge weather with survey
d <- left_join(df, averain, by=c("date"))




##########################
#Unadjusted analysis
##########################

colnames(d)
#Mean h2s prevalence by quartile:
d_rain<-subset(d, rain.ave7.lag7>0)
low_hi_cut<-quantile(d_rain$rain.ave7.lag7, prob = seq(0, 1, length = 4), type = 5, na.rm=TRUE)[2:3]

d$raincat.7<-d$raincat.14<-d$raincat.21<-"No rain"
d$raincat.7[d$rain.ave7.lag7>0]<-"Low rain"
d$raincat.14[d$rain.ave7.lag14>0]<-"Low rain"
d$raincat.21[d$rain.ave7.lag21>0]<-"Low rain"
d$raincat.7[d$rain.ave7.lag7>low_hi_cut[1]]<-"Med rain"
d$raincat.14[d$rain.ave7.lag14>low_hi_cut[1]]<-"Med rain"
d$raincat.21[d$rain.ave7.lag21>low_hi_cut[1]]<-"Med rain"
d$raincat.7[d$rain.ave7.lag7>low_hi_cut[2]]<-"High rain"
d$raincat.14[d$rain.ave7.lag14>low_hi_cut[2]]<-"High rain"
d$raincat.21[d$rain.ave7.lag21>low_hi_cut[2]]<-"High rain"
d$raincat.7<-factor(d$raincat.7, levels=c("No rain", "Low rain", "Med rain", "High rain"))
d$raincat.14<-factor(d$raincat.14, levels=c("No rain", "Low rain", "Med rain", "High rain"))
d$raincat.21<-factor(d$raincat.21, levels=c("No rain", "Low rain", "Med rain", "High rain"))



#########################
#Summary tables.
#########################

d%>%
  group_by(raincat.7)%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meanH2s=mean(h2s), Nh2s=sum(h2s))
d%>%
  group_by(raincat.14)%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meanh2s=mean(h2s), Nh2s=sum(h2s))
d%>%
  group_by(raincat.21)%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meanh2s=mean(h2s), Nh2s=sum(h2s))


d%>%
  group_by(raincat.7, quartile = ntile(LT8, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meanh2s=mean(h2s), Nh2s=sum(h2s))
d%>%
  group_by(raincat.14, quartile = ntile(LT15, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meanh2s=mean(h2s), Nh2s=sum(h2s))
d%>%
  group_by(raincat.21, quartile = ntile(LT22, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meanh2s=mean(h2s), Nh2s=sum(h2s))



d%>%
  group_by(quartile = ntile(rain.ave7.lag7, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meanh2s=mean(h2s))
d%>%
  group_by(quartile = ntile(rain.ave7.lag14, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meanh2s=mean(h2s))
d%>%
  group_by(quartile = ntile(rain.ave7.lag21, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meanh2s=mean(h2s))





################
#Heavy rain events
################

#Set to date format for merge
PriorHeavyRain$date<-as.Date(paste(PriorHeavyRain$month,PriorHeavyRain$day,PriorHeavyRain$year, sep="/"),"%m/%d/%Y")
PriorHeavyRain<-subset(PriorHeavyRain, select=-c(year,month,day))

#merge in long term rain trends
LT <- averain %>% subset(., select=c(date, LT1, LT8, LT15, LT22))

#merge suvery and weather
dim(df)
dim(PriorHeavyRain)
dHR<-merge(df,PriorHeavyRain, by="date", all.x = T, all.y = F) 
dHR<-merge(dHR,LT, by="date", all.x = T, all.y = F)  
dim(dHR)


#Examine mean rain exposure in heavy rain and non heavy rain groups
rainmean.d <- merge(d, dHR,  by=c("date","id"), all.x = T, all.y = T)
HR.rainmean7<-rainmean.d %>% group_by(HeavyRain.lag7) %>% summarize(mean.rainfall=mean(rain.ave7.lag7, na.rm=T)) %>% as.data.frame()
HR.rainmean14<-rainmean.d %>% group_by(HeavyRain.lag14) %>% summarize(mean.rainfall=mean(rain.ave7.lag14, na.rm=T)) %>% as.data.frame()
HR.rainmean21<-rainmean.d %>% group_by(HeavyRain.lag21) %>% summarize(mean.rainfall=mean(rain.ave7.lag21, na.rm=T)) %>% as.data.frame()

save(HR.rainmean7, 
     HR.rainmean14,
     HR.rainmean21,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/villageh2s.HRmean.Results.Rdata")






##############
#Effect of heavy rain on 
#hydrogen sulfide in the 
#stored drinking water
##############


H2STMLE<-function(dat=d,i,Wvars=NULL,sl.library=library, seed=12345){
  W<-NULL
  res<-matrix(NA, nrow=1, ncol=8)
  set.seed(seed)
 colnames(dat)[which(colnames(dat)==paste0("HeavyRain.lag",i))]="HeavyRain"
 df<-dat

  a=round(sum(df$HeavyRain==1 & df$h2s==1, na.rm=T),0)
  b=round(sum(df$HeavyRain==1 & df$h2s==0, na.rm=T),0)
  c=round(sum(df$HeavyRain==0 & df$h2s==1, na.rm=T),0)
  D=round(sum(df$HeavyRain==0 & df$h2s==0, na.rm=T),0)
  cat(a,"\n",b,"\n",c,"\n",D,"\n")
  
if(!is.null(Wvars)){
  Wall<-subset(df, select=Wvars)
}
if(is.null(Wvars)){
    res<-as.numeric(washb_glm(Y=df$h2s, tr=df$HeavyRain, id=df$id, contrast=c(0,1), family=poisson(link='log'))$TR)
    res.vec<-c(a, b, c, D, res)
    names(res.vec) <- c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR", "Z", "p")

}else{
  Wall<-Wall[!is.na(df$h2s),]
  Wall<-droplevels(Wall)

  Wscreen <- washb_prescreen(Y=df$h2s[!is.na(df$h2s)],Ws=Wall,family="binomial", pval = 0.2 ,print=T)
  d.W<-subset(Wall, select=Wscreen)
  d.W<-design_matrix(d.W)
  #remove near zero variance columns
  preproc = caret::preProcess(d.W, method = c("zv", "nzv"))
  d.W = predict(preproc, d.W)
  d.W<-droplevels(d.W)
  fit<-tmle(Y=df$h2s[!is.na(df$h2s)],A=df$HeavyRain[!is.na(df$h2s)],W=d.W, family="binomial",Q.SL.library=sl.library,g.SL.library=sl.library,id=df$id[!is.na(df$h2s)])
  res[1,c(1:4)]<-c(a,b,c,D)
  res[1,5]<-fit$estimates$RR$psi
  res[1,c(6,7)]<-fit$estimates$RR$CI
  res[1,8]<-fit$estimates$RR$pvalue
  res.vec<-as.vector(res)
}
  return(res.vec)
}


H2S.Hrain.unadj<-matrix(NA, nrow=3, ncol=11)
H2S.Hrain.unadj[1,]<-H2STMLE(dat=dHR,i=1)
H2S.Hrain.unadj[2,]<-H2STMLE(dat=dHR,i=7)
H2S.Hrain.unadj[3,]<-H2STMLE(dat=dHR,i=14)
colnames(H2S.Hrain.unadj)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
H2S.Hrain.unadj





##############
#Stratified by long-term rain trends
##############



#heavy rain
  LT<-dHR$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
rainmean.d <- merge(d, dHR,  by=c("date","id"), all.x = T, all.y = T)
HR.strat.rainmean7<-rainmean.d %>% group_by(LTQ.y, HeavyRain.lag7) %>% summarize(mean.rainfall=mean(rain.ave7.lag7, na.rm=T), mean.h2s=mean(h2s.x, na.rm=T)) %>% filter(!is.na(mean.rainfall)) %>% as.data.frame()

  LT<-dHR$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
rainmean.d <- merge(d, dHR,  by=c("date","id"), all.x = T, all.y = T)
HR.strat.rainmean14<-rainmean.d %>% group_by(LTQ.y, HeavyRain.lag14) %>% summarize(mean.rainfall=mean(rain.ave7.lag14, na.rm=T), mean.h2s=mean(h2s.x, na.rm=T)) %>% filter(!is.na(mean.rainfall)) %>% as.data.frame()

  LT<-dHR$LT22
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
rainmean.d <- merge(d, dHR,  by=c("date","id"), all.x = T, all.y = T)
HR.strat.rainmean21<-rainmean.d %>% group_by(LTQ.y, HeavyRain.lag21) %>% summarize(mean.rainfall=mean(rain.ave7.lag21, na.rm=T), mean.h2s=mean(h2s.x, na.rm=T)) %>% filter(!is.na(mean.rainfall)) %>%  as.data.frame()


save(HR.strat.rainmean7, 
     HR.strat.rainmean14,
     HR.strat.rainmean21,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/villageh2s.HRmean-stratified.Results.Rdata")



#hydrogen sulfide
H2s.unadjLT3<-H2s.unadjLT2<-H2s.unadjLT1<-matrix(NA, nrow=3, ncol=11)
  LT<-dHR$LT1
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  H2s.unadjLT1[1,]<-H2STMLE(dat=dHR[dHR$LTQ==1,],i=1)
  H2s.unadjLT2[1,]<-H2STMLE(dat=dHR[dHR$LTQ==2,],i=1)
  H2s.unadjLT3[1,]<-H2STMLE(dat=dHR[dHR$LTQ==3,],i=1)
  
  LT<-dHR$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  H2s.unadjLT1[2,]<-H2STMLE(dat=dHR[dHR$LTQ==1,],i=7)
  H2s.unadjLT2[2,]<-H2STMLE(dat=dHR[dHR$LTQ==2,],i=7)
  H2s.unadjLT3[2,]<-H2STMLE(dat=dHR[dHR$LTQ==3,],i=7)
 
  
  LT<-dHR$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  H2s.unadjLT1[3,]<-H2STMLE(dat=dHR[dHR$LTQ==1,],i=14)
  H2s.unadjLT2[3,]<-H2STMLE(dat=dHR[dHR$LTQ==2,],i=14)
  H2s.unadjLT3[3,]<-H2STMLE(dat=dHR[dHR$LTQ==3,],i=14)

colnames(H2s.unadjLT3)<-colnames(H2s.unadjLT2)<-colnames(H2s.unadjLT1)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
H2S.Hrain.unadj
H2s.unadjLT1
H2s.unadjLT2
H2s.unadjLT3


#save unadjusted objects
save(
H2S.Hrain.unadj,
H2s.unadjLT1,
H2s.unadjLT2,
H2s.unadjLT3,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/villageH2S.Results.unadjusted.Rdata")

#Save dataset
save(d, 
     dHR,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/villageH2S.datasets.Rdata")




