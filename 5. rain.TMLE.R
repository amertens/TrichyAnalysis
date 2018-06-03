
rm(list=ls())

library(dplyr)
library(SuperLearner)
library(tmle)
library(washb)
library(caret)


setwd("C:/Users/andre/Documents/Trichy analysis/Cleaned Data")
load("bl_covariates.Rdata")
load("LaggedWeather.Rdata")



###################
#Weather data processing
###################
head(weather)

rainvars<-sapply(1:100, function(x) paste0("rain.",x))
rain<-subset(weather, select=c("year","month","day",rainvars))

#Change "trace" to 0.1 in rainfall variables
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

#Save heavy rain threshold
save(HeavyRainThres, 
     file="C:/Users/andre/Documents/Trichy analysis/Results/HRthres.Rdata")



HeavyRain<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres,1,0))
table(HeavyRain)

 

HeavyRain.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) ifelse(sum(x)>0,1,0)))
  return(windowmean)
}


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
#Survey data processing
###################

# Expand out factors into indicators
dim(df)
df<-df[,-5] #remove duplicate "round" variable
head(df)

length(unique(df$vilid)) #25 villages
length(unique(df$hhid)) #900 households
length(unique(df$individ)) # 1284 children
num_obs<-df %>% group_by(individ) %>% summarize(n=n())
mean(num_obs$n, na.rm=T)

h<-df %>% filter(!is.na(h2s))
length(unique(h$vilid)) #25 villages
length(unique(h$hhid)) #900 households
num_obs<-h %>% group_by(hhid) %>% summarize(n=n())
mean(num_obs$n, na.rm=T)

Y<-df$diar7d
h2s<-df$h2s
id<-df$hhid #SE clustering on hhid 
intdate<-df$intdate
Wfac<-subset(df,select= -c(vilid,hhid,individ,intdate,bdate,mdiar7d, diar2d, diar7d, diar14d, hcgi7d, diardays,stdywk, round, h2s ))
colnames(Wfac)

#Check missingness of W
table(is.na(Wfac))

#Save variable names
Wvars<-colnames(Wfac)



survey<-cbind(intdate,Y,id,h2s,Wfac)


##########################
#Unadjusted analysis
##########################

#merge suvery and weather
dim(survey)
dim(averain)
d<-merge(survey,averain, by="intdate", all.x = F, all.y = F) 
dim(d)

colnames(d)



#Mean diarrhea by quartile:
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
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meandiar=mean(Y), Ndiar=sum(Y))
d%>%
  group_by(raincat.14)%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meandiar=mean(Y), Ndiar=sum(Y))
d%>%
  group_by(raincat.21)%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meandiar=mean(Y), Ndiar=sum(Y))


d%>%
  group_by(raincat.7, quartile = ntile(LT8, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meandiar=mean(Y), Ndiar=sum(Y))
d%>%
  group_by(raincat.14, quartile = ntile(LT15, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meandiar=mean(Y), Ndiar=sum(Y))
d%>%
  group_by(raincat.21, quartile = ntile(LT22, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meandiar=mean(Y), Ndiar=sum(Y))



d%>%
  group_by(quartile = ntile(rain.ave7.lag7, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meandiar=mean(Y))
d%>%
  group_by(quartile = ntile(rain.ave7.lag14, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meandiar=mean(Y))
d%>%
  group_by(quartile = ntile(rain.ave7.lag21, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meandiar=mean(Y))





####################
#Unadjusted MH tests
####################

rain_tmle<-function(d,Wvars=NULL,weather="rain.ave7.lag7",q=c(1,4),cutoff=4, library=NULL){
colnames(d)[which( colnames(d)==paste0(weather))]="weather"

#Quantile non-zero rain days
ntiles<-as.vector(quantile(d$weather[d$weather!=0], probs= seq(0, 1, 0.3333))[c(2,3)])

d$Q<-1
d$Q[d$weather>0]<-2
d$Q[d$weather>ntiles[1]]<-3
d$Q[d$weather>ntiles[2]]<-4
table(d$Q)

df<-d %>%
  subset(Q==q[1]|Q==q[2]) %>%
  mutate(A=ifelse(Q==q[2],1,0))

  a=round(sum(df$A==1 & df$Y==1, na.rm=T),0)
  b=round(sum(df$A==1 & df$Y==0, na.rm=T),0)
  c=round(sum(df$A==0 & df$Y==1, na.rm=T),0)
  d=round(sum(df$A==0 & df$Y==0, na.rm=T),0)
  
if(is.null(Wvars)){
    res<-as.numeric(washb_glm(Y=df$Y, tr=df$A, id=df$id, contrast=c(0,1), family=binomial(link='log'))$TR)
    res<-c(a, b, c, d, res)
    names(res) <- c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR", "Z", "p")
}else{
      res<-tmle(Y=df$Y,A=df$A,W=W, family="binomial",Q.SL.library=library,g.SL.library=library,id=df$id)
  }  
  return(res)
}

rain.unadj1v4<-rain.unadj2v4<-rain.unadj1v2<-rain.unadj1v3<-rain.unadj3v4<-matrix(NA, nrow=3, ncol=11)
for(j in 1:3){
  i<-j*7
  set.seed(12345)
  rain.unadj1v2[j,]<-rain_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,2),cutoff=4)
  rain.unadj1v3[j,]<-rain_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,3),cutoff=4)
  rain.unadj1v4[j,]<-rain_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,4),cutoff=4)
  rain.unadj2v4[j,]<-rain_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(2,4),cutoff=4)
  rain.unadj3v4[j,]<-rain_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(3,4),cutoff=4)
}

rain.unadj1v2
rain.unadj1v3
rain.unadj1v4
rain.unadj2v4
rain.unadj3v4

rain.unadj<-NULL
for(i in 1:3){
  rain.unadj<-rbind(rain.unadj,rain.unadj1v2[i,],rain.unadj1v3[i,],rain.unadj1v4[i,])
}
rownames(rain.unadj)<-c("lag7 1v2","lag7 1v3","lag7 1v4",
                        "lag14 1v2","lag14 1v3","lag14 1v4",
                        "lag21 1v2","lag21 1v3","lag21 1v4")
colnames(rain.unadj)<- c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
rain.unadj




################
#Heavy rain events
################

#Set to date format for merge
PriorHeavyRain$intdate<-as.Date(paste(PriorHeavyRain$month,PriorHeavyRain$day,PriorHeavyRain$year, sep="/"),"%m/%d/%Y")
PriorHeavyRain<-subset(PriorHeavyRain, select=-c(year,month,day))

#merge in long term rain trends
LT <- averain %>% subset(., select=c(intdate, LT1, LT8, LT15, LT22))

#merge suvery and weather
dim(survey)
dim(PriorHeavyRain)
dHR<-merge(survey,PriorHeavyRain, by="intdate", all.x = F, all.y = F)  
dHR<-merge(dHR,LT, by="intdate", all.x = T, all.y = F) 
dim(dHR)


#Examine mean rain exposure in heavy rain and non heavy rain groups

rainmean.d <- merge(d, dHR,  by=c("intdate","id"), all.x = T, all.y = T)
HR.rainmean7<-rainmean.d %>% group_by(HeavyRain.lag7) %>% summarize(mean.rainfall=mean(rain.ave7.lag7, na.rm=T)) %>% as.data.frame()
HR.rainmean14<-rainmean.d %>% group_by(HeavyRain.lag14) %>% summarize(mean.rainfall=mean(rain.ave7.lag14, na.rm=T)) %>% as.data.frame()
HR.rainmean21<-rainmean.d %>% group_by(HeavyRain.lag21) %>% summarize(mean.rainfall=mean(rain.ave7.lag21, na.rm=T)) %>% as.data.frame()

save(HR.rainmean7, 
     HR.rainmean14,
     HR.rainmean21,
     file="C:/Users/andre/Documents/Trichy analysis/Results/HRmean.Results.Rdata")




#Heavy rain TMLE function
HeavyRainTMLE<-function(dat=d,i,Wvars=NULL,sl.library=library){
  W<-NULL
  res<-matrix(NA, nrow=1, ncol=8)
  set.seed(12345)
 colnames(dat)[which(colnames(dat)==paste0("HeavyRain.lag",i))]="HeavyRain"

 df<-dat

  a=round(sum(df$HeavyRain==1 & df$Y==1, na.rm=T),0)
  b=round(sum(df$HeavyRain==1 & df$Y==0, na.rm=T),0)
  c=round(sum(df$HeavyRain==0 & df$Y==1, na.rm=T),0)
  D=round(sum(df$HeavyRain==0 & df$Y==0, na.rm=T),0)
  
if(!is.null(Wvars)){
  Wall<-subset(df, select=Wvars)
}
if(is.null(Wvars)){
  
      res<-as.numeric(washb_glm(Y=df$Y, tr=df$HeavyRain, id=df$id, contrast=c(0,1), family=poisson(link='log'))$TR)
    res.vec<-c(a, b, c, D, res)
    names(res.vec) <- c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR", "Z", "p")

}else{
  Keep<-rep(T, ncol(Wall))
  for(i in 1:ncol(Wall)){
    if(length(levels(Wall[,i]))==1){
      print(table(Wall$wqsource))
      Keep[i]<-F
    }
  }
  
  Wall<-droplevels(Wall)
  Wscreen <- washb_prescreen(Y=df$Y,Ws=Wall,family="binomial", pval = 0.2 ,print=T)
  d.W<-subset(Wall, select=Wscreen)
  d.W<-design_matrix(d.W)
  #remove near zero variance columns
  preproc = caret::preProcess(d.W, method = c("zv", "nzv"))
  d.W = predict(preproc, d.W)
  fit<-tmle(Y=df$Y,A=df$HeavyRain,W=d.W, family="binomial",Q.SL.library=sl.library,g.SL.library=sl.library,id=df$id)
  res[1,c(1:4)]<-c(a,b,c,D)
  res[1,5]<-fit$estimates$RR$psi
  res[1,c(6,7)]<-fit$estimates$RR$CI
  res[1,8]<-fit$estimates$RR$pvalue
  res.vec<-as.vector(res)
}
  return(res.vec)
}
  

Hrain.unadj<-matrix(NA, nrow=3, ncol=11)
Hrain.unadj[1,]<-HeavyRainTMLE(dat=dHR,i=7)
Hrain.unadj[2,]<-HeavyRainTMLE(dat=dHR,i=14)
Hrain.unadj[3,]<-HeavyRainTMLE(dat=dHR,i=21)
colnames(Hrain.unadj)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.unadj





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

rain.unadj1v4LT1<-rain.unadj2v4LT1<-rain.unadj1v2LT1<-rain.unadj1v3LT1<-rain.unadj3v4LT1<-matrix(NA, nrow=3, ncol=11)
rain.unadj1v4LT2<-rain.unadj2v4LT2<-rain.unadj1v2LT2<-rain.unadj1v3LT2<-rain.unadj3v4LT2<-matrix(NA, nrow=3, ncol=11)
rain.unadj1v4LT3<-rain.unadj2v4LT3<-rain.unadj1v2LT3<-rain.unadj1v3LT3<-rain.unadj3v4LT3<-matrix(NA, nrow=3, ncol=11)

for(j in 1:3){
  LT<-d$LT8
  if(j==2){LT<-d$LT15}
  if(j==3){LT<-d$LT22}

  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  d$LTQ<-2
  d$LTQ[LT < ntiles[1]]<-1
  d$LTQ[LT >= ntiles[2]]<-3
  
  
  i<-j*7
  set.seed(12345)
  rain.unadj1v2LT1[j,]<-rain_tmle(d=d[d$LTQ==1,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,2),cutoff=4)
  rain.unadj1v3LT1[j,]<-rain_tmle(d=d[d$LTQ==1,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,3),cutoff=4)
  rain.unadj1v4LT1[j,]<-rain_tmle(d=d[d$LTQ==1,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,4),cutoff=4)
  rain.unadj2v4LT1[j,]<-rain_tmle(d=d[d$LTQ==1,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(2,4),cutoff=4)
  rain.unadj3v4LT1[j,]<-rain_tmle(d=d[d$LTQ==1,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(3,4),cutoff=4)
  
  rain.unadj1v2LT2[j,]<-rain_tmle(d=d[d$LTQ==2,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,2),cutoff=4)
  rain.unadj1v3LT2[j,]<-rain_tmle(d=d[d$LTQ==2,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,3),cutoff=4)
  rain.unadj1v4LT2[j,]<-rain_tmle(d=d[d$LTQ==2,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,4),cutoff=4)
  rain.unadj2v4LT2[j,]<-rain_tmle(d=d[d$LTQ==2,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(2,4),cutoff=4)
  rain.unadj3v4LT2[j,]<-rain_tmle(d=d[d$LTQ==2,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(3,4),cutoff=4)
  
  rain.unadj1v2LT3[j,]<-rain_tmle(d=d[d$LTQ==3,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,2),cutoff=4)
  rain.unadj1v3LT3[j,]<-rain_tmle(d=d[d$LTQ==3,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,3),cutoff=4)
  rain.unadj1v4LT3[j,]<-rain_tmle(d=d[d$LTQ==3,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,4),cutoff=4)
  rain.unadj2v4LT3[j,]<-rain_tmle(d=d[d$LTQ==3,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(2,4),cutoff=4)
  rain.unadj3v4LT3[j,]<-rain_tmle(d=d[d$LTQ==3,],W=NULL,weather=paste0("rain.ave7.lag",i),q=c(3,4),cutoff=4)
}


rain.unadj1v4LT1
rain.unadj1v4LT2
rain.unadj1v4LT3


rain.unadjLT3<-rain.unadjLT2<-rain.unadjLT1<-NULL
for(i in 1:3){
  rain.unadjLT1<-rbind(rain.unadjLT1,rain.unadj1v2LT1[i,],rain.unadj1v3LT1[i,],rain.unadj1v4LT1[i,])
  rain.unadjLT2<-rbind(rain.unadjLT2,rain.unadj1v2LT2[i,],rain.unadj1v3LT2[i,],rain.unadj1v4LT2[i,])
  rain.unadjLT3<-rbind(rain.unadjLT3,rain.unadj1v2LT3[i,],rain.unadj1v3LT3[i,],rain.unadj1v4LT3[i,])

}
rownames(rain.unadjLT3)<-rownames(rain.unadjLT2)<-rownames(rain.unadjLT1)<-c(
                        "lag7 1v2","lag7 1v3","lag7 1v4",
                        "lag14 1v2","lag14 1v3","lag14 1v4",
                        "lag21 1v2","lag21 1v3","lag21 1v4")
colnames(rain.unadjLT3)<-colnames(rain.unadjLT2)<-colnames(rain.unadjLT1)<- c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")

rain.unadjLT1
rain.unadjLT2
rain.unadjLT3




#heavy rain
Hrain.unadjLT3<-Hrain.unadjLT2<-Hrain.unadjLT1<-matrix(NA, nrow=3, ncol=11)
  LT<-dHR$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
rainmean.d <- merge(d, dHR,  by=c("intdate","id"), all.x = T, all.y = T)
HR.strat.rainmean7<-rainmean.d %>% group_by(LTQ.y, HeavyRain.lag7) %>% summarize(mean.rainfall=mean(rain.ave7.lag7, na.rm=T), mean.h2s=mean(h2s.x, na.rm=T)) %>% as.data.frame()
  Hrain.unadjLT1[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==1 & (d$rain.ave7.lag7==0 | dHR$HeavyRain.lag7==1),],i=7)
  Hrain.unadjLT2[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==2 & (d$rain.ave7.lag7==0 | dHR$HeavyRain.lag7==1),],i=7)
  Hrain.unadjLT3[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==3 & (d$rain.ave7.lag7==0 | dHR$HeavyRain.lag7==1),],i=7)
  
  LT<-dHR$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
rainmean.d <- merge(d, dHR,  by=c("intdate","id"), all.x = T, all.y = T)
HR.strat.rainmean14<-rainmean.d %>% group_by(LTQ.y, HeavyRain.lag14) %>% summarize(mean.rainfall=mean(rain.ave7.lag14, na.rm=T), mean.h2s=mean(h2s.x, na.rm=T)) %>% as.data.frame()
  Hrain.unadjLT1[2,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==1,],i=14)
  Hrain.unadjLT2[2,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==2,],i=14)
  Hrain.unadjLT3[2,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==3,],i=14)
  
  LT<-dHR$LT22
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
rainmean.d <- merge(d, dHR,  by=c("intdate","id"), all.x = T, all.y = T)
HR.strat.rainmean21<-rainmean.d %>% group_by(LTQ.y, HeavyRain.lag21) %>% summarize(mean.rainfall=mean(rain.ave7.lag21, na.rm=T), mean.h2s=mean(h2s.x, na.rm=T)) %>% as.data.frame()
  Hrain.unadjLT1[3,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==1 ,],i=21)
  Hrain.unadjLT2[3,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==2 ,],i=21)
  Hrain.unadjLT3[3,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==3 ,],i=21)

colnames(Hrain.unadjLT3)<-colnames(Hrain.unadjLT2)<-colnames(Hrain.unadjLT1)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.unadjLT1
Hrain.unadjLT2
Hrain.unadjLT3

save(HR.strat.rainmean7, 
     HR.strat.rainmean14,
     HR.strat.rainmean21,
     file="C:/Users/andre/Documents/Trichy analysis/Results/HRmean-stratified.Results.Rdata")





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
H2s.unadjLT1
H2s.unadjLT2
H2s.unadjLT3


#save unadjusted objects
save(rain.unadj, 
     rain.unadjLT1,
rain.unadjLT2,
rain.unadjLT3,
Hrain.unadj,
Hrain.unadjLT1,
Hrain.unadjLT2,
Hrain.unadjLT3,
H2S.Hrain.unadj,
H2s.unadjLT1,
H2s.unadjLT2,
H2s.unadjLT3,
     file="C:/Users/andre/Documents/Trichy analysis/Results/rain.Results.unadjusted.Rdata")

#Save dataset
save(d, 
     dHR,
     file="C:/Users/andre/Documents/Trichy analysis/Results/rain.datasets.Rdata")





##############
#Adjusted analysis
##############



####################
#Adjusted TMLE
####################


library=c("SL.mean","SL.glm","SL.glmnet", "SL.gam", "SL.bayesglm")


dHR$intdate<-as.numeric(dHR$intdate)


HR.adj<-matrix(NA, nrow<-3, ncol=8)
  set.seed(12345)
  HR.adj[1,]<-HeavyRainTMLE(d=dHR,Wvars=Wvars,i=7,sl.library=library)
  HR.adj[2,]<-HeavyRainTMLE(d=dHR,Wvars=Wvars,i=14,sl.library=library)
  HR.adj[3,]<-HeavyRainTMLE(d=dHR,Wvars=Wvars,i=21,sl.library=library)
colnames(HR.adj)<-c("a","b","c","d","PR","ci.lb","ci.ub","pvalue")
  rownames(HR.adj)<-c("1week.lag", "2week.lag", "3week.lag")
Hrain.adj<-HR.adj
Hrain.adj



H2S.Hrain.adj<-matrix(NA, nrow=3, ncol=8)
H2S.Hrain.adj[1,]<-H2STMLE(dat=dHR[!is.na(dHR$h2s),],Wvars=Wvars,i=1, sl.library=library)
H2S.Hrain.adj[2,]<-H2STMLE(dat=dHR[!is.na(dHR$h2s),],Wvars=Wvars,i=7, sl.library=library)
H2S.Hrain.adj[3,]<-H2STMLE(dat=dHR[!is.na(dHR$h2s),],Wvars=Wvars,i=14, sl.library=library)
colnames(H2S.Hrain.adj)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "pval")
H2S.Hrain.adj




####################
#Adjusted subgroup analysis
####################
#heavy rain
Hrain.adjLT3<-Hrain.adjLT2<-Hrain.adjLT1<-matrix(NA, nrow=3, ncol=8)
  LT<-dHR$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  Hrain.adjLT1[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==1,],i=7,Wvars=Wvars,sl.library=library)
  Hrain.adjLT2[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==2,],i=7,Wvars=Wvars,sl.library=library)
  Hrain.adjLT3[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==3,],i=7,Wvars=Wvars,sl.library=library)
  
  LT<-dHR$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  Hrain.adjLT1[2,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==1,],i=14,Wvars=Wvars,sl.library=library)
  Hrain.adjLT2[2,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==2,],i=14,Wvars=Wvars,sl.library=library)
  Hrain.adjLT3[2,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==3,],i=14,Wvars=Wvars,sl.library=library)
  
  LT<-dHR$LT22
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  Hrain.adjLT1[3,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==1,],i=21,Wvars=Wvars,sl.library=library)
  Hrain.adjLT2[3,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==2,],i=21,Wvars=Wvars,sl.library=library)
  Hrain.adjLT3[3,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==3,],i=21,Wvars=Wvars,sl.library=library)

colnames(Hrain.adjLT3)<-colnames(Hrain.adjLT2)<-colnames(Hrain.adjLT1)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub",   "p")
Hrain.adjLT1
Hrain.adjLT2
Hrain.adjLT3


#h2s
H2S.adjLT3<-H2S.adjLT2<-H2S.adjLT1<-matrix(NA, nrow=3, ncol=8)
  LT<-dHR$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  H2S.adjLT1[1,]<-H2STMLE(dat=dHR[dHR$LTQ==1,],i=1,Wvars=Wvars,sl.library=library)
  H2S.adjLT2[1,]<-H2STMLE(dat=dHR[dHR$LTQ==2,],i=1,Wvars=Wvars,sl.library=library)
  H2S.adjLT3[1,]<-H2STMLE(dat=dHR[dHR$LTQ==3,],i=1,Wvars=Wvars,sl.library=library)
  
  LT<-dHR$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  H2S.adjLT1[2,]<-H2STMLE(dat=dHR[dHR$LTQ==1,],i=7,Wvars=Wvars,sl.library=library)
  H2S.adjLT2[2,]<-H2STMLE(dat=dHR[dHR$LTQ==2,],i=7,Wvars=Wvars,sl.library=library)
  H2S.adjLT3[2,]<-H2STMLE(dat=dHR[dHR$LTQ==3,],i=7,Wvars=Wvars,sl.library=library)
  
  LT<-dHR$LT22
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  dHR$LTQ<-2
  dHR$LTQ[LT < ntiles[1]]<-1
  dHR$LTQ[LT >= ntiles[2]]<-3
  H2S.adjLT1[3,]<-H2STMLE(dat=dHR[dHR$LTQ==1,],i=14,Wvars=Wvars,sl.library=library)
  H2S.adjLT2[3,]<-H2STMLE(dat=dHR[dHR$LTQ==2,],i=14,Wvars=Wvars,sl.library=library)
  H2S.adjLT3[3,]<-H2STMLE(dat=dHR[dHR$LTQ==3,],i=14,Wvars=Wvars,sl.library=library)

colnames(H2S.adjLT3)<-colnames(H2S.adjLT2)<-colnames(H2S.adjLT1)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub",   "p")
H2S.adjLT1
H2S.adjLT2
H2S.adjLT3


save(
Hrain.adj,
Hrain.adjLT1,
Hrain.adjLT2,
Hrain.adjLT3,
H2S.Hrain.adj,
     file="C:/Users/andre/Documents/Trichy analysis/Results/rain.Results.adjusted.Rdata")

save(
  H2S.adjLT1,
H2S.adjLT2,
H2S.adjLT3,
     file="C:/Users/andre/Documents/Trichy analysis/Results/H2S.Results.adjusted.Rdata")








