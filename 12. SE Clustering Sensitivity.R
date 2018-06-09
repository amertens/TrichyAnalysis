

rm(list=ls())
if(!require("devtools")){
  install.packages("devtools", repos = "http://cran.us.r-project.org")
  devtools::install_github("hadley/devtools")
  library("devtools")}

if(!require("washb")){install_github("ben-arnold/washb")  ; library("washb")}
if(!require("dplyr")){install.packages("dplyr", repos = "http://cran.us.r-project.org"); library("dplyr")}
if(!require("SuperLearner")){install.packages("SuperLearner", repos = "http://cran.us.r-project.org"); library("SuperLearner")}
if(!require("tmle")){install.packages("tmle", repos = "http://cran.us.r-project.org"); library("tmle")}
if(!require("caret")){install.packages("caret", repos = "http://cran.us.r-project.org"); library("caret")}
if(!require("zoo")){install.packages("zoo", repos = "http://cran.us.r-project.org"); library("zoo")}
if(!require("cowplot")){install.packages("cowplot", repos = "http://cran.us.r-project.org"); library("cowplot")}


try(setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data"))
load("bl_covariates.Rdata")
load("LaggedWeather.Rdata")



###################
#Temperature data processing
###################
head(weather)


#Get quartiles cutpoints of overall weekly mean
weekly_temp<-rollmean(weather$avetemp,7,fill=NA, align="right")
tempQ <- quantile(weekly_temp, na.rm=T)[2:4]

#Calculate moving average with different lag times

tempvars<-sapply(1:100, function(x) paste0("at.",x))
temp<-subset(weather, select=c("year","month","day",tempvars))

#change temp vars to numeric:
temp<-apply(temp, 2, function(x) as.numeric(as.character(x))) 


ave.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) mean(x)))
  return(windowmean)
}

avetemp<-matrix(data=NA,nrow=nrow(temp),ncol=21)
ave.names<-rep(NA, 21)
for(i in 1:21){
  avetemp[,i]<-ave.window(i, 7, "at", temp)
  ave.names[i]<- paste0("temp.ave7.lag",i)
}
avetemp<-cbind(temp[,1:3],avetemp)
colnames(avetemp)[4:ncol(avetemp)]<-ave.names
avetemp<-as.data.frame(avetemp)


#Set to date format for merge
avetemp$intdate<-as.Date(paste(avetemp$month,avetemp$day,avetemp$year, sep="/"),"%m/%d/%Y")
avetemp<-subset(avetemp, select=-c(year,month,day))





###################
#Rain data processing
###################
head(weather)

#Get quartiles cutpoints of overall weekly mean
weekly_rain<-rollmean(weather$avetemp,7,fill=NA, align="right")
rainQ <- quantile(weekly_rain, na.rm=T)[2:4]

#Calculate moving average with different lag times

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

#Set date to date format for merge
averain$intdate<-as.Date(paste(averain$month,averain$day,averain$year, sep="/"),"%m/%d/%Y")
averain<-subset(averain, select=c(intdate, LT1, LT8, LT15, LT22))


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

Y<-df$diar7d
table(is.na(df$individ))
table(is.na(df$hhid))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
id<-df$vilid  #SE's clusted on village id
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

intdate<-df$intdate
Wfac<-subset(df,select= -c(vilid,hhid,individ,intdate,bdate,mdiar7d, diar2d, diar7d, diar14d, hcgi7d, diardays,stdywk, round, h2s ))
colnames(Wfac)

H2S<-df$h2s

#Check missingness of W
table(is.na(Wfac))

#Save variable names
Wvars<-colnames(Wfac)

survey<-cbind(intdate,Y,id,H2S,Wfac)



#merge suvery and weather
dim(survey)
dim(avetemp)
d<-merge(survey,avetemp, by="intdate", all.x = F, all.y = F) 
dim(d)
d<-merge(d,averain, by="intdate", all.x = F, all.y = F) 
dim(d)




####################
#Quartile analysis function
####################


quart_tmle<-function(d,Wvars=NULL,weather="temp.ave7.lag7",q=c(1,4),cutoff=4, library=NULL){
colnames(d)[which( colnames(d)==paste0(weather))]="weather"

d$quartile<-cut(d$weather, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))

df <- d[d$quartile==levels(d$quartile)[q[1]]|
           d$quartile==levels(d$quartile)[q[2]],]
df$A=ifelse(df$quartile==levels(df$quartile)[q[2]],1,0)

print(table(df$Y, df$quartile))

if(!is.null(Wvars)){
  Wall<-subset(df, select=Wvars)
}


if(is.null(Wvars)){
  res<-as.numeric(washb_glm(Y=df$Y, tr=df$A, id=df$id, contrast=c(0,1), family=poisson(link='log'))$TR)
        
}else{
 
  Wscreen <- washb_prescreen(Y=df$Y,Ws=Wall,family="binomial", pval = 0.2 ,print=T)
  W<-subset(Wall, select=Wscreen)
  W<-design_matrix(W)
  #remove near zero variance columns
  preproc = caret::preProcess(W, method = c("zv", "nzv"))
  W = predict(preproc, W)
  W<-droplevels(W)
  d<-cbind(df$Y,df$A,W)
  d<-na.omit(d)


      res<-tmle(Y=d[,1],A=d[,2],W=d[,-c(1:2)], family="binomial",Q.SL.library=library,g.SL.library=library,id=df$id)
          }  

  return(res)
}



####################
#Quartile analysis function
####################


quart_tmle<-function(d,Wvars=NULL,weather="temp.ave7.lag7",q=c(1,4),cutoff=4, library=NULL){
colnames(d)[which( colnames(d)==paste0(weather))]="weather"

# df<-d%>%
#   mutate(quartile = ntile(weather, cutoff))%>%
#   subset(quartile==q[1]|quartile==q[2])%>%
#   mutate(A=ifelse(quartile==q[2],1,0))

d$quartile<-cut(d$weather, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))

df <- d[d$quartile==levels(d$quartile)[q[1]]|
           d$quartile==levels(d$quartile)[q[2]],]
df$A=ifelse(df$quartile==levels(df$quartile)[q[2]],1,0)

print(table(df$Y, df$quartile))

if(!is.null(Wvars)){
  Wall<-subset(df, select=Wvars)
}


if(is.null(Wvars)){
  res<-as.numeric(washb_glm(Y=df$Y, tr=df$A, id=df$id, contrast=c(0,1), family=poisson(link='log'))$TR)
        
}else{
 
  Wscreen <- washb_prescreen(Y=df$Y,Ws=Wall,family="binomial", pval = 0.2 ,print=T)
  W<-subset(Wall, select=Wscreen)
  W<-design_matrix(W)
  #remove near zero variance columns
  preproc = caret::preProcess(W, method = c("zv", "nzv"))
  W = predict(preproc, W)
  W<-droplevels(W)
  d<-cbind(df$Y,df$A,W)
  d<-na.omit(d)


      res<-tmle(Y=d[,1],A=d[,2],W=d[,-c(1:2)], family="binomial",Q.SL.library=library,g.SL.library=library,id=df$id)
          }  

  return(res)
}



####################
#Unadjusted MH tests
####################

temp.unadj1v4<-temp.unadj2v4<-temp.unadj1v2<-temp.unadj1v3<-temp.unadj3v4<-matrix(NA, nrow=3, ncol=7)
for(j in 1:3){
  i<-j*7
  set.seed(12345)
  temp.unadj1v4[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,4),cutoff=4)
  temp.unadj2v4[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(2,4),cutoff=4)
  temp.unadj1v2[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,2),cutoff=4)
  temp.unadj1v3[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,3),cutoff=4)
  temp.unadj3v4[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(3,4),cutoff=4)
}
temp.unadj<-NULL
for(i in 1:3){
  temp.unadj<-rbind(temp.unadj,temp.unadj1v2[i,],temp.unadj1v3[i,],temp.unadj1v4[i,])
}
rownames(temp.unadj)<-c("lag7 1v2","lag7 1v3","lag7 1v4",
                        "lag14 1v2","lag14 1v3","lag14 1v4",
                        "lag21 1v2","lag21 1v3","lag21 1v4")
colnames(temp.unadj)<- c("PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
temp.unadj

  
  
  
##########################
#Unadjusted analysis -Heavy rain
##########################
  


#merge survey and weather
dim(survey)
dim(averain)
d<-merge(survey,averain, by="intdate", all.x = F, all.y = F) 
dim(d)

colnames(d)


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
  


Hrain<-matrix(NA, nrow=3, ncol=11)
Hrain[1,]<-HeavyRainTMLE(dat=dHR,i=7)
Hrain[2,]<-HeavyRainTMLE(dat=dHR,i=14)
Hrain[3,]<-HeavyRainTMLE(dat=dHR,i=21)

colnames(Hrain) <- c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")



#heavy rain - stratified
Hrain.unadjLT3<-Hrain.unadjLT2<-Hrain.unadjLT1<-matrix(NA, nrow=3, ncol=11)
  LT<-dHR$LT8
   ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
   dHR$LTQ<-2
   dHR$LTQ[LT < ntiles[1]]<-1
   dHR$LTQ[LT >= ntiles[2]]<-3
  Hrain.unadjLT1[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==1,],i=7)
  Hrain.unadjLT2[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==2,],i=7)
  Hrain.unadjLT3[1,]<-HeavyRainTMLE(dat=dHR[dHR$LTQ==3,],i=7)
  
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





##############
#Effect of heavy rain on 
#hydrogen sulfide in the 
#stored drinking water
##############

dHR$h2s <- dHR$H2S
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











#Plot parameters
theme_set(theme_bw())

#Temp pallete
tempcol <- c("#000000","#ED8021", "#E15324", "#D62728")
  
#Rain pallete
raincol <- c("#000000",  "#56B4E9", "#4C8FE5", "#426AE3")



scaleFUN <- function(x) sprintf("%.2f", x)




raindf <- data.frame(exposure=rep("Heavy rain",12),
             strata=c(rep("Unstratified",3),
             rep("Low",3),
             rep("Medium",3),
             rep("High",3)),
             x=as.character(c(1,1,1,2,2,2,3,3,3,4,4,4)),
             lag=rep(c("One week lag", "Two week lag", "Three week lag"),4),
              rbind(
              Hrain,
              Hrain.unadjLT1,
              Hrain.unadjLT2,
              Hrain.unadjLT3))

raindf$strata <- factor(raindf$strata, levels=unique(raindf$strata))
raindf$lag <- factor(raindf$lag, levels=unique(raindf$lag))

raindf<- subset(raindf, select=-c(a,b,c,d, logPR,se.logPR,Z,p))

#Add blank rows for reference levels
temp.unadj <- rbind(NA,temp.unadj[1:3,], NA,temp.unadj[4:6,], NA,temp.unadj[7:9,])
rownames(temp.unadj) <- NULL
tempdf <- data.frame(exposure=rep("Temperature",12),
             strata=rep(c("(ref.)","1v2","1v3","1v4"),3),
                          x=as.character(c(1,2,3,4,1,2,3,4,1,2,3,4)),
             lag=c(rep("One week lag",4), rep("Two week lag",4), rep("Three week lag",4)),
                         temp.unadj) %>% 
  subset(., select=-c(logPR,se.logPR,Z,p))
tempdf$strata <- factor(tempdf$strata, levels=unique(tempdf$strata))
tempdf$lag <- factor(tempdf$lag, levels=unique(tempdf$lag))


H2Sdf <- data.frame(exposure=rep("Heavy rain",12),
             strata=c(rep("Unstratified",3),
             rep("Low",3),
             rep("Medium",3),
             rep("High",3)),
             x=as.character(c(1,1,1,2,2,2,3,3,3,4,4,4)),
             lag=rep(c("One week lag", "Two week lag", "Three week lag"),4),
              rbind(
              H2S.Hrain.unadj,
              H2s.unadjLT1,
              H2s.unadjLT2,
              H2s.unadjLT3))
H2Sdf<- H2Sdf %>% rename(PR=PR.)
H2Sdf$strata <- factor(H2Sdf$strata, levels=unique(H2Sdf$strata))
H2Sdf$lag <- factor(H2Sdf$lag, levels=unique(H2Sdf$lag))

H2Sdf<- subset(H2Sdf, select=-c(a,b,c,d, logPR,se.logPR,Z,p))



#--------------------------------------
#plot rain
#--------------------------------------


  
    yticks=c(0.125,0.25,0.5,1,2,4,8)
  
  
prain<-ggplot(raindf, aes(x=strata)) + 
  geom_point(aes(y=PR, fill=strata, color=strata), size = 4) +
  geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strata),
                 alpha=0.5, size = 3) +
  labs(x = "Long-term rainfall strata", y = "Prevalence Ratio") +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim=range(yticks)) +
  scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
  scale_fill_manual(values=rep(raincol,4)) +
  scale_colour_manual(values=rep(raincol,4)) +
  scale_size_continuous(range = c(0.5, 1))+
  theme(strip.background = element_blank(),
    legend.position="none",
    strip.text.x = element_text(size=12),
    axis.text.x = element_text(size=12)) +
  facet_wrap(~lag,  scales = "fixed") +
  ggtitle("Rain")
prain



#--------------------------------------
#plot H2S-rain
#--------------------------------------
  
    yticks=c(0.125,0.25,0.5,1,2,4,8)
  
  
pH2S<-ggplot(H2Sdf, aes(x=strata)) + 
  geom_point(aes(y=PR, fill=strata, color=strata), size = 4) +
  geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strata),
                 alpha=0.5, size = 3) +
  labs(x = "Long-term rainfall strata", y = "Prevalence Ratio") +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim=range(yticks)) +
  scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
  scale_fill_manual(values=rep(raincol,4)) +
  scale_colour_manual(values=rep(raincol,4)) +
  scale_size_continuous(range = c(0.5, 1))+
  theme(strip.background = element_blank(),
    legend.position="none",
    strip.text.x = element_text(size=12),
    axis.text.x = element_text(size=12)) +
  facet_wrap(~lag,  scales = "fixed") +
  ggtitle("H2S")
pH2S

#--------------------------------------
#plot temperature
#--------------------------------------


    yticks=c(0.125,0.25,0.5,1,2,4,8)
  
  
ptemp<-ggplot(tempdf, aes(x=strata)) + 
  geom_point(aes(y=PR, fill=strata, color=strata), size = 4) +
  geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strata),
                 alpha=0.5, size = 3) +
  labs(x = "Quartile contrast", y = "Prevalence Ratio") +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim=range(yticks)) +
  scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
  scale_fill_manual(values=rep(tempcol,4)) +
  scale_colour_manual(values=rep(tempcol,4)) +
  scale_size_continuous(range = c(0.5, 1))+
  theme(strip.background = element_blank(),
    legend.position="none",
    strip.text.x = element_text(size=12),
    axis.text.x = element_text(size=12)) +
  facet_wrap(~lag,  scales = "fixed") +
  ggtitle("Temperature")
ptemp




  
#------------------------------
# Combine plots into one figure
#------------------------------
  
  
  p <- plot_grid(ptemp, prain, pH2S, labels = c("A", "B","C"), ncol = 1)
  p
  
  ggsave(p, file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Clusterid_sensitivity_figure.pdf")
  
  
      