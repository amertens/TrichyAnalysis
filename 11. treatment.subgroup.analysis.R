


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
id<-df$hhid  #SE's clusted on household
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



#Subset dataset into control and intervention arms
dc <- d %>% filter(wpi==0)
di <- d %>% filter(wpi==1)
head(d)

summary(glm(d$Y~ d$wpi))

#Effect modification?
summary(glm(d$Y~ d$wpi * d$temp.ave7.lag7, family="binomial"))
summary(glm(d$Y~ d$wpi * d$temp.ave7.lag14, family="binomial"))
summary(glm(d$Y~ d$wpi * d$temp.ave7.lag21, family="binomial"))
summary(glm(d$Y~ d$wpi * d$rain.ave7.lag7, family="binomial"))
summary(glm(d$Y~ d$wpi * d$rain.ave7.lag14, family="binomial"))
summary(glm(d$Y~ d$wpi * d$rain.ave7.lag21, family="binomial"))


 wpi_lrtest<- function (dat, A, Y="Y", family = "binomial", print=F){
  require(lmtest)
   dat$A <- dat[,A]

    fit1 <- glm(Y ~ A*wpi, data = dat, family = family)
    
    if(print==T){print(summary(fit1))}
    
    fit0 <- glm(Y ~ A, data = dat, family = family)
    LRp <- lrtest(fit1, fit0)[2, 5]

    return(LRp)
}

#Continious exposures
# wpi_lrtest(d, A="temp.ave7.lag7")
# wpi_lrtest(d, A="temp.ave7.lag14")
# wpi_lrtest(d, A="temp.ave7.lag21")
# 
# wpi_lrtest(d, A="rain.ave7.lag7")
# wpi_lrtest(d, A="rain.ave7.lag14")
# wpi_lrtest(d, A="rain.ave7.lag21")


#Categorized temperature exposures
d$temp.ave7.lag7Q<-cut(d$temp.ave7.lag7, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
d$temp.ave7.lag14Q<-cut(d$temp.ave7.lag14, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
d$temp.ave7.lag21Q<-cut(d$temp.ave7.lag21, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))

wpi_lrtest(d, A="temp.ave7.lag7Q")
wpi_lrtest(d, A="temp.ave7.lag14Q")
wpi_lrtest(d, A="temp.ave7.lag21Q")

#Heavy rain
wpi_lrtest(dHR, A="HeavyRain.lag7")
wpi_lrtest(dHR, A="HeavyRain.lag14")
wpi_lrtest(dHR, A="HeavyRain.lag21")


    fit1 <- glm(Y ~ temp.ave7.lag7Q*wpi + temp.ave7.lag14Q*wpi + temp.ave7.lag21Q*wpi, data = d, family = "binomial")
    fit0 <- glm(Y ~ temp.ave7.lag7Q + temp.ave7.lag14Q + temp.ave7.lag21Q, data = d, family = "binomial")
    LRp <- lrtest(fit1, fit0)[2, 5]
    LRp
    
    fit1 <- glm(Y ~ HeavyRain.lag7*wpi + HeavyRain.lag14*wpi + HeavyRain.lag21*wpi, data = dHR, family = "binomial")
    fit0 <- glm(Y ~ HeavyRain.lag7 + HeavyRain.lag14 + HeavyRain.lag21, data = dHR, family = "binomial")
    LRp <- lrtest(fit1, fit0)[2, 5]
   LRp

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




##########################
#Unadjusted analysis -Temperature
##########################



temp.C1v4<-temp.C2v4<-temp.C1v2<-temp.C1v3<-temp.C3v4<-matrix(NA, nrow=3, ncol=7)
temp.I1v4<-temp.I2v4<-temp.I1v2<-temp.I1v3<-temp.I3v4<-matrix(NA, nrow=3, ncol=7)

for(j in 1:3){
  i<-j*7
  set.seed(12345)
  temp.C1v4[j,]<-quart_tmle(d=dc,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,4),cutoff=4)
  temp.C2v4[j,]<-quart_tmle(d=dc,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(2,4),cutoff=4)
  temp.C1v2[j,]<-quart_tmle(d=dc,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,2),cutoff=4)
  temp.C1v3[j,]<-quart_tmle(d=dc,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,3),cutoff=4)
  temp.C3v4[j,]<-quart_tmle(d=dc,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(3,4),cutoff=4)

  temp.I1v4[j,]<-quart_tmle(d=di,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,4),cutoff=4)
  temp.I2v4[j,]<-quart_tmle(d=di,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(2,4),cutoff=4)
  temp.I1v2[j,]<-quart_tmle(d=di,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,2),cutoff=4)
  temp.I1v3[j,]<-quart_tmle(d=di,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,3),cutoff=4)
  temp.I3v4[j,]<-quart_tmle(d=di,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(3,4),cutoff=4)
  
}

temp.C <- temp.I <- NULL
for(i in 1:3){
  temp.C<-rbind(temp.C,temp.C1v2[i,],temp.C1v3[i,],temp.C1v4[i,])
  temp.I<-rbind(temp.I,temp.I1v2[i,],temp.I1v3[i,],temp.I1v4[i,])
}
rownames(temp.C)<-rownames(temp.I)<-c("lag7 1v2","lag7 1v3","lag7 1v4",
                        "lag14 1v2","lag14 1v3","lag14 1v4",
                        "lag21 1v2","lag21 1v3","lag21 1v4")
colnames(temp.C)<-colnames(temp.I)<- c("PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
temp.C <- data.frame(exposure="Temperature",contrast=paste0(rownames(temp.C)," Control"), arm="Control", temp.C)
temp.I <- data.frame(exposure="Temperature",contrast=paste0(rownames(temp.I),arm=" Intervention"), arm="Intervention", temp.I)
temp.C$n <- 1:nrow(temp.C)
temp.I$n <- 1:nrow(temp.I)
  
  

  
  
  
##########################
#Unadjusted analysis -Heavy rain
##########################
  
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



##########################
#Unadjusted analysis
##########################

#merge suvery and weather
dim(survey)
dim(averain)
d<-merge(survey,averain, by="intdate", all.x = F, all.y = F) 
dim(d)

colnames(d)


#Set to date format for merge
PriorHeavyRain$intdate<-as.Date(paste(PriorHeavyRain$month,PriorHeavyRain$day,PriorHeavyRain$year, sep="/"),"%m/%d/%Y")
PriorHeavyRain<-subset(PriorHeavyRain, select=-c(year,month,day))

#merge suvery and weather
dim(survey)
dim(PriorHeavyRain)
dHR<-merge(survey,PriorHeavyRain, by="intdate", all.x = F, all.y = F)  
dHR<-merge(dHR,LT, by="intdate", all.x = T, all.y = F) 
dim(dHR)

rainmean.d <- merge(d, dHR,  by=c("intdate","id"), all.x = T, all.y = T)


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
  

dHR.C <- dHR %>% filter(wpi==0)
dHR.I <- dHR %>% filter(wpi==1)


Hrain.C<-Hrain.I<-matrix(NA, nrow=3, ncol=11)
Hrain.C[1,]<-HeavyRainTMLE(dat=dHR.C,i=7)
Hrain.C[2,]<-HeavyRainTMLE(dat=dHR.C,i=14)
Hrain.C[3,]<-HeavyRainTMLE(dat=dHR.C,i=21)

Hrain.I[1,]<-HeavyRainTMLE(dat=dHR.I,i=7)
Hrain.I[2,]<-HeavyRainTMLE(dat=dHR.I,i=14)
Hrain.I[3,]<-HeavyRainTMLE(dat=dHR.I,i=21)

colnames(Hrain.C)<-colnames(Hrain.I)<-c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.C
Hrain.I

Hrain.C <- data.frame(exposure="Rain", contrast=paste0(rownames(Hrain.C)," Control"), arm="Control", Hrain.C)
Hrain.I <- data.frame(exposure="Rain", contrast=paste0(rownames(Hrain.I),arm=" Intervention"), arm="Intervention", Hrain.I)
Hrain.C$n <- 1:nrow(Hrain.C)
Hrain.I$n <- 1:nrow(Hrain.I)


temp.C <- temp.C %>% subset(., select=c(exposure, contrast, arm, PR, ci.lb, ci.ub, n))
temp.I <- temp.I %>% subset(., select=c(exposure, contrast, arm, PR, ci.lb, ci.ub, n)) 
Hrain.C <- Hrain.C %>% subset(., select=c(exposure, contrast, arm, PR, ci.lb, ci.ub, n)) 
Hrain.I <- Hrain.I %>% subset(., select=c(exposure, contrast, arm, PR, ci.lb, ci.ub, n)) 


#Create plotting dataframe
tdf <- rbind(temp.C, temp.I, Hrain.C, Hrain.I)
tdf <- tdf %>% arrange(n)
tdf$contrast <- factor(as.character(tdf$contrast), levels=unique(tdf$contrast))



#Plot effect

  #plot
  RRplot<-ggplot(data=tdf) + 
    labs( x = "Treatment arm", y = "Risk ratio") +
    geom_hline(yintercept = 1) +
    scale_y_continuous(breaks=c(0.125,0.25,0.5,1,2,4,8), trans='log10') +
    coord_flip(ylim = c(0.25, 4,8)) +
    geom_pointrange( mapping=aes(x=contrast, y= PR, ymin=ci.lb, ymax=ci.ub , colour=arm)) +
    facet_wrap(~exposure, scales = "free_x") +
    theme(panel.border = element_blank(), 
          strip.background = element_blank(), legend.position = "none") + 
    ggtitle("Effect of temperature on diarrhea by intervention arm") +theme(legend.position="none")
  
  RRplot



