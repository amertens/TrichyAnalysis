
rm(list=ls())

library(dplyr)
library(foreign)
library(zoo)
library(usdm)

###Generating lagged weather variables
setwd("C:/Users/andre/Dropbox/Trichy analysis/Data")
df<-read.dta("Trichy_Weather_Dec07-Apr09_formatted.dta")

head(df)

#Order data
df<-arrange(df, year, month, day) %>% 
    subset(., select=c(year, month, day, maxtemp, mintemp, rain))


#Generate average temp
df$avetemp<-(df$maxtemp+df$mintemp)/2


# examine correlation between temperature and rainfall
cor.test(df$avetemp, df$rain)
vif(data.frame(df$avetemp, df$rain))


#Generate lagged ave temp
lag.n <- function(df,name,var, n) {
    varname <- paste(name, n , sep=".")
    df[[varname]] <- with(df,  lag(var, n))
    return(df)
}


dayslag<-100

for(i in 1:dayslag){
  df<-lag.n(df=df,name="at",var=df$avetemp,n=i)
}


#Generate lagged min and max temp
for(i in 1:dayslag){
  df<-lag.n(df=df,name="maxt",var=df$maxtemp,n=i)
}
for(i in 1:dayslag){
  df<-lag.n(df=df,name="mint",var=df$mintemp,n=i)
}



#Generate lagged rain
for(i in 1:dayslag){
  df<-lag.n(df=df,name="rain",var=df$rain,n=i)
}






#save data
weather<-df

save(weather, file="C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data/LaggedWeather.Rdata")




#--------------------------------------
#  Clean weather data into exposure datasets
#--------------------------------------




###################
#Temperature data processing
###################
head(weather)


#Get quartiles cutpoints of overall weekly mean
weekly_temp<-rollmean(weather$avetemp,7,fill=NA, align="right")
tempQ <- as.numeric(quantile(weekly_temp, na.rm=T)[2:4])

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
avetemp<-subset(avetemp, select=c(intdate, temp.ave7.lag1, temp.ave7.lag7, temp.ave7.lag14, temp.ave7.lag21))

#Quartile temperatures
avetemp$tempQ1<-cut(avetemp$temp.ave7.lag1, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
avetemp$tempQ7<-cut(avetemp$temp.ave7.lag7, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
avetemp$tempQ14<-cut(avetemp$temp.ave7.lag14, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
avetemp$tempQ21<-cut(avetemp$temp.ave7.lag21, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))


# Minimum temp
min.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) min(x)))
  return(windowmean)
}

mintemp<-matrix(data=NA,nrow=nrow(temp),ncol=21)
min.names<-rep(NA, 21)
for(i in 1:21){
  mintemp[,i]<-min.window(i, 7, "at", temp)
  min.names[i]<- paste0("temp.min7.lag",i)
}
mintemp<-cbind(temp[,1:3],mintemp)
colnames(mintemp)[4:ncol(mintemp)]<-min.names
mintemp<-as.data.frame(mintemp)


#Set to date format for merge
mintemp$intdate<-as.Date(paste(mintemp$month,mintemp$day,mintemp$year, sep="/"),"%m/%d/%Y")
mintemp<-subset(mintemp, select=c(intdate, temp.min7.lag7, temp.min7.lag14, temp.min7.lag21))

#Quartile temperatures
mintemp$mintempQ1<-cut(mintemp$temp.min7.lag1, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
mintemp$mintempQ7<-cut(mintemp$temp.min7.lag7, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
mintemp$mintempQ14<-cut(mintemp$temp.min7.lag14, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
mintemp$mintempQ21<-cut(mintemp$temp.min7.lag21, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))

# Maximum temp
max.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) max(x)))
  return(windowmean)
}

maxtemp<-matrix(data=NA,nrow=nrow(temp),ncol=21)
max.names<-rep(NA, 21)
for(i in 1:21){
  maxtemp[,i]<-max.window(i, 7, "at", temp)
  max.names[i]<- paste0("temp.max7.lag",i)
}
maxtemp<-cbind(temp[,1:3],maxtemp)
colnames(maxtemp)[4:ncol(maxtemp)]<-max.names
maxtemp<-as.data.frame(maxtemp)


#Set to date format for merge
maxtemp$intdate<-as.Date(paste(maxtemp$month,maxtemp$day,maxtemp$year, sep="/"),"%m/%d/%Y")
maxtemp<-subset(maxtemp, select=c(intdate, temp.max7.lag7, temp.max7.lag14, temp.max7.lag21))

#Quartile temperatures
maxtemp$maxtempQ1<-cut(maxtemp$temp.max7.lag1, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
maxtemp$maxtempQ7<-cut(maxtemp$temp.max7.lag7, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
maxtemp$maxtempQ14<-cut(maxtemp$temp.max7.lag14, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
maxtemp$maxtempQ21<-cut(maxtemp$temp.max7.lag21, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))


###################
#Rain data processing
###################
head(weather)

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
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) mean(x))) #NOTE: This does not computes averages when there are NA's
    return(windowmean)
}

averain<-matrix(data=NA,nrow=nrow(rain),ncol=21)
ave.names<-rep(NA, 28)
for(i in 1:28){
  averain[,i]<-ave.window(i, 7, "rain", rain)
  ave.names[i]<- paste0("rain.ave7.lag",i)
}
averain<-cbind(rain[,1:3],averain)
colnames(averain)[4:ncol(averain)]<-ave.names
averain<-as.data.frame(averain)

#Create long term rain average

#60 day mean
bimonth_rain<-rollmean(weather$rain,60,fill=NA, align="right")
LTrainQ <- as.numeric(quantile(bimonth_rain,probs = seq(0, 1, 1/3), na.rm=T)[2:3])

#set i to 7 days after the associated 7-day lagged variable
averain$LT1<-ave.window(8, 60, "rain", rain)
averain$LT8<-ave.window(15, 60, "rain", rain)
averain$LT15<-ave.window(22, 60, "rain", rain)
averain$LT22<-ave.window(29, 60, "rain", rain)

averain$LT1_T<-cut(averain$LT1, breaks = c(0,LTrainQ,999), labels=c("T1","T2","T3")) 
averain$LT8_T<-cut(averain$LT8, breaks = c(0,LTrainQ,999), labels=c("T1","T2","T3")) 
averain$LT15_T<-cut(averain$LT15, breaks = c(0,LTrainQ,999), labels=c("T1","T2","T3")) 
averain$LT22_T<-cut(averain$LT22, breaks = c(0,LTrainQ,999), labels=c("T1","T2","T3")) 

table(averain$LT8_T)
table(averain$LT15_T)
table(averain$LT22_T)


#Set date to date format for merge
averain$intdate<-as.Date(paste(averain$month,averain$day,averain$year, sep="/"),"%m/%d/%Y")
LT <- subset(averain, select=c(intdate, rain.ave7.lag1, rain.ave7.lag7, rain.ave7.lag14, rain.ave7.lag21, LT1, LT8, LT15, LT22, LT1_T, LT8_T, LT15_T, LT22_T))


#Calculate heavy rainfall events
#80th percentile of rainfall on rainy days
HeavyRainThres<-as.numeric(quantile(rain[rain[,4]!=0,4],probs = seq(0, 1, 0.1) ,na.rm=T)[9])

#70th percentile of rainfall on rainy days
HeavyRainThres70<-as.numeric(quantile(rain[rain[,4]!=0,4],probs = seq(0, 1, 0.1) ,na.rm=T)[8])

#90th percentile of rainfall on rainy days
HeavyRainThres90<-as.numeric(quantile(rain[rain[,4]!=0,4],probs = seq(0, 1, 0.1) ,na.rm=T)[10])


#90th percentile of rainfall on any days
HeavyRainThres90a<-as.numeric(quantile(rain,probs = seq(0, 1, 0.1) ,na.rm=T)[10])


HeavyRain<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres,1,0))
table(HeavyRain)

HeavyRain70<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres70,1,0))
table(HeavyRain70)

HeavyRain90<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres90,1,0))
table(HeavyRain90)

HeavyRain90a<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres90a,1,0))
table(HeavyRain90a)

 

HeavyRain.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) ifelse(sum(x)>0,1,0)))
  return(windowmean)
}


PriorHeavyRain<-matrix(data=NA,nrow=nrow(rain),ncol=28)
PriorHeavyRain.names<-rep(NA, 28)
for(i in 1:28){
  PriorHeavyRain[,i]<-HeavyRain.window(i, 7, "rain", HeavyRain)
  PriorHeavyRain.names[i]<- paste0("HeavyRain.lag",i)
}
PriorHeavyRain<-cbind(rain[,1:3],PriorHeavyRain)
colnames(PriorHeavyRain)[4:ncol(PriorHeavyRain)]<-PriorHeavyRain.names
PriorHeavyRain<-as.data.frame(PriorHeavyRain)
head(PriorHeavyRain,30)


#Set to date format for merge
PriorHeavyRain$intdate<-as.Date(paste(PriorHeavyRain$month,PriorHeavyRain$day,PriorHeavyRain$year, sep="/"),"%m/%d/%Y")
PriorHeavyRain_sens<-subset(PriorHeavyRain, select=-c(year, month, day)) #Keep all lags for sensitivity analysis
PriorHeavyRain<-subset(PriorHeavyRain, select=c(intdate, HeavyRain.lag1, HeavyRain.lag7,HeavyRain.lag14,HeavyRain.lag21))

#Sensitivity rain dataframes
PriorHeavyRain70<-matrix(data=NA,nrow=nrow(rain),ncol=28)
PriorHeavyRain70.names<-rep(NA, 28)
for(i in 1:28){
  PriorHeavyRain70[,i]<-HeavyRain.window(i, 7, "rain", HeavyRain70)
  PriorHeavyRain70.names[i]<- paste0("HeavyRain70.lag",i)
}
PriorHeavyRain70<-cbind(rain[,1:3],PriorHeavyRain70)
colnames(PriorHeavyRain70)[4:ncol(PriorHeavyRain70)]<-PriorHeavyRain70.names
PriorHeavyRain70<-as.data.frame(PriorHeavyRain70)
PriorHeavyRain70$intdate<-as.Date(paste(PriorHeavyRain70$month,PriorHeavyRain70$day,PriorHeavyRain70$year, sep="/"),"%m/%d/%Y")
PriorHeavyRain70<-subset(PriorHeavyRain70, select=-c(year, month, day))
                            
PriorHeavyRain90<-matrix(data=NA,nrow=nrow(rain),ncol=28)
PriorHeavyRain90.names<-rep(NA, 28)
for(i in 1:28){
  PriorHeavyRain90[,i]<-HeavyRain.window(i, 7, "rain", HeavyRain90)
  PriorHeavyRain90.names[i]<- paste0("HeavyRain90.lag",i)
}
PriorHeavyRain90<-cbind(rain[,1:3],PriorHeavyRain90)
colnames(PriorHeavyRain90)[4:ncol(PriorHeavyRain90)]<-PriorHeavyRain90.names
PriorHeavyRain90<-as.data.frame(PriorHeavyRain90)
PriorHeavyRain90$intdate<-as.Date(paste(PriorHeavyRain90$month,PriorHeavyRain90$day,PriorHeavyRain90$year, sep="/"),"%m/%d/%Y")
PriorHeavyRain90<-subset(PriorHeavyRain90, select=-c(year, month, day))
                            
PriorHeavyRain90a<-matrix(data=NA,nrow=nrow(rain),ncol=28)
PriorHeavyRain90a.names<-rep(NA, 28)
for(i in 1:28){
  PriorHeavyRain90a[,i]<-HeavyRain.window(i, 7, "rain", HeavyRain90a)
  PriorHeavyRain90a.names[i]<- paste0("HeavyRain90a.lag",i)
}
PriorHeavyRain90a<-cbind(rain[,1:3],PriorHeavyRain90a)
colnames(PriorHeavyRain90a)[4:ncol(PriorHeavyRain90a)]<-PriorHeavyRain90a.names
PriorHeavyRain90a<-as.data.frame(PriorHeavyRain90a)
PriorHeavyRain90a$intdate<-as.Date(paste(PriorHeavyRain90a$month,PriorHeavyRain90a$day,PriorHeavyRain90a$year, sep="/"),"%m/%d/%Y")
PriorHeavyRain90a<-subset(PriorHeavyRain90a, select=-c(year, month, day))

#--------------------------------------
#Load and merge in survey data
#--------------------------------------


 load("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data/survey_dataset.Rdata")

#merge survey and temp
dim(survey)
dim(avetemp)
d<-merge(survey, avetemp, by="intdate", all.x = F, all.y = F) 
dim(d)
d<-merge(d, mintemp, by="intdate", all.x = F, all.y = F) 
d<-merge(d, maxtemp, by="intdate", all.x = F, all.y = F) 
dim(d)
colnames(d)

#merge suvery and rain
dim(survey)
dim(PriorHeavyRain)
d<-merge(d,PriorHeavyRain, by="intdate", all.x = F, all.y = F)  
d<-merge(d,PriorHeavyRain70, by="intdate", all.x = F, all.y = F)  
d<-merge(d,PriorHeavyRain90, by="intdate", all.x = F, all.y = F)  
d<-merge(d,PriorHeavyRain90a, by="intdate", all.x = F, all.y = F)  
d<-merge(d,LT, by="intdate", all.x = T, all.y = F) 
dim(d)

colnames(d)



#Make HR sensitivity analysis dataframe
d_HRsens<-merge(survey, PriorHeavyRain_sens, by="intdate", all.x = F, all.y = F) 
d_HRsens<-merge(d_HRsens, PriorHeavyRain70, by="intdate", all.x = F, all.y = F)  
d_HRsens<-merge(d_HRsens, PriorHeavyRain90, by="intdate", all.x = F, all.y = F)  
d_HRsens<-merge(d_HRsens, PriorHeavyRain90a, by="intdate", all.x = F, all.y = F)  
d_HRsens<-merge(d_HRsens, LT, by="intdate", all.x = T, all.y = F) 
dim(d_HRsens)


save(d, d_HRsens, tempQ, LT, Wvars, file="C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data/analysis_datasets.Rdata")



#Make weather dataset to merge with village H2S measurements

weather <- merge(PriorHeavyRain, avetemp, by="intdate")
weather <- merge(weather, LT, by="intdate")
save(weather, file="C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data/prior_weather_dataset.Rdata")



#Get descriptive statistics of the mean weather by quartile for the supplimentary tables
d %>% group_by(tempQ7) %>% filter(!is.na(Y)) %>% summarize(round(mean(temp.ave7.lag7),2))
d %>% group_by(tempQ14) %>% filter(!is.na(Y)) %>% summarize(round(mean(temp.ave7.lag14),2))
d %>% group_by(tempQ21) %>% filter(!is.na(Y)) %>% summarize(round(mean(temp.ave7.lag21),2))

d %>% group_by(tempQ1) %>% filter(!is.na(H2S)) %>% summarize(round(mean(temp.ave7.lag1),2))
d %>% group_by(tempQ7) %>% filter(!is.na(H2S)) %>% summarize(round(mean(temp.ave7.lag7),2))
d %>% group_by(tempQ14) %>% filter(!is.na(H2S)) %>% summarize(round(mean(temp.ave7.lag14),2))
