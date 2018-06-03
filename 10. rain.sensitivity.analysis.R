
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


#Set to date format for merge
averain$intdate<-as.Date(paste(averain$month,averain$day,averain$year, sep="/"),"%m/%d/%Y")
averain<-subset(averain, select=-c(year,month,day))


#Calculate heavy rainfall events
#80th percentile of rainfall
HeavyRainThres<-as.numeric(quantile(rain[rain[,4]!=0,4],probs = seq(0, 1, 0.1) ,na.rm=T)[9])
HeavyRainThres

#90th percentile of rainfall
HeavyRainThres90<-as.numeric(quantile(rain[rain[,4]!=0,4],probs = seq(0, 1, 0.1) ,na.rm=T)[10])
HeavyRainThres90

#70th percentile of rainfall
HeavyRainThres70<-as.numeric(quantile(rain[rain[,4]!=0,4],probs = seq(0, 1, 0.1) ,na.rm=T)[8])
HeavyRainThres70

#80th percentile of all days (including 0 rain)
HeavyRainThresAnyDays<-as.numeric(quantile(rain[,4],probs = seq(0, 1, 0.1) ,na.rm=T)[9])
HeavyRainThresAnyDays

HeavyRain<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres,1,0))
table(HeavyRain)

HeavyRain90<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres90,1,0))
table(HeavyRain90)

HeavyRain70<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThres70,1,0))
table(HeavyRain70)

HeavyRainAnyDays<-apply(rain[,-c(1:3)],2, function(x) ifelse(x>=HeavyRainThresAnyDays,1,0))
table(HeavyRainAnyDays)

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

PriorHeavyRainAnyDays<-matrix(data=NA,nrow=nrow(rain),ncol=28)
PriorHeavyRain.namesAnyDays<-rep(NA, 28)
for(i in 1:28){
  PriorHeavyRainAnyDays[,i]<-HeavyRain.window(i, 7, "rain", HeavyRainAnyDays)
  PriorHeavyRain.namesAnyDays[i]<- paste0("HeavyRainAny.lag",i)
}

PriorHeavyRain90<-matrix(data=NA,nrow=nrow(rain),ncol=28)
PriorHeavyRain.names90<-rep(NA, 28)
for(i in 1:28){
  PriorHeavyRain90[,i]<-HeavyRain.window(i, 7, "rain", HeavyRain90)
  PriorHeavyRain.names90[i]<- paste0("HeavyRain90.lag",i)
}


PriorHeavyRain70<-matrix(data=NA,nrow=nrow(rain),ncol=28)
PriorHeavyRain.names70<-rep(NA, 28)
for(i in 1:28){
  PriorHeavyRain70[,i]<-HeavyRain.window(i, 7, "rain", HeavyRain70)
  PriorHeavyRain.names70[i]<- paste0("HeavyRain70.lag",i)
}

PriorHeavyRain<-cbind(rain[,1:3],PriorHeavyRain, PriorHeavyRain90, PriorHeavyRainAnyDays, PriorHeavyRain70)
colnames(PriorHeavyRain)[4:ncol(PriorHeavyRain)]<-c(PriorHeavyRain.names,PriorHeavyRain.names90, PriorHeavyRain.namesAnyDays, PriorHeavyRain.names70)
PriorHeavyRain<-as.data.frame(PriorHeavyRain)
head(PriorHeavyRain, 30)


###################
#Survey data processing
###################

# Expand out factors into indicators 
dim(df)
df<-df[,-5] #remove duplicate "round" variable
head(df)

Y<-df$diar7d
h2s<-df$h2s
id<-df$hhid #SE clustered on hhid 
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
dHR<-merge(survey,PriorHeavyRain, by="intdate", all.x = F, all.y = F)  #Do I need to impute missing weather data 
dHR<-merge(dHR,LT, by="intdate", all.x = T, all.y = F)  
dim(dHR)



#Heavy rain TMLE function
HeavyRainTMLE<-function(dat=d,i,var="HeavyRain.lag",Wvars=NULL,sl.library=library){
  W<-NULL
  res<-matrix(NA, nrow=1, ncol=8)
  set.seed(12345)
 colnames(dat)[which(colnames(dat)==paste0(var,i))]="HeavyRain"

 df<-dat

  a=round(sum(df$HeavyRain==1 & df$Y==1, na.rm=T),0)
  b=round(sum(df$HeavyRain==1 & df$Y==0, na.rm=T),0)
  c=round(sum(df$HeavyRain==0 & df$Y==1, na.rm=T),0)
  D=round(sum(df$HeavyRain==0 & df$Y==0, na.rm=T),0)
  
if(!is.null(Wvars)){
  Wall<-subset(df, select=Wvars)
}
if(is.null(Wvars)){
   
      res<-as.numeric(washb_glm(Y=df$Y, tr=df$HeavyRain, id=df$id, contrast=c(0,1), family=binomial(link='log'))$TR)
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
  

#time lag sensitivity analysis, 90% heavy rain
Hrain.unadj.sense80<-matrix(NA, nrow=28, ncol=11)
for(i in 1:28){
Hrain.unadj.sense80[i,]<-HeavyRainTMLE(dat=dHR,i=i)
}
colnames(Hrain.unadj.sense80)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.unadj.sense80


Hrain.unadj.sense90<-matrix(NA, nrow=28, ncol=11)
for(i in 1:28){
Hrain.unadj.sense90[i,]<-HeavyRainTMLE(dat=dHR, var="HeavyRain90.lag" ,i=i)
}
colnames(Hrain.unadj.sense90)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.unadj.sense90

Hrain.unadj.senseAnyDays<-matrix(NA, nrow=28, ncol=11)
for(i in 1:28){
Hrain.unadj.senseAnyDays[i,]<-HeavyRainTMLE(dat=dHR, var="HeavyRainAny.lag" ,i=i)
}
colnames(Hrain.unadj.senseAnyDays)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.unadj.senseAnyDays

Hrain.unadj.sense70<-matrix(NA, nrow=28, ncol=11)
for(i in 1:28){
Hrain.unadj.sense70[i,]<-HeavyRainTMLE(dat=dHR, var="HeavyRain70.lag" ,i=i)
}
colnames(Hrain.unadj.sense70)<-c("a","b","c","d","PR,", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
Hrain.unadj.sense70


HrainSens<-cbind(Hrain.unadj.sense70[,5:7],
                 Hrain.unadj.sense80[,5:7],
                 Hrain.unadj.sense90[,5:7],
                 Hrain.unadj.senseAnyDays[,5:7])
colnames(HrainSens)<-c("PR70","ci.lb70","ci.ub70",
                       "PR80","ci.lb80","ci.ub80",
                       "PR90","ci.lb90","ci.ub90",
                       "PRAny","ci.lbAny","ci.ubAny")
save(HrainSens, file="Rain_sensitivity.Rdata")
HrainSens<-as.data.frame(HrainSens)
mean(HrainSens$PR70)
mean(HrainSens$PR80)
mean(HrainSens$PR90)
mean(HrainSens$PRAny)

mean(HrainSens$ci.lb70)
mean(HrainSens$ci.lb80)
mean(HrainSens$ci.lb90)
mean(HrainSens$ci.lbAny)

mean(HrainSens$ci.ub70)
mean(HrainSens$ci.ub80)
mean(HrainSens$ci.ub90)
mean(HrainSens$ci.ubAny)

library(xtable)
xtable(HrainSens)

tab<-matrix(NA, nrow=22, ncol=5)

tab[,1]<-c(0:21)
tab[,2]<-paste(sprintf("%1.2f",HrainSens[7:28,1])," (",sprintf("%1.2f",HrainSens[7:28,2]),", ",sprintf("%1.2f",HrainSens[7:28,3]),")", sep="")
tab[,3]<-paste(sprintf("%1.2f",HrainSens[7:28,4])," (",sprintf("%1.2f",HrainSens[7:28,5]),", ",sprintf("%1.2f",HrainSens[7:28,6]),")", sep="")
tab[,4]<-paste(sprintf("%1.2f",HrainSens[7:28,7])," (",sprintf("%1.2f",HrainSens[7:28,8]),", ",sprintf("%1.2f",HrainSens[7:28,9]),")", sep="")
tab[,5]<-paste(sprintf("%1.2f",HrainSens[7:28,10])," (",sprintf("%1.2f",HrainSens[7:28,11]),", ",sprintf("%1.2f",HrainSens[7:28,12]),")", sep="")

tab<-as.data.frame(tab)
rownames(tab)<-NULL

# Table functions
cleantable <- function(x,digits) {
 print( xtable(x,digits=digits),
        sanitize.text.function=function(y){y},
        floating=FALSE,
        include.rownames=FALSE,
        include.colnames=FALSE,
        only.contents=TRUE,
        hline.after=NULL
 )
}



cleantable(tab[,-5], 2)


#Create supplimentary figure
library("ggplot2")
library("ggthemes")
theme_set(theme_bw())

d<-HrainSens[7:28,]
d$x<-0:21

plotdf<-data.frame(Percentile = c(rep("70th percentile", 22), rep("80th percentile", 22), rep("90th percentile", 22)),
                   PR = rep(NA, 66),
                   ci.lb = rep(NA, 66),
                   ci.ub = rep(NA, 66))

plotdf[1:22, 2:4] <- d[,1:3]
plotdf[23:44, 2:4] <- d[,4:6]
plotdf[45:66, 2:4] <- d[,7:9]
plotdf$lag<-plotdf$lagjitter<-0:21
plotdf$lagjitter[1:22] <- plotdf$lagjitter[1:22] - 0.3
plotdf$lagjitter[45:66] <- plotdf$lagjitter[45:66] + 0.3


plotdf %>% group_by(Percentile) %>%
  summarize(mean( PR), var(PR))

#Add 7 days to the lag to match up plots with the 1,2,3 wk 
#terminology in manuscript
plotdf$lag <- plotdf$lag+7



setwd("C:/Users/andre/Documents/Trichy analysis/Figures and Tables/")
#pdf("HRain_sensitivity_figure.pdf",width=10,height=5)
png("HRain_sensitivity_figure.png",width=8.5,height=4.249999, units="in", res=800)

ggplot(plotdf, aes(x = lag)) +
    geom_pointrange(aes(y=PR, ymin = ci.lb, ymax = ci.ub), color="grey30") +
    geom_hline(yintercept = 1) +
    facet_grid( . ~ Percentile) +
    scale_y_log10(breaks=c(0.5, 0.75, 1, 1.5, 2), labels=c(0.5, 0.75, 1, 1.5, 2)) +
    labs(y = "Prevalence ratio",
         x = "Lag period (days)") +
    scale_x_continuous(breaks=c(7,14,21,28)) +
    theme(strip.background = element_blank())

dev.off()


