


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


source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")

try(setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data"))

load("analysis_datasets.Rdata")

####################
#Set ID to village ID
####################

d$id <- d$vilid

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
wpi_lrtest(d, A="HeavyRain.lag7")
wpi_lrtest(d, A="HeavyRain.lag14")
wpi_lrtest(d, A="HeavyRain.lag21")


    fit1 <- glm(Y ~ temp.ave7.lag7Q*wpi + temp.ave7.lag14Q*wpi + temp.ave7.lag21Q*wpi, data = d, family = "binomial")
    fit0 <- glm(Y ~ temp.ave7.lag7Q + temp.ave7.lag14Q + temp.ave7.lag21Q, data = d, family = "binomial")
    LRp <- lrtest(fit1, fit0)[2, 5]
    LRp
    
    fit1 <- glm(Y ~ HeavyRain.lag7*wpi + HeavyRain.lag14*wpi + HeavyRain.lag21*wpi, data = d, family = "binomial")
    fit0 <- glm(Y ~ HeavyRain.lag7 + HeavyRain.lag14 + HeavyRain.lag21, data = d, family = "binomial")
    LRp <- lrtest(fit1, fit0)[2, 5]
   LRp


#-------------------------------------------------
#Subset dataset into control and intervention arms
#-------------------------------------------------
dc <- d %>% filter(wpi==0)
di <- d %>% filter(wpi==1)   
   
   
   
   
############################################################
#Control group results
############################################################

d <- dc


####################
#Unadjusted MH tests
####################

temp.unadj1v4<-temp.unadj1v2<-temp.unadj1v3<-matrix(NA, nrow=3, ncol=7)
for(j in 1:3){
  i<-j*7
  set.seed(12345)
  temp.unadj1v2[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("tempQ",i),q=c(1,2),cutoff=4)
  temp.unadj1v3[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("tempQ",i),q=c(1,3),cutoff=4)
  temp.unadj1v4[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("tempQ",i),q=c(1,4),cutoff=4)
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


Hrain<-matrix(NA, nrow=3, ncol=11)
Hrain[1,]<-HeavyRainTMLE(dat=d,i=7)
Hrain[2,]<-HeavyRainTMLE(dat=d,i=14)
Hrain[3,]<-HeavyRainTMLE(dat=d,i=21)

colnames(Hrain) <- c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")



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


ctrl_subgroup_res <- list(temp.unadj, Hrain, Hrain.unadjLT1,Hrain.unadjLT2, Hrain.unadjLT3, H2S.Hrain.unadj, H2s.unadjLT1, H2s.unadjLT2, H2s.unadjLT3)





############################################################
#Intervention group results
############################################################

d <- di


####################
#Unadjusted MH tests
####################

temp.unadj1v4<-temp.unadj1v2<-temp.unadj1v3<-matrix(NA, nrow=3, ncol=7)
for(j in 1:3){
  i<-j*7
  set.seed(12345)
  temp.unadj1v2[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("tempQ",i),q=c(1,2),cutoff=4)
  temp.unadj1v3[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("tempQ",i),q=c(1,3),cutoff=4)
  temp.unadj1v4[j,]<-quart_tmle(d=d,Wvars=NULL,weather=paste0("tempQ",i),q=c(1,4),cutoff=4)
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


Hrain<-matrix(NA, nrow=3, ncol=11)
Hrain[1,]<-HeavyRainTMLE(dat=d,i=7)
Hrain[2,]<-HeavyRainTMLE(dat=d,i=14)
Hrain[3,]<-HeavyRainTMLE(dat=d,i=21)

colnames(Hrain) <- c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")



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





trt_subgroup_res <- list(temp.unadj, Hrain, Hrain.unadjLT1,Hrain.unadjLT2, Hrain.unadjLT3, H2S.Hrain.unadj, H2s.unadjLT1, H2s.unadjLT2, H2s.unadjLT3)



# Save results
save(ctrl_subgroup_res, trt_subgroup_res,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/treatment_subgroup_results.Rdata")




  
  
      