


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





################################
#Calculate LR p-value of sanitation
#Coverage as an interaction term
################################

summary(glm(d$Y~ d$sanlat))

#Effect modification?
summary(glm(d$Y~ d$sanlat * d$temp.ave7.lag7, family="binomial"))
summary(glm(d$Y~ d$sanlat * d$temp.ave7.lag14, family="binomial"))
summary(glm(d$Y~ d$sanlat * d$temp.ave7.lag21, family="binomial"))
summary(glm(d$Y~ d$sanlat * d$rain.ave7.lag7, family="binomial"))
summary(glm(d$Y~ d$sanlat * d$rain.ave7.lag14, family="binomial"))
summary(glm(d$Y~ d$sanlat * d$rain.ave7.lag21, family="binomial"))


 sanlat_lrtest<- function (dat, A, Y="Y", family = "binomial", print=F){
  require(lmtest)
   dat$A <- dat[,A]

    fit1 <- glm(Y ~ A*sanlat, data = dat, family = family)   #Note, update this with SE clustering
    
    if(print==T){print(summary(fit1))}
    
    fit0 <- glm(Y ~ A, data = dat, family = family)
    LRp <- lrtest(fit1, fit0)[2, 5]

    return(LRp)
}

#Continious exposures
# sanlat_lrtest(d, A="temp.ave7.lag7")
# sanlat_lrtest(d, A="temp.ave7.lag14")
# sanlat_lrtest(d, A="temp.ave7.lag21")
# 
# sanlat_lrtest(d, A="rain.ave7.lag7")
# sanlat_lrtest(d, A="rain.ave7.lag14")
# sanlat_lrtest(d, A="rain.ave7.lag21")


#Categorized temperature exposures
d$temp.ave7.lag7Q<-cut(d$temp.ave7.lag7, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
d$temp.ave7.lag14Q<-cut(d$temp.ave7.lag14, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))
d$temp.ave7.lag21Q<-cut(d$temp.ave7.lag21, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))

sanlat_lrtest(d, A="temp.ave7.lag7Q")
sanlat_lrtest(d, A="temp.ave7.lag14Q")
sanlat_lrtest(d, A="temp.ave7.lag21Q")

#Heavy rain
sanlat_lrtest(d, A="HeavyRain.lag7")
sanlat_lrtest(d, A="HeavyRain.lag14")
sanlat_lrtest(d, A="HeavyRain.lag21")


    fit1 <- glm(Y ~ temp.ave7.lag7Q*sanlat + temp.ave7.lag14Q*sanlat + temp.ave7.lag21Q*sanlat, data = d, family = "binomial")
    fit0 <- glm(Y ~ temp.ave7.lag7Q + temp.ave7.lag14Q + temp.ave7.lag21Q, data = d, family = "binomial")
    LRp <- lrtest(fit1, fit0)[2, 5]
    LRp
    
    fit1 <- glm(Y ~ HeavyRain.lag7*sanlat + HeavyRain.lag14*sanlat + HeavyRain.lag21*sanlat, data = d, family = "binomial")
    fit0 <- glm(Y ~ HeavyRain.lag7 + HeavyRain.lag14 + HeavyRain.lag21, data = d, family = "binomial")
    LRp <- lrtest(fit1, fit0)[2, 5]
   LRp


   
   
############################################################
# GAMM analysis
############################################################
   
   m <- gamm(Y ~ tempQ7*wpi , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log'))
  res<-summary(m$gam)
  
  m0 <- gamm(Y ~ tempQ7 , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log'))
  res0<-summary(m0$gam)
anova(m$gam, m0$gam)

print(x, digits = max(3, getOption("digits") - 3),...)
"Chisq", "F" or "Cp"

a <- m$lme
b <- m0$lme
anova.lme(a,b)

gamm_anova<- function (dat, A, V, print=T){
   dat$A <- dat[,A]
   dat$V <- dat[,V]

    fit1 <- gamm(Y ~ A*V , data=dat, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log'))
    
    if(print==T){print(summary(fit1$gam))}
    
    fit0 <-  gamm(Y ~ A , data=dat, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log'))
    anova.res <- anova.gam(fit1$gam, fit0$gam)
    if(print==T){print(anova.res)}

    return(anova.res)
}


temp1_wpi<-gamm_anova(d, A="tempQ7", V="wpi")
temp2_wpi<-gamm_anova(d, A="tempQ14", V="wpi")
temp3_wpi<-gamm_anova(d, A="tempQ21", V="wpi")

HR1_wpi<-gamm_anova(d, A="HeavyRain.lag7", V="wpi")
HR2_wpi<-gamm_anova(d, A="HeavyRain.lag14", V="wpi")
HR3_wpi<-gamm_anova(d, A="HeavyRain.lag21", V="wpi")


temp1_sanlat<-gamm_anova(d, A="tempQ7", V="sanlat")
temp2_sanlat<-gamm_anova(d, A="tempQ14", V="sanlat")
temp3_sanlat<-gamm_anova(d, A="tempQ21", V="sanlat")

HR1_sanlat<-gamm_anova(d, A="HeavyRain.lag7", V="sanlat")
HR2_sanlat<-gamm_anova(d, A="HeavyRain.lag14", V="sanlat")
HR3_sanlat<-gamm_anova(d, A="HeavyRain.lag21", V="sanlat")



#test	what sort of test to perform for a multi-model call. One of "Chisq", "F" or "LRT".


   m <- gam(Y ~ tempQ7*wpi , data=d, random = list(vilid=~1), family=binomial(link='log'))

  m0 <- gam(Y ~ tempQ7 , data=d, random = list(vilid=~1),  family=binomial(link='log'))
anova.gam(m, m0,  test = "Chisq")
   

gamm_anova<- function (dat, A, V, print=T){
   dat$A <- dat[,A]
   dat$V <- dat[,V]

    fit1 <- gamm(Y ~ A*V , data=dat, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log'))
    
    if(print==T){print(summary(fit1$gam))}
    
    fit0 <-  gamm(Y ~ A , data=dat, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log'))
    anova.res <- anova.gam(fit1$gam, fit0$gam)
    if(print==T){print(anova.res)}

    return(anova.res)
}



   m <- gamm(Y ~ HeavyRain.lag14*sanlat , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log'))

  m0 <- gamm(Y ~ HeavyRain.lag14 , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log'))

  
  glmmPQL(Y ~ HeavyRain.lag14*sanlat , data=d, random = ~1|vilid, correlation = corAR1(form = ~ stdywk|individ), family=binomial)
  
a <- logLik(m$lme)[1] 
b <- logLik(m0$lme)[1]

DF <- attr(logLik(m$lme), "df") -  attr(logLik(m0$lme), "df")
dcxfv ]

AIC(m$lme) - AIC(m0$lme)


compareML( m0$gam, m$gam)

#-------------------------------------------------
#Subset dataset into coverad and not-covered 
#-------------------------------------------------
dc <- d %>% filter(sanlat==0)
di <- d %>% filter(sanlat==1)   
   
   
   
   
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


HR.strat.rainmean7<- d %>% group_by(LT8_T, HeavyRain.lag7) %>% summarize(mean.rainfall=mean(rain.ave7.lag7, na.rm=T), mean.h2s=mean(H2S, na.rm=T)) %>% as.data.frame()
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




  
  
      