
rm(list=ls())

library(tidyverse)
library(washb)
library(caret)
library(haven)
library(lubridate)
library(mgcv)

source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")



setwd("C:/Users/andre/Dropbox/Trichy analysis/Data")
load("Cleaned data/prior_weather_dataset.Rdata")

df<-read_dta("trichy_fu_vilwat.dta")

bl<-read_dta("trichy_bl.dta")
bl <- bl %>% subset(., select=c(vilid, chhid))
bl$chhid <- substr(bl$chhid, 1,3)
bl <- bl %>% group_by(vilid) %>% slice(1) 


#extract covariates:
load("Cleaned data/analysis_datasets.Rdata")

covariates <- d %>% group_by(vilid) %>% slice(1) %>% subset(., select=c(vilid, wpi,villOD )) %>% ungroup 

###################
#village Survey data processing
###################


#clean and select needed variables
df <- df %>% subset(., select=c(sid, svy, vdatecol, vdatearr, v10)) %>% rename(Y=v10, individ=sid) %>% filter(!is.na(Y))
head(df)

#get village id
df$chhid <- substr(df$individ, 1,3)
df$individ <- as.numeric(as.factor(df$individ))

#Village names Melpathu and Thudaiyur are the same village ID/intervention village with the same coordinates, so combine here
df$chhid[df$chhid=="TDR"] <- "MLP"

df <- left_join(df, bl, by="chhid")
table(df$vilid)
table(is.na(df$individ))

#Mean and SD of water sources per village
df %>% group_by(chhid,individ) %>% slice(1) %>% ungroup() %>% summarize(n())
df %>% group_by(chhid,individ) %>% slice(1) %>% group_by(chhid) %>% summarize(N=n()) %>% 
        ungroup() %>% summarize(mean(N), sd(N), min(N), max(N))
df %>% group_by(chhid,individ) %>% summarize(N=n()) %>% 
  ungroup() %>% summarize(mean(N), sd(N), min(N), max(N))


#fix errors in date input
df$vdatecol[31] <- "03/06/08"
df$vdatecol[87] <- "20/08/08"
df$vdatecol[84] <- "17/05/08"

df$intdate <- as.Date(df$vdatecol,format='%d/%m/%y')


#merge weather with survey
d <- merge(df, weather, by=c("intdate"), all.x=T, all.y=F) %>% arrange(vilid, intdate)


d <- d %>% filter(!is.na(vilid))

#Get study week
d$stdywk <- week(d$intdate)



#merge in covariates
d <- left_join(d, covariates, by="vilid")
table(d$wpi)

Wvars<-c("wpi","villOD")

d$id <- d$individ

# N's and prevalences by group
d %>% filter(!is.na(HeavyRain.lag1)) %>% group_by(HeavyRain.lag1) %>% summarise(N=n(), mean=mean(Y))
d %>% filter(!is.na(HeavyRain.lag7)) %>% group_by(HeavyRain.lag7) %>% summarise(N=n(), mean=mean(Y))
d %>% filter(!is.na(HeavyRain.lag14)) %>% group_by(HeavyRain.lag14) %>% summarise(N=n(), mean=mean(Y))

d %>% filter(!is.na(HeavyRain.lag1)) %>% group_by(LT1_T, HeavyRain.lag1) %>% summarise(N=n(), mean=mean(Y))
d %>% filter(!is.na(HeavyRain.lag7)) %>% group_by(LT8_T, HeavyRain.lag7) %>% summarise(N=n(), mean=mean(Y))
d %>% filter(!is.na(HeavyRain.lag14)) %>% group_by(LT15_T, HeavyRain.lag14) %>% summarise(N=n(), mean=mean(Y))



#Summary statistics
summary_fun <- function(d, A, weather, strat=NULL,strat_weathervar=NULL){
  colnames(d)[colnames(d)==A] <- "Avar"
  colnames(d)[colnames(d)==weather] <- "weathervar"
  if(is.null(strat)){
    res<-d %>% group_by(Avar) %>% filter(!is.na(Y)&!is.na(Avar)) %>% summarize(N=n(), count = sum(Y), prev=round(mean(Y)*100,2), weather_ave=round(mean(weathervar),2)) %>%
      as.data.frame()
    print(res)
    
  }else{
    colnames(d)[colnames(d)==strat] <- "strat"
    colnames(d)[colnames(d)==strat_weathervar] <- "strat_weathervar"
    res<-d %>% group_by(strat,Avar) %>% filter(!is.na(Y)&!is.na(Avar)&!is.na(strat)) %>% summarize(N=n(), count = sum(Y), 
          prev=round(mean(Y)*100,2), weather_ave=round(mean(weathervar),2)) %>%
      as.data.frame()
    res2<-d %>% group_by(strat) %>% filter(!is.na(Y)&!is.na(Avar)&!is.na(strat)) %>% summarize(strat_ave=round(mean(strat_weathervar),2)) %>%
      as.data.frame()
    print(res)
    print(res2)
  }
}


summary_fun(d, "HeavyRain.lag1","rain.ave7.lag1")
summary_fun(d, "HeavyRain.lag7","rain.ave7.lag8")
summary_fun(d, "HeavyRain.lag14","rain.ave7.lag15")

summary_fun(d, "HeavyRain.lag1","rain.ave7.lag1", "LT1_T","LT1")
summary_fun(d, "HeavyRain.lag7","rain.ave7.lag8", "LT8_T","LT8")
summary_fun(d, "HeavyRain.lag14","rain.ave7.lag15", "LT15_T","LT15")

#Prevalence by sampling round
d %>% group_by(svy, vilid) %>% summarise(mn=mean(Y)) %>% ungroup %>% summarise(min(mn), mean(mn)) 
d %>% group_by(svy) %>% summarise(mn=mean(Y), N=n()) %>% ungroup %>% summarise(min(mn), mean(mn)) 
d %>% group_by(svy) %>% summarise(mn=mean(Y), N=n()) 


##########################
#Unadjusted analysis
##########################

h2s.HR1<-trichy_gamm(d, A="HeavyRain.lag1")
h2s.HR1

h2s.HR2<-trichy_gamm(d, A="HeavyRain.lag7")
h2s.HR2

h2s.HR3<-trichy_gamm(d, A="HeavyRain.lag14")
h2s.HR3

h2s.HR1_strat<-trichy_gamm(d, A="HeavyRain.lag1", strat="LT1_T")
h2s.HR1_strat

h2s.HR2_strat<-trichy_gamm(d, A="HeavyRain.lag7", strat="LT8_T")
h2s.HR2_strat

h2s.HR3_strat<-trichy_gamm(d, A="HeavyRain.lag14", strat="LT15_T")
h2s.HR3_strat


#Adjusted
h2s.HR1_adj<-trichy_gamm(d, Y="H2S", A="HeavyRain.lag1", weathervar="temp.ave7.lag1", Wvars = Wvars)
h2s.HR1_adj

h2s.HR2_adj<-trichy_gamm(d, Y="H2S", A="HeavyRain.lag7", weathervar="temp.ave7.lag7", Wvars = Wvars)
h2s.HR2_adj

h2s.HR3_adj<-trichy_gamm(d, Y="H2S", A="HeavyRain.lag14", weathervar="temp.ave7.lag14", Wvars = Wvars)
h2s.HR3_adj

#stratified
h2s.HR1_strat_adj <- h2s.HR2_strat_adj <- h2s.HR3_strat_adj <- NULL

h2s.HR1_strat_adj<-trichy_gamm(d, Y="H2S", A="HeavyRain.lag1", strat="LT1_T", weathervar="temp.ave7.lag1", Wvars = Wvars)
h2s.HR1_strat_adj

h2s.HR2_strat_adj<-trichy_gamm(d, Y="H2S", A="HeavyRain.lag7", strat="LT8_T", weathervar="temp.ave7.lag7", Wvars = Wvars)
h2s.HR2_strat_adj

h2s.HR3_strat_adj<-trichy_gamm(d, Y="H2S", A="HeavyRain.lag14", strat="LT15_T", weathervar="temp.ave7.lag14", Wvars = Wvars)
# No covariates selected in adjusted analysis so converges to unadjusted:
h2s.HR3_strat_adj<-trichy_gamm(d, Y="H2S", A="HeavyRain.lag14", strat="LT15_T", weathervar="temp.ave7.lag14", Wvars = NULL)
h2s.HR3_strat_adj$adjusted <- "Y"
h2s.HR3_strat_adj




vilH2S_resdf <- rbind( 
 h2s.HR1, h2s.HR2, h2s.HR3,
 h2s.HR1_strat, h2s.HR2_strat, h2s.HR3_strat)

vilH2S_resdf_adj <- rbind( 
  h2s.HR1_adj, h2s.HR2_adj, h2s.HR3_adj,
  h2s.HR1_strat_adj, h2s.HR2_strat_adj, h2s.HR3_strat_adj)


save(vilH2S_resdf, vilH2S_resdf_adj, file="C:/Users/andre/Dropbox/Trichy analysis/Results/villageh2s.HRmean-gamm-results.Rdata")

#Print results in table format
vilH2S_resdf$PR <- round(vilH2S_resdf$PR, 2)
vilH2S_resdf$ci.lb <- round(vilH2S_resdf$ci.lb, 2)
vilH2S_resdf$ci.ub <- round(vilH2S_resdf$ci.ub, 2)
vilH2S_resdf$y <- paste0(vilH2S_resdf$PR," (", vilH2S_resdf$ci.lb,", ", vilH2S_resdf$ci.ub,")")

df <- vilH2S_resdf[,c(2,8)]
df
noquote(df[,2])



vilH2S_resdf_adj$PR <- round(vilH2S_resdf_adj$PR, 2)
vilH2S_resdf_adj$ci.lb <- round(vilH2S_resdf_adj$ci.lb, 2)
vilH2S_resdf_adj$ci.ub <- round(vilH2S_resdf_adj$ci.ub, 2)
vilH2S_resdf_adj$y <- paste0(vilH2S_resdf_adj$PR," (", vilH2S_resdf_adj$ci.lb,", ", vilH2S_resdf_adj$ci.ub,")")

df <- vilH2S_resdf_adj[,c(2,8)]
df
noquote(df[,2])

