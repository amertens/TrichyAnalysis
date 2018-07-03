
rm(list=ls())

library(tidyverse)
library(washb)
library(caret)
library(foreign)
library(lubridate)
library(mgcv)

source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")



setwd("C:/Users/andre/Dropbox/Trichy analysis/Data")
load("Cleaned data/prior_weather_dataset.Rdata")

df<-read.dta("trichy_fu_vilwat.dta")

bl<-read.dta("trichy_bl.dta")
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
df <- left_join(df, bl, by="chhid")
table(df$vilid)

#fix errors in date input
df$vdatecol[31] <- "03/06/08"
df$vdatecol[87] <- "20/08/08"

df$intdate <- as.Date(df$vdatecol,format='%d/%m/%y')




#merge weather with survey
d <- merge(df, weather, by=c("intdate"), all.x=T, all.y=F) %>% arrange(vilid, intdate)


#Get study week
d$stdywk <- week(d$intdate)

(d$intdate[1:10])
week(d$intdate[1:10])

#
#Make unique studyweek so AR1 works
table(d$stdywk)
table(d$individ, d$stdywk)

# d <- d %>%  mutate(stdywk2 = stdywk) %>% arrange(individ, stdywk) %>% group_by(individ, stdywk) %>%
#   mutate(stdywk2 = stdywk2 + (row_number()-0.5)) %>% ungroup() %>% mutate(stdywk=stdywk2) %>%
#   subset(., select=-c(stdywk2)) %>% as.data.frame()
# table(d$individ, d$stdywk)

#merge in covariates
d <- left_join(d, covariates, by="vilid")
table(d$wpi)

Wvars<-c("wpi","villOD")

d$id <- d$individ

##########################
#Unadjusted analysis
##########################


# T1<-trichy_gamm(d, Y="H2S", A="tempQ7")
# T1

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
h2s.HR3_strat_adj





vilH2S_resdf <- rbind( 
 h2s.HR1, h2s.HR2, h2s.HR3,
 h2s.HR1_strat, h2s.HR2_strat, h2s.HR3_strat)


save(vilH2S_resdf, file="C:/Users/andre/Dropbox/Trichy analysis/Results/villageh2s.HRmean-gamm-results.Rdata")

