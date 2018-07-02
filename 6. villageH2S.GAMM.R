
rm(list=ls())

library(tidyverse)
library(washb)
library(caret)
library(foreign)
library(lubridate)

source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")



setwd("C:/Users/andre/Dropbox/Trichy analysis/Data")
df<-read.dta("trichy_fu_vilwat.dta")


#load("Cleaned data/analysis_datasets.Rdata")
load("Cleaned data/prior_weather_dataset.Rdata")


###################
#village Survey data processing
###################


#clean and select needed variables
df <- df %>% subset(., select=c(sid, svy, vdatecol, vdatearr, v10)) %>% rename(Y=v10, id=sid) %>% filter(!is.na(Y))
head(df)

#fix errors in date input
df$vdatecol[31] <- "03/06/08"
df$vdatecol[87] <- "20/08/08"

df$intdate <- as.Date(df$vdatecol,format='%d/%m/%y')
# df$year <- as.numeric(year(df$date))
# df$month <- as.numeric(month(df$date))
# df$day <- as.numeric(days(df$date))


#merge weather with survey
d <- merge(df, weather, by=c("intdate"), all.x=T, all.y=F) %>% arrange(id, intdate)


#Get study week
d$stdywk <- week(d$intdate)
d$individ <- d$vilid <- d$id

(d$intdate[1:10])
week(d$intdate[1:10])

#
#Make unique studyweek so AR1 works
table(d$stdywk)
table(d$individ, d$stdywk)

d <- d %>%  mutate(stdywk2 = stdywk) %>% arrange(individ, stdywk) %>% group_by(individ, stdywk) %>%
  mutate(stdywk2 = stdywk2 + (row_number()-0.5)) %>% ungroup() %>% mutate(stdywk=stdywk2) %>%
  subset(., select=-c(stdywk2)) %>% as.data.frame()
table(d$stdywk)

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




vilH2S_resdf <- rbind( 
 h2s.HR1, h2s.HR2, h2s.HR3,
 h2s.HR1_strat, h2s.HR2_strat, h2s.HR3_strat)


save(vilH2S_resdf, file="C:/Users/andre/Dropbox/Trichy analysis/Results/villageh2s.HRmean-gamm-results.Rdata")

