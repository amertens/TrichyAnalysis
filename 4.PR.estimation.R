

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
library(mgcv)


source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")

try(setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data"))

load("analysis_datasets.Rdata")

#Check missing weather data
table(is.na(d$tempQ21))
table(is.na(d$HeavyRain.lag21))
table(is.na(d$LT8_T))
table(is.na(d$LT15_T))
table(is.na(d$LT22_T))


#Check for duplicate study weeks
d %>% group_by(individ, round, stdywk) %>% summarise(n=n()) %>% ungroup() %>% summarise(mean(n))


#Make unique studyweek so AR1 works
d <- d %>%  mutate(stdywk2 = stdywk) %>% arrange(individ, stdywk) %>% group_by(individ, stdywk) %>%
  mutate(stdywk2 = stdywk2 + (row_number()-1)) %>% ungroup() %>% mutate(stdywk=stdywk2) %>%
  subset(., select=-c(stdywk2)) %>% as.data.frame()




#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#Diarrhea
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



#---------------------------
#Temperature
#---------------------------

#Unadjusted
T1<-trichy_gamm(d, A="tempQ7")
T1

T2<-trichy_gamm(d, A="tempQ14")
T2

T3<-trichy_gamm(d, A="tempQ21")
T3

#Adjusted
T1_adj<-trichy_gamm(d, A="tempQ7", weathervar="rain.ave7.lag7", Wvars = Wvars)
T1_adj

T2_adj<-trichy_gamm(d, A="tempQ14", weathervar="rain.ave7.lag14", Wvars = Wvars)
T2_adj

T3_adj<-trichy_gamm(d, A="tempQ21", weathervar="rain.ave7.lag21", Wvars = Wvars)
T3_adj

#---------------------------
#Rainfall
#---------------------------

#Unadjusted
HR1<-trichy_gamm(d, A="HeavyRain.lag7")
HR1

HR2<-trichy_gamm(d, A="HeavyRain.lag14")
HR2

HR3<-trichy_gamm(d, A="HeavyRain.lag21")
HR3

#stratified

#NOTE: maybe change function to include interaction term and extract results from that
HR1_strat<-trichy_gamm(d, A="HeavyRain.lag7", strat="LT8_T")
HR1_strat

HR2_strat<-trichy_gamm(d, A="HeavyRain.lag14", strat="LT15_T")
HR2_strat

HR3_strat<-trichy_gamm(d, A="HeavyRain.lag21", strat="LT22_T")
HR3_strat



#Adjusted
HR1_adj<-trichy_gamm(d, A="HeavyRain.lag7", weathervar="temp.ave7.lag7", Wvars = Wvars)
HR1_adj

HR2_adj<-trichy_gamm(d, A="HeavyRain.lag14", weathervar="temp.ave7.lag14", Wvars = Wvars)
HR2_adj

HR3_adj<-trichy_gamm(d, A="HeavyRain.lag21", weathervar="temp.ave7.lag21", Wvars = Wvars)
HR3_adj

#stratified
HR1_strat_adj<-trichy_gamm(d, A="HeavyRain.lag7", strat="LT8_T", weathervar="temp.ave7.lag7", Wvars = Wvars)
HR1_strat_adj

HR2_strat_adj<-trichy_gamm(d, A="HeavyRain.lag14", strat="LT15_T", weathervar="temp.ave7.lag14", Wvars = Wvars)
HR2_strat_adj

HR3_strat_adj<-trichy_gamm(d, A="HeavyRain.lag21", strat="LT22_T", weathervar="temp.ave7.lag21", Wvars = Wvars)
HR3_strat_adj




#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H2S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

dH2S <- d %>% subset(., select=-Y) %>% filter(!is.na(H2S)) %>% rename(Y=H2S) %>% as.data.frame()


#---------------------------
#Temperature
#---------------------------

#Unadjusted
h2s.T1<-trichy_gamm(dH2S, Y="H2S", A="tempQ1")
h2s.T1

h2s.T2<-trichy_gamm(dH2S, Y="H2S", A="tempQ7")
h2s.T2

h2s.T3<-trichy_gamm(dH2S, Y="H2S", A="tempQ14")
h2s.T3


#---------------------------
#Rainfall
#---------------------------

#Unadjusted
h2s.HR1<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag1")
h2s.HR1

h2s.HR2<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag7")
h2s.HR2

h2s.HR3<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag14")
h2s.HR3

#stratified

#NOTE: maybe change function to include interaction term and extract results from that
h2s.HR1_strat<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag1", strat="LT1_T")
h2s.HR1_strat

h2s.HR2_strat<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag7", strat="LT8_T")
h2s.HR2_strat

h2s.HR3_strat<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag14", strat="LT15_T")
h2s.HR3_strat



#Adjusted
h2s.HR1_adj<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag1", weathervar="temp.ave7.lag1", Wvars = Wvars)
h2s.HR1_adj

h2s.HR2_adj<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag7", weathervar="temp.ave7.lag7", Wvars = Wvars)
h2s.HR2_adj

h2s.HR3_adj<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag14", weathervar="temp.ave7.lag14", Wvars = Wvars)
h2s.HR3_adj

#stratified
h2s.HR1_strat_adj <- h2s.HR2_strat_adj <- h2s.HR3_strat_adj <- NULL

h2s.HR1_strat_adj<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag1", strat="LT1_T", weathervar="temp.ave7.lag1", Wvars = Wvars)
h2s.HR1_strat_adj

h2s.HR2_strat_adj<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag7", strat="LT8_T", weathervar="temp.ave7.lag7", Wvars = Wvars)
h2s.HR2_strat_adj

h2s.HR3_strat_adj<-trichy_gamm(dH2S, Y="H2S", A="HeavyRain.lag14", strat="LT15_T", weathervar="temp.ave7.lag14", Wvars = Wvars)
h2s.HR3_strat_adj


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Sensitivity Analyses
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Subgroup Analyses
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



# Min and max temp
minT1<-trichy_gamm(d, A="mintempQ7")
minT1

minT2<-trichy_gamm(d, A="mintempQ14")
minT2

minT3<-trichy_gamm(d, A="mintempQ21")
minT3


maxT1<-trichy_gamm(d, A="maxtempQ7")
maxT1

maxT2<-trichy_gamm(d, A="maxtempQ14")
maxT2

maxT3<-trichy_gamm(d, A="maxtempQ21")
maxT3


  
  
  

#stratified by intervention
#NOTE: maybe change function to include interaction term and extract results from that

T1_wpi<-trichy_gamm(d, A="tempQ7", strat="wpi")
T1_wpi

T2_wpi<-trichy_gamm(d, A="tempQ14", strat="wpi")
T2_wpi

T3_wpi<-trichy_gamm(d, A="tempQ21", strat="wpi")
T3_wpi

T1_wpi_adj<-trichy_gamm(d, A="tempQ7", strat="wpi", weathervar="rain.ave7.lag7", Wvars = Wvars)
T1_wpi_adj

T2_wpi_adj<-trichy_gamm(d, A="tempQ14", strat="wpi", weathervar="rain.ave7.lag14", Wvars = Wvars)
T2_wpi_adj

T3_wpi_adj<-trichy_gamm(d, A="tempQ21", strat="wpi", weathervar="rain.ave7.lag21", Wvars = Wvars)
T3_wpi_adj


HR1_wpi<-trichy_gamm(d, A="HeavyRain.lag7", strat="wpi")
HR1_wpi

HR2_wpi<-trichy_gamm(d, A="HeavyRain.lag14", strat="wpi")
HR2_wpi

HR3_wpi<-trichy_gamm(d, A="HeavyRain.lag21", strat="wpi")
HR3_wpi

HR1_wpi_adj<-trichy_gamm(d, A="HeavyRain.lag7", strat="wpi", weathervar="temp.ave7.lag7", Wvars = Wvars)
HR1_wpi_adj

HR2_wpi_adj<-trichy_gamm(d, A="HeavyRain.lag14", strat="wpi", weathervar="temp.ave7.lag14", Wvars = Wvars)
HR2_wpi_adj

HR3_wpi_adj<-trichy_gamm(d, A="HeavyRain.lag21", strat="wpi", weathervar="temp.ave7.lag21", Wvars = Wvars)
HR3_wpi_adj


#Stratified
d$LT8_wpi <- factor(paste0(d$LT8_T, "_", d$wpi))
d$LT15_wpi <- paste0(d$LT15_T, "_", d$wpi)
d$LT22_wpi <- paste0(d$LT22_T, "_", d$wpi)
d$LT8_wpi[d$LT8_wpi=="NA_0" | d$LT8_wpi=="NA_1"] <- NA
d$LT15_wpi[d$LT15_wpi=="NA_0" | d$LT15_wpi=="NA_1"] <- NA
d$LT22_wpi[d$LT22_wpi=="NA_0" | d$LT22_wpi=="NA_1"] <- NA


HR1_wpi_strat<-trichy_gamm(d, A="HeavyRain.lag7", strat="LT8_wpi")
HR1_wpi_strat

HR2_wpi_strat<-trichy_gamm(d, A="HeavyRain.lag14", strat="LT15_wpi")
HR2_wpi_strat

HR3_wpi_strat<-trichy_gamm(d, A="HeavyRain.lag21", strat="LT22_wpi")
HR3_wpi_strat

  
  
#90th percentile of all days (lower threshold than 80th percentile of rainy days)
HR1_90<-trichy_gamm(d, A="HeavyRain90.lag7")
HR1_90

HR2_90<-trichy_gamm(d, A="HeavyRain90.lag14")
HR2_90

HR3_90<-trichy_gamm(d, A="HeavyRain90.lag21")
HR3_90

#stratified
HR1_strat_90<-trichy_gamm(d, A="HeavyRain90.lag7", strat="LT8_T")
HR1_strat_90

HR2_strat_90<-trichy_gamm(d, A="HeavyRain90.lag14", strat="LT15_T")
HR2_strat_90

HR3_strat_90<-trichy_gamm(d, A="HeavyRain90.lag21", strat="LT22_T")
HR3_strat_90

  







#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# save results
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ls()
res <- ls()
res <- paste(res[!(res %in% c("d","dH2S", "gammAR1","gammAR1_adj"))],sep="")
save(list=res, file="C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_results.Rdata")


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Compile results for plotting
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

ls()

#primary unadjusted
res_prim <- rbind(T1, T2, T3, HR1, HR2, HR3, HR1_strat, HR2_strat, HR3_strat)

#primary adjusted
res_prim_adj <- rbind(T1_adj, T2_adj, T3_adj, 
                      HR1_adj, HR2_adj, HR3_adj,
                      HR1_strat_adj, HR2_strat_adj, HR3_strat_adj)

#h2s
res_H2S <- rbind( h2s.T1, h2s.T2, h2s.T3,
 h2s.HR1, h2s.HR2, h2s.HR3,
 h2s.HR1_strat, h2s.HR2_strat, h2s.HR3_strat)

res_H2S_adj <- rbind(h2s.HR1_adj, h2s.HR2_adj, h2s.HR3_adj,
h2s.HR1_strat_adj, h2s.HR2_strat_adj, h2s.HR3_strat_adj)

#subgroups



res_minmax_temp <- rbind(maxT1, maxT2, maxT3, minT1,  minT2, minT3)
res_wpi_strat <- rbind(HR1_wpi, HR2_wpi, HR3_wpi, HR1_wpi_adj, HR2_wpi_adj, HR3_wpi_adj)

#90% rain
res_HR90 <- rbind(HR1_90, HR2_90, HR3_90, HR1_strat_90, HR2_strat_90, HR3_strat_90)



save(res_prim, res_prim_adj, res_H2S, res_H2S_adj,
     res_minmax_temp, res_wpi_strat, res_HR90,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_plot_dfs.Rdata")

# res_H2S2 <- res_H2S
# res_H2S_adj2 <- res_H2S_adj
# load("C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_plot_dfs.Rdata")
# res_H2S <- res_H2S2
# res_H2S_adj <- res_H2S_adj2
