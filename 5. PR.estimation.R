

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
d <- d %>%  mutate(stdywk2 = stdywk) %>% arrange(individ, stdywk, age) %>% group_by(individ, stdywk) %>%
  mutate(stdywk2 = stdywk2 + (row_number()-1)) %>% ungroup() %>% mutate(stdywk=stdywk2) %>%
  subset(., select=-c(stdywk2)) %>% as.data.frame()




#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#Diarrhea
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


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
summary_fun(d, "tempQ7", "temp.ave7.lag7")
summary_fun(d, "tempQ14", "temp.ave7.lag14")
summary_fun(d, "tempQ21", "temp.ave7.lag21")

summary_fun(d, "HeavyRain.lag7","rain.ave7.lag8")
summary_fun(d, "HeavyRain.lag14","rain.ave7.lag15")
summary_fun(d, "HeavyRain.lag21","rain.ave7.lag22")

summary_fun(d, "HeavyRain.lag7","rain.ave7.lag8", "LT8_T", "LT8")
summary_fun(d, "HeavyRain.lag14","rain.ave7.lag15", "LT15_T", "LT15")
summary_fun(d, "HeavyRain.lag21","rain.ave7.lag22", "LT22_T", "LT22")


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
T1_adj<-trichy_gamm(d, A="tempQ7", weathervar="rain.ave7.lag8", Wvars = Wvars)
T1_adj

T2_adj<-trichy_gamm(d, A="tempQ14", weathervar="rain.ave7.lag15", Wvars = Wvars)
T2_adj

T3_adj<-trichy_gamm(d, A="tempQ21", weathervar="rain.ave7.lag22", Wvars = Wvars)
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

#Number of households samples by round
d %>% group_by(hhid, round) %>% slice(1) %>% mutate(hasSample=!is.na(H2S)) %>% group_by(round) %>% summarize(N=n(), mn=round(mean(hasSample)*100,1), perctot=round(sum(hasSample)/900*100,1))

#Grab one obs per household per round to avoid double counting the samples from households with multiple children.
dH2S <- d %>% group_by(hhid, round) %>% slice(1) %>% subset(., select=-Y) %>% filter(!is.na(H2S)) %>% rename(Y=H2S) %>% as.data.frame()

#summary statistics
summary_fun(dH2S, "tempQ1", "temp.ave7.lag1")
summary_fun(dH2S, "tempQ7", "temp.ave7.lag7")
summary_fun(dH2S, "tempQ14", "temp.ave7.lag14")

summary_fun(dH2S, "HeavyRain.lag1","rain.ave7.lag1")
summary_fun(dH2S, "HeavyRain.lag7","rain.ave7.lag8")
summary_fun(dH2S, "HeavyRain.lag14","rain.ave7.lag15")

summary_fun(dH2S, "HeavyRain.lag1","rain.ave7.lag1", "LT1_T", "LT1")
summary_fun(dH2S, "HeavyRain.lag7","rain.ave7.lag8", "LT8_T", "LT8")
summary_fun(dH2S, "HeavyRain.lag14","rain.ave7.lag15", "LT15_T", "LT15")

#Prevalence by sampling round
dH2S %>% group_by(round, vilid) %>% summarise(mn=mean(Y)) %>% ungroup %>% summarise(min(mn), mean(mn)) 
dH2S %>% group_by(round) %>% summarise(mn=mean(Y)) %>% ungroup %>% summarise(min(mn), mean(mn)) 

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

#Adjusted
h2s.T1_adj<-trichy_gamm(dH2S, Y="H2S", A="tempQ1", weathervar="rain.ave7.lag1",  Wvars = Wvars)
h2s.T1_adj

h2s.T2_adj<-trichy_gamm(dH2S, Y="H2S", A="tempQ7", weathervar="rain.ave7.lag8", Wvars = Wvars)
h2s.T2_adj

h2s.T3_adj<-trichy_gamm(dH2S, Y="H2S", A="tempQ14", weathervar="rain.ave7.lag15", Wvars = Wvars)
h2s.T3_adj


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



# save H2S results
save(h2s.T1, h2s.T2, h2s.T3,
     h2s.T1_adj, h2s.T2_adj, h2s.T3_adj,
     h2s.HR1, h2s.HR2, h2s.HR3,
     h2s.HR1_strat, h2s.HR2_strat, h2s.HR3_strat,
     h2s.HR1_adj, h2s.HR2_adj, h2s.HR3_adj,
     h2s.HR1_strat_adj, h2s.HR2_strat_adj, h2s.HR3_strat_adj, 
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_H2S_results.Rdata")



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

res_H2S_adj <- rbind(h2s.T1_adj, h2s.T2_adj, h2s.T3_adj,
                     h2s.HR1_adj, h2s.HR2_adj, h2s.HR3_adj,
 h2s.HR1_strat_adj, h2s.HR2_strat_adj, h2s.HR3_strat_adj)





save(res_prim, res_prim_adj, res_H2S, res_H2S_adj,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_plot_dfs.Rdata")

res_prim$PR <- round(res_prim$PR, 2)
res_prim$ci.lb <- round(res_prim$ci.lb, 2)
res_prim$ci.ub <- round(res_prim$ci.ub, 2)
res_prim$y <- paste0(res_prim$PR," (", res_prim$ci.lb,", ", res_prim$ci.ub,")")

df <- res_prim[,c(2,8)]
df
noquote(df[,2])

res_prim_adj$PR <- round(res_prim_adj$PR, 2)
res_prim_adj$ci.lb <- round(res_prim_adj$ci.lb, 2)
res_prim_adj$ci.ub <- round(res_prim_adj$ci.ub, 2)
res_prim_adj$y <- paste0(res_prim_adj$PR," (", res_prim_adj$ci.lb,", ", res_prim_adj$ci.ub,")")

df <- res_prim_adj[,c(2,8)]
df
noquote(df[,2])

res_H2S$PR <- round(res_H2S$PR, 2)
res_H2S$ci.lb <- round(res_H2S$ci.lb, 2)
res_H2S$ci.ub <- round(res_H2S$ci.ub, 2)
res_H2S$y <- paste0(res_H2S$PR," (", res_H2S$ci.lb,", ", res_H2S$ci.ub,")")

df <- res_H2S[,c(2,8)]
df
noquote(df[,2])

res_H2S_adj$PR <- round(res_H2S_adj$PR, 2)
res_H2S_adj$ci.lb <- round(res_H2S_adj$ci.lb, 2)
res_H2S_adj$ci.ub <- round(res_H2S_adj$ci.ub, 2)
res_H2S_adj$y <- paste0(res_H2S_adj$PR," (", res_H2S_adj$ci.lb,", ", res_H2S_adj$ci.ub,")")

df <- res_H2S_adj[,c(2,8)]
df
noquote(df[,2])

