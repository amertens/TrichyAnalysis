

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
sink("./PR.CC.results.txt")

load("analysis_datasets_CC.Rdata")
d<-dcc

load("C:/Users/andre/Dropbox/Trichy analysis/Results/adjustment_sets.Rdata")


#Make unique studyweek so AR1 works
d <- d %>%  mutate(stdywk2 = stdywk) %>% arrange(individ, stdywk, age) %>% group_by(individ, stdywk) %>%
  mutate(stdywk2 = stdywk2 + (row_number()-1)) %>% ungroup() %>% mutate(stdywk=stdywk2) %>%
  subset(., select=-c(stdywk2)) %>% as.data.frame()




#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#Diarrhea
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


#---------------------------
#Temperature
#---------------------------


#Adjusted
T1_adj<-trichy_gammCC(d, A="tempQ8", adj_set=res_prim_adj$T1_adj , weathervar="rain.ave7.lag8", Wvars = Wvars)
T1_adj

T2_adj<-trichy_gammCC(d, A="tempQ15", adj_set=res_prim_adj$T2_adj , weathervar="rain.ave7.lag15", Wvars = Wvars)
T2_adj

T3_adj<-trichy_gammCC(d, A="tempQ22", adj_set=res_prim_adj$T3_adj , weathervar="rain.ave7.lag22", Wvars = Wvars)
T3_adj

#---------------------------
#Rainfall
#---------------------------



#Adjusted
HR1_adj<-trichy_gammCC(d, A="HeavyRain.lag8", adj_set=res_prim_adj$HR1_adj , weathervar="temp.ave7.lag8", Wvars = Wvars)
HR1_adj

HR2_adj<-trichy_gammCC(d, A="HeavyRain.lag15", adj_set=res_prim_adj$HR2_adj , weathervar="temp.ave7.lag15", Wvars = Wvars)
HR2_adj

HR3_adj<-trichy_gammCC(d, A="HeavyRain.lag22", adj_set=res_prim_adj$HR3_adj , weathervar="temp.ave7.lag22", Wvars = Wvars)
HR3_adj

#stratified
adj_list=list(res_prim_adj$HR1_strat_adjT1, res_prim_adj$HR1_strat_adjT2, res_prim_adj$HR1_strat_adjT3)
HR1_strat_adj<-trichy_gammCC(d, A="HeavyRain.lag8", strat="LT8_T", adj_set=adj_list , weathervar="temp.ave7.lag8", Wvars = Wvars)
HR1_strat_adj

adj_list=list(res_prim_adj$HR2_strat_adjT1, res_prim_adj$HR2_strat_adjT2, res_prim_adj$HR2_strat_adjT3)
HR2_strat_adj<-trichy_gammCC(d, A="HeavyRain.lag15", strat="LT15_T", adj_set=adj_list , weathervar="temp.ave7.lag15", Wvars = Wvars)
HR2_strat_adj

adj_list=list(res_prim_adj$HR3_strat_adjT1, res_prim_adj$HR3_strat_adjT2, res_prim_adj$HR3_strat_adjT3)
HR3_strat_adj<-trichy_gammCC(d, A="HeavyRain.lag22", strat="LT22_T", adj_set=adj_list , weathervar="temp.ave7.lag22", Wvars = Wvars)
HR3_strat_adj



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H2S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#Number of households samples by round
d %>% group_by(hhid, round) %>% slice(1) %>% mutate(hasSample=!is.na(H2S)) %>% group_by(round) %>% summarize(N=n(), mn=round(mean(hasSample)*100,1), perctot=round(sum(hasSample)/900*100,1))

#Grab one obs per household per round to avoid double counting the samples from households with multiple children.
dH2S <- d %>% group_by(hhid, round) %>% slice(1) %>% subset(., select=-Y) %>% filter(!is.na(H2S)) %>% rename(Y=H2S) %>% as.data.frame()

#---------------------------
#Temperature
#---------------------------


#Adjusted
h2s.T1_adj<-trichy_gammCC(dH2S, Y="H2S", A="tempQ1", adj_set=res_prim_adj$h2s.T1_adj , weathervar="rain.ave7.lag1",  Wvars = Wvars)
h2s.T1_adj

h2s.T2_adj<-trichy_gammCC(dH2S, Y="H2S", A="tempQ8", adj_set=res_prim_adj$h2s.T2_adj , weathervar="rain.ave7.lag8", Wvars = Wvars)
h2s.T2_adj

h2s.T3_adj<-trichy_gammCC(dH2S, Y="H2S", A="tempQ15", adj_set=res_prim_adj$h2s.T3_adj , weathervar="rain.ave7.lag15", Wvars = Wvars)
h2s.T3_adj


#---------------------------
#Rainfall
#---------------------------


#Adjusted
h2s.HR1_adj<-trichy_gammCC(dH2S, Y="H2S", A="HeavyRain.lag1", adj_set=res_prim_adj$h2s.HR1_adj , weathervar="temp.ave7.lag1", Wvars = Wvars)
h2s.HR1_adj

h2s.HR2_adj<-trichy_gammCC(dH2S, Y="H2S", A="HeavyRain.lag8", adj_set=res_prim_adj$h2s.HR2_adj , weathervar="temp.ave7.lag8", Wvars = Wvars)
h2s.HR2_adj

h2s.HR3_adj<-trichy_gammCC(dH2S, Y="H2S", A="HeavyRain.lag15", adj_set=res_prim_adj$h2s.HR3_adj , weathervar="temp.ave7.lag15", Wvars = Wvars)
h2s.HR3_adj

#stratified
h2s.HR1_strat_adj <- h2s.HR2_strat_adj <- h2s.HR3_strat_adj <- NULL

adj_list=list(res_prim_adj$h2s.HR1_strat_adjT1, res_prim_adj$h2s.HR1_strat_adjT2, res_prim_adj$h2s.HR1_strat_adjT3)
h2s.HR1_strat_adj<-trichy_gammCC(dH2S, Y="H2S", A="HeavyRain.lag1", strat="LT1_T", adj_set=adj_list , weathervar="temp.ave7.lag1", Wvars = Wvars)
h2s.HR1_strat_adj

adj_list=list(res_prim_adj$h2s.HR1_strat_adjT1, res_prim_adj$h2s.HR1_strat_adjT2, res_prim_adj$h2s.HR1_strat_adjT3)
h2s.HR2_strat_adj<-trichy_gammCC(dH2S, Y="H2S", A="HeavyRain.lag8", strat="LT8_T", adj_set=adj_list , weathervar="temp.ave7.lag8", Wvars = Wvars)
h2s.HR2_strat_adj

adj_list=list(res_prim_adj$h2s.HR1_strat_adjT1, res_prim_adj$h2s.HR1_strat_adjT2, res_prim_adj$h2s.HR1_strat_adjT3)
h2s.HR3_strat_adj<-trichy_gammCC(dH2S, Y="H2S", A="HeavyRain.lag15", strat="LT15_T", adj_set=adj_list , weathervar="temp.ave7.lag15", Wvars = Wvars)
h2s.HR3_strat_adj




#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# save results
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ls()
res <- ls()
res <- paste(res[!(res %in% c("d","dH2S", "gammAR1","gammAR1_adj"))],sep="")
save(list=res, file="C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_resultsCC.Rdata")


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Compile results for printing
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

ls()

#primary 
res_prim_adj <- rbind(T1_adj$`resdf`, T2_adj$`resdf`, T3_adj$`resdf`, 
                      HR1_adj$`resdf`, HR2_adj$`resdf`, HR3_adj$`resdf`,
                      HR1_strat_adj[[1]],HR1_strat_adj[[6]],HR1_strat_adj[[11]],
                      HR2_strat_adj[[1]],HR2_strat_adj[[6]],HR2_strat_adj[[11]],
                      HR3_strat_adj[[1]],HR3_strat_adj[[6]],HR3_strat_adj[[11]])

#h2s
res_H2S_adj <- rbind(h2s.T1_adj$`resdf`, h2s.T2_adj$`resdf`, h2s.T3_adj$`resdf`,
                     h2s.HR1_adj$`resdf`, h2s.HR2_adj$`resdf`, h2s.HR3_adj$`resdf`,
                     h2s.HR1_strat_adj[[1]],h2s.HR1_strat_adj[[6]],h2s.HR1_strat_adj[[11]],
                     h2s.HR2_strat_adj[[1]],h2s.HR2_strat_adj[[6]],h2s.HR2_strat_adj[[11]],
                     h2s.HR3_strat_adj[[1]],h2s.HR3_strat_adj[[6]],h2s.HR3_strat_adj[[11]])





#Print results
res_prim_adj$PR <- round(res_prim_adj$PR, 2)
res_prim_adj$ci.lb <- round(res_prim_adj$ci.lb, 2)
res_prim_adj$ci.ub <- round(res_prim_adj$ci.ub, 2)
res_prim_adj$y <- paste0(res_prim_adj$PR," (", res_prim_adj$ci.lb,", ", res_prim_adj$ci.ub,")")

df <- res_prim_adj[,c(2,7)]
df
noquote(df[,2])



res_H2S_adj$PR <- round(res_H2S_adj$PR, 2)
res_H2S_adj$ci.lb <- round(res_H2S_adj$ci.lb, 2)
res_H2S_adj$ci.ub <- round(res_H2S_adj$ci.ub, 2)
res_H2S_adj$y <- paste0(res_H2S_adj$PR," (", res_H2S_adj$ci.lb,", ", res_H2S_adj$ci.ub,")")

df <- res_H2S_adj[,c(2,7)]
df
noquote(df[,2])


sink()
