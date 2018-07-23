

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


#Check mean rain and heavy rain events bly lag
mean(d$rain.ave7.lag7>0, na.rm=T)
mean(d$rain.ave7.lag14>0, na.rm=T)
mean(d$rain.ave7.lag21>0, na.rm=T)

mean(d$HeavyRain.lag7)
mean(d$HeavyRain.lag14)
mean(d$HeavyRain.lag21)

#Check correlation of 
 cor(d$temp.ave7.lag7, d$rain.ave7.lag7) 
 cor(d$temp.ave7.lag14, d$rain.ave7.lag14) 
 cor(d$temp.ave7.lag21, d$rain.ave7.lag21) 

 
 model <- lm(mpg ~ disp + hp + wt + qsec, data = mtcars)
ols_vif_tol(model)


#Check for duplicate study weeks
d %>% group_by(individ, round, stdywk) %>% summarise(n=n()) %>% ungroup() %>% summarise(mean(n))


#Make unique studyweek so AR1 works
# d <- d %>%  arrange(individ, round) %>% group_by(individ, round) %>%
#   mutate(stdywk = stdywk + (row_number()-1)/2)

#Temp
d <- d %>% arrange(individ, round, stdywk, rev(Y)) %>% distinct(., individ, stdywk, .keep_all=T)



#Temperature
T1<-trichy_gamm(d, A="tempQ7")
T1

T2<-trichy_gamm(d, A="tempQ14")
T2

T3<-trichy_gamm(d, A="tempQ21")
T3

#Rainfall
HR1<-trichy_gamm(d, A="HeavyRain.lag7")
HR1

HR2<-trichy_gamm(d, A="HeavyRain.lag14")
HR2

HR3<-trichy_gamm(d, A="HeavyRain.lag21")
HR3

#stratified

#NOTE: change function to include interaction term and extract results from that
HR1_strat<-trichy_gamm(d, A="HeavyRain.lag7", strat="LT8_T")
HR1_strat

HR2_strat<-trichy_gamm(d, A="HeavyRain.lag14", strat="LT15_T")
HR2_strat

HR3_strat<-trichy_gamm(d, A="HeavyRain.lag21", strat="LT22_T")
HR3_strat



#H2S



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


#Possion log link for PR
  
  m1 <- gamm(Y ~ tempQ7 + HeavyRain.lag7, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res1<-summary(m1$gam)
  res1

    m1 <- gamm(Y ~ tempQ7 + HeavyRain.lag7, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res1<-summary(m1$gam)
  res1


#model with both rain and temp

  m1 <- gamm(Y ~ tempQ7 + HeavyRain.lag7, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res1<-summary(m1$gam)
  res1
  
  m2 <- gamm(Y ~ tempQ14 + HeavyRain.lag14, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res2<-summary(m2$gam)
  res2
  
  m3 <- gamm(Y ~ tempQ21 + HeavyRain.lag21, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res3<-summary(m3$gam)
  res3
  
  
  # sm1 <- gamm(Y ~ tempQ7 + HeavyRain.lag7 + s(stdywk), data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  # sres1<-summary(sm1$gam)
  # sres1
  # 
  # sm2 <- gamm(Y ~ tempQ14 + HeavyRain.lag14 + s(stdywk), data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  # sres2<-summary(sm2$gam)
  # sres2
  # 
  # sm3 <- gamm(Y ~ tempQ21 + HeavyRain.lag21 + s(stdywk), data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  # sres3<-summary(sm3$gam)
  # sres3
  
  
#stratified by long-term rain trends
  m1 <- gamm(Y ~ tempQ7 + HeavyRain.lag7*LT8_T, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res1_strat<-summary(m1$gam)
  res1_strat
  
  m2 <- gamm(Y ~ tempQ14 + HeavyRain.lag14*LT15_T, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res2_strat<-summary(m2$gam)
  res2_strat
  
  m3 <- gamm(Y ~ tempQ21 + HeavyRain.lag21*LT22_T, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res3_strat<-summary(m3$gam)
  res3_strat
  
  

 #stratified by intervention
  m1 <- gamm(Y ~ tempQ7*wpi + HeavyRain.lag7*wpi, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res1_int<-summary(m1$gam)
  res1_int
  
  m2 <- gamm(Y ~ tempQ14*wpi + HeavyRain.lag14*wpi, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res2_int<-summary(m2$gam)
  res2_int
  
  m3 <- gamm(Y ~ tempQ21*wpi + HeavyRain.lag21*wpi, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res3_int<-summary(m3$gam)
  res3_int
  
  
  
#90th percentile of all days (lower threshold than 80th percentile of rainy days)
HR1<-trichy_gamm(d, A="HeavyRain90.lag7")
HR1

HR2<-trichy_gamm(d, A="HeavyRain90.lag14")
HR2

HR3<-trichy_gamm(d, A="HeavyRain90.lag21")
HR3

#stratified
HR1_strat<-trichy_gamm(d, A="HeavyRain90.lag7", strat="LT8_T")
HR1_strat

HR2_strat<-trichy_gamm(d, A="HeavyRain90.lag14", strat="LT15_T")
HR2_strat

HR3_strat<-trichy_gamm(d, A="HeavyRain90.lag21", strat="LT22_T")
HR3_strat

  

#subgroup analyses


Temp1_trt<-trichy_gamm(d, A="tempQ7", strat="wpi")
Temp1_trt

HR1_trt<-trichy_gamm(d, A="HeavyRain.lag7", strat="wpi")
HR1_trt


Temp1_sanlat<-trichy_gamm(d, A="tempQ7", strat="sanlat")
Temp1_sanlat

HR1_sanlat<-trichy_gamm(d, A="HeavyRain.lag7", strat="sanlat")
HR1_sanlat



#Check interaction effects
  require(lmtest)

m2 <- gamm(Y ~ tempQ7 * wpi , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
m1 <- gamm(Y ~ tempQ7 , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
anova(m2$gam, m1$gam)

m2 <- gamm(Y ~ HeavyRain.lag7 * wpi , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
m1 <- gamm(Y ~ HeavyRain.lag7 , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
anova(m2$gam, m1$gam)


m4 <- gamm(Y ~ tempQ7 * sanlat , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
m3 <- gamm(Y ~ tempQ7 , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
anova(m4$gam, m3$gam)




m2 <- gamm(Y ~ tempQ7 * wpi + HeavyRain.lag7*wpi , data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
m1 <- gamm(Y ~ tempQ7 + HeavyRain.lag7, data=d, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
anova(m2$gam, m1$gam)


# 
# #Adding in interaction extraction
# 
# 
# HR1_strat<-trichy_gamm(d, 
#                        A="HeavyRain.lag7"
#                        strat="LT8_T"
# 
#       
#   df <-d %>% subset(., select=c("Y","stdywk", "individ", "vilid",Avar, strat))
#   colnames(df)[colnames(df)==Avar] <- "A"
#   colnames(df)[colnames(df)==strat] <- "V"
#   df <-df[!is.na(df$A),]
#   nlevels <- length(unique(df$A))
#               
#                        
#                        
#   m <- gamm(Y ~ A*V , data=df, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
#  
#   res<-summary(m$gam)
#      res 
#     
#      
#         est <- exp(t(lc) %*% x[, 4])
#         se <- sqrt(t(lc) %*% vcv %*% lc)
#         Z <- log(est)/se
#         P <- 2 * pnorm(-abs(Z))
#         lb <- exp(log(est) - 1.96 * se)
#         ub <- exp(log(est) + 1.96 * se)
#         res <- matrix(c(est, se, lb, ub, Z, P), nrow = 1)
#         res[1:6] <- round(res[1:6], 4)
#         colnames(res) <- c("RR", "se.RR", "RR.lb", "RR.ub", "Z", 
#             "P")
# 
# #Make GAMM analysis function
# gammAR1 <- function(data, Avar, Yvar="Diarrhea", strat=NULL){
#   df <-data %>% subset(., select=c("Y","stdywk", "individ", "vilid",Avar))
#   colnames(df)[colnames(df)==Avar] <- "A"
#   df <-df[!is.na(df$A),]
#   nlevels <- length(unique(df$A))
#   
#     m <- gamm(Y ~ A , data=df, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
#  
#   res<-summary(m$gam)
#   #print(res)
#   resdf <- data.frame(outcome=rep(Yvar,nlevels-1), exposure=rep(Avar,nlevels-1), 
#                       PR=exp(res$p.coeff)[2:nlevels], 
#                       ci.lb=exp(res$p.coeff[2:nlevels] - res$se[2:nlevels] * 1.96), 
#                       ci.ub=exp(res$p.coeff[2:nlevels] + res$se[2:nlevels] * 1.96))
#   if(is.null(strat)){
#     resdf$strat="Unstratified"
#   }else{
#     resdf$strat=strat
#   }
#   return(resdf)
# }
# 
# trichy_gamm <- function(d, A, Y="Diarrhea", strat=NULL){
#   if(is.null(strat)){
#     res<-gammAR1(data=d, Avar=A, Yvar=Y, strat=NULL)
#   }else{
#     V= d[,strat]
#     d<-d[!is.na(V),]
#     V<-V[!is.na(V)]
#     res<-NULL
#     for(i in 1:length(unique(V))){
#       dsub= d[V==levels(V)[i],]
#       res_strat <- gammAR1(data=dsub, Avar=A, Yvar=Y, strat=levels(V)[i], nospline=nospline)
#       res<-rbind(res, res_strat)
#     }
#     
#   }
#   return(res)
# }
# 
