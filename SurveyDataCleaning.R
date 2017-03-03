rm(list=ls())

library(dplyr)
library(foreign)
library(dplyr)
library(h2o)
library(SuperLearner)
library(tmle)
library(tmle.npvi)
library(washb)


###Cleaning survey data
#read in data and subset to relevant variables (excluding covariates to prescreen)

setwd("C:/Users/andre/Documents/Trichy analysis")
svy<-read.dta("trichy_long_old.dta") %>%
  subset(select=c("vilid","hhid","individ","wpi","round","intdate","stdywk","bdate","mdiar7d", "diar2d", "diar7d", "diar14d","hcgi7d","diardays")) %>%
  filter(!is.na(diar7d))


#tabulate the outcomes
table(svy$diar2d)
table(svy$diar7d)
table(svy$diar14d)

table(svy$mdiar7d)
table(svy$hcgi7d)

table(svy$diardays)


preW<-read.dta("trichy_long_old.dta") %>%
  filter(!is.na(diar7d)) %>%
  subset(select=c("round","age","sex",
                  "bfcur","sanOD","sanlat","primwat",
                  "soil","thatch","radio","tv","cell","bank","moto","bike","mosq","stove","fuel1","fuel2", #assets
                  "anygrp","shgw","agri","elec","land","ownhome",
                  "pceduc","pclit","momage","momwork", # education
                  "rooms","totp","sc","kitchen","kitchvent", #house
                  "hwflies","hwwater","hwsoap","hwash","hwtow","hwsink","hwcntwat","hwcntsoap", #handwashing station
                  "hw1","hw2","hw3","hw4","hw5","hw6","hw7","hw8","hw9","hw10","hw11","hw12", #handwash w/ water
                  "hws1","hws2","hws3","hws4","hws5","hws6","hws7","hws8","hws9","hws10","hws11","hws12", #handwash w/soap
                  "boil_freq","boil_ready","boil_min","boil_present","wqsource","wqcol","wstore","wboil", #watertreat
                  #"ecoli","log10ecoli","tc","log10tc","dilution","h2s", #water contamination, 
                         #(don't include, because on causal pathway between weather and diar)
                  "latage","latfinance","latused","lathole","latwat","latsoap","latflies","latfec", #latrine quality
                  "an_buff","an_cow","an_ox","an_calf","an_goat","an_chick", #household animals
                  "an_dogcat","stockliv","goatliv","chickliv","dogcatliv", #household animals
                  "fecesliv","fecessm", #visible fecal contamination
                  "chobs_hand","chobs_nails","chobs_face","chobs_cloth","chobs_nocloth","chobs_shoes",
                  "qF1")) #asset composite index



#Check if any variable is  missing >90% 
too.missing<-which(apply(preW, 2, function(x) sum(is.na(x))/length(x))>0.9)
too.missing

preW<-subset(preW, select=-too.missing)

#Check if variables are missing any variation
apply(preW, 2, function(x) length(table(x)))

#Replace missingness with missing category 99 if the variable is binary or a factor
#set aside continious variables
preWcont<-subset(preW, select=c(age,momage))
W<-subset(preW, select= -c(age,momage))


table(is.na(W))

for(i in 1:ncol(W)){
  W[,i]<-factor(W[,i])
  if(length(is.na(W[,i]))>0){
    W[,i]<-addNA(W[,i])
    levels(W[,i])[length(levels(W[,i]))]<-"miss"
  }
}

table(is.na(W))

#Add back in the continious variables
preW <- cbind(W,preWcont)
head(preW)



#Two variables are causing error in prescreen function due to sparsity: "wqsource"  "an_dogcat"
#Collapse rare factor levels here

table(preW$wqsource)
levels(preW$wqsource)
to.replace<-preW$wqsource==levels(preW$wqsource)[1]|preW$wqsource==levels(preW$wqsource)[6]|preW$wqsource==levels(preW$wqsource)[7]|is.na(preW$wqsource) #-99 and 95 are missing vars, and 5 is drink from river- set as missing here
preW[to.replace,"wqsource"]<-levels(preW$wqsource)[8]
table(preW$wqsource)

table(preW$an_dogcat) #change all number of animal variables to indicators for presence

animal.vars<-
transmute(preW,
  buff=as.numeric(an_buff!=0),  
  cow=as.numeric(an_cow!=0), 
  ox=as.numeric(an_ox!=0), 
  calf=as.numeric(an_calf!=0), 
  goat=as.numeric(an_goat!=0), 
  chick=as.numeric(an_chick!=0), 
  dogcat=as.numeric(an_dogcat!=0))

preW<-cbind(preW,animal.vars) %>%
  select(-(an_buff:an_dogcat))

head(preW)

#prescreen via likelihood ratio function
Wscreen <- washb_prescreen(Y=svy$diar7d,Ws=preW,family="binomial", pval = 0.2 ,print=T)


#subset W dataframe to prescreen selected variables
W <- subset(preW,select=Wscreen)

#Add covariates into dataframe with the outcome
df<-cbind(svy,W)
glimpse(df)

names.svy<-colnames(svy)
names.W<-colnames(W)
save(df,names.svy,names.W,file="prescreened_blcov.Rdata")