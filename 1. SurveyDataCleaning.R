



rm(list=ls())

library(tidyverse)
library(foreign)
library(washb)
library(caret)


###Cleaning survey data
#read in data and subset to relevant variables (excluding covariates to prescreen)
setwd("C:/Users/andre/Dropbox/Trichy analysis/Data")
svy<-read.dta("trichy_long.dta") %>%
  subset(select=c("vilid","hhid","individ","wpi","round","intdate","stdywk","bdate","mdiar7d", "diar2d", "diar7d", "diar14d","hcgi7d","diardays","h2s")) %>%
  filter(!is.na(diar7d))



#tabulate the outcomes
table(svy$diar2d)
table(svy$diar7d)
table(svy$diar14d)

table(svy$mdiar7d)
table(svy$hcgi7d)

table(svy$diardays)


preW<-read.dta("trichy_long.dta") %>%
  filter(!is.na(diar7d)) %>%
  subset(select=c("vilid",
                  "round","age","sex", 
                  "bfcur","sanOD","sanlat","primwat",
                  "soil","thatch","radio","tv","cell","bank","moto","bike","mosq","stove","fuel1","fuel2", #assets
                  "anygrp","shgw","agri","elec","land","ownhome",
                  "pceduc","pclit","momage","momwork", # education
                  "rooms","totp","sc","kitchen","kitchvent", #house
                  "hwflies","hwwater","hwsoap","hwash","hwtow","hwsink", #handwashing station
                  #"hwcntwat","hwcntsoap", Don't include count of HW behavior. Great missingness, and redundant with HW indicators
                  "hw1","hw2","hw3","hw4","hw5","hw6","hw7","hw8","hw9","hw10","hw11","hw12", #handwash w/ water
                  "hws1","hws2","hws3","hws4","hws5","hws6","hws7","hws8","hws9","hws10","hws11","hws12", #handwash w/soap
                  "boil_freq","boil_ready","boil_min","boil_present","wqsource","wqcol","wstore","wboil", #watertreat
                  #"ecoli","log10ecoli","tc","log10tc","dilution","h2s", #water contamination, 
                         #(don't include, because on causal pathway between weather and diar)
                  "latage","latfinance","latused","lathole","latwat","latsoap","latflies","latfec", #latrine quality
                  "an_buff","an_cow","an_ox","an_calf","an_goat","an_chick", #household animals
                  "an_dogcat","stockliv","goatliv","chickliv","dogcatliv", #household animals
                  #"fecesliv","fecessm", #visible fecal contamination; don't include because may vary with rainfall/muddiness
                  #"chobs_hand","chobs_nails","chobs_face","chobs_cloth","chobs_nocloth","chobs_shoes", #Don't include because may vary with rainfall/muddiness
                  "qF1")) #asset composite index

#Calculate village level open defication
villOD <- preW %>% group_by(vilid) %>% summarize(villOD=mean(sanOD)) %>% as.data.frame()
summary(villOD$villOD)
preW <- left_join(preW, villOD, by="vilid")
preW <- preW %>% subset(., select = -c(vilid))

#table 1 summary
tab1 <- read.dta("trichy_long.dta") %>% filter(round==0) %>% 
    filter(!is.na(diar7d))%>% 
  mutate(animals_owned=an_buff+an_cow+an_ox+an_calf+an_goat+an_chick+an_dogcat) %>%
  subset(select=c("age","ageyrs","wpi","sex", "bfcur",
                  "momage","pclit", "pceduc","momwork",
                  "primwat",
                  "totp", "elec",
                  "soil","land",
                  "rooms","animals_owned",
                  "sc","kitchen","kitchvent", #house
                  "hwsoap","latsoap",
                  "sanOD","sanlat","fecessm","fecesliv",
                   "agri"
                   )) 
tab1 %>% mutate_all(as.numeric) %>% summarize_all(funs(mean), na.rm=T) 
prop.table(table(tab1$pceduc))
prop.table(table(tab1$primwat))

fivenum(tab1$age)
fivenum(tab1$ageyrs)

#Check and drop if any variable is  missing >80% 
too.missing<-which(apply(preW, 2, function(x) sum(is.na(x))/length(x))>0.5)
too.missing

preW<-subset(preW, select=-too.missing)



table(preW$an_dogcat) #change all number of animal variables to indicators for presence

#Convert animal variables from (sparse) counts to indicators
animal.vars<-
transmute(preW,
  buff=as.numeric(an_buff!=0),  
  cow=as.numeric(an_cow!=0), 
  ox=as.numeric(an_ox!=0), 
  calf=as.numeric(an_calf!=0), 
  goat=as.numeric(an_goat!=0), 
  chick=as.numeric(an_chick!=0), 
  dogcat=as.numeric(an_dogcat!=0))

animal.vars$buff[is.na(preW$an_buff)]<-NA
animal.vars$cow[is.na(preW$an_cow)]<-NA
animal.vars$ox[is.na(preW$an_ox)]<-NA
animal.vars$calf[is.na(preW$an_calf)]<-NA
animal.vars$goat[is.na(preW$an_goat)]<-NA
animal.vars$chick[is.na(preW$an_chick)]<-NA
animal.vars$dogcat[is.na(preW$an_dogcat)]<-NA

preW<-cbind(preW,animal.vars) %>%
  subset(., select= -c(an_buff:an_dogcat))



#Check missingness of continious variables
table(is.na(preW$age))
table(is.na(preW$momage))
table(is.na(preW$villOD))


#Replace missingness with missing category 99 if the variable is binary or a factor
#set aside continious variables
preWcont<-subset(preW, select=c(age,momage, villOD))
W<-subset(preW, select= -c(age,momage, villOD))


table(is.na(W))

#Add missing category to categorical variables 
for(i in 1:ncol(W)){
  W[,i]<-factor(W[,i])
  if(length(is.na(W[,i]))>0){
    W[,i]<-addNA(W[,i])
    levels(W[,i])[length(levels(W[,i]))]<-"miss"
  }
}

table(is.na(W))

#Median impute continious variables (after adding missinginess indicator)

#Add in missingess indicator variables for covariates with missingness
  #Create dataframe for missingness indicators
  miss.ind<-data.frame(matrix(NA, nrow = dim(preWcont)[1], ncol = dim(preWcont)[2]))
  for(i in 1:ncol(preWcont)){
    miss.ind[,i] <- as.numeric(is.na(preWcont[,i]))
    preWcont[is.na(preWcont[,i]),i] <- median(preWcont[,i], na.rm=T)
  }
  colnames(miss.ind)<-paste0(colnames(preWcont),".miss")
  
  
#Add back in the continious variables
preW <- cbind(W,preWcont,miss.ind)
table(is.na(preW))
head(preW)


#Two variables are causing error in prescreen function due to sparsity: "wqsource"  "an_dogcat"
#Collapse rare factor levels here

table(preW$wqsource)
levels(preW$wqsource)
to.replace<-preW$wqsource==levels(preW$wqsource)[1]|preW$wqsource==levels(preW$wqsource)[6]|preW$wqsource==levels(preW$wqsource)[7]|is.na(preW$wqsource) #-99 and 95 are missing vars, and 5 is drink from river- set as missing here
preW[to.replace,"wqsource"]<-levels(preW$wqsource)[8]
table(preW$wqsource)

#drop empty levels
preW<-droplevels(preW)



#Remove W's with zero variance

#Check if variables are missing any variation
apply(preW, 2, function(x) length(table(x)))

dim(preW)
  preproc = caret::preProcess(preW, method = c("zv"))
  preW = predict(preproc, preW)
dim(preW)



#Add covariates into dataframe with the outcome
df<-cbind(svy,preW)
glimpse(df)

names.svy<-colnames(svy)
names.W<-colnames(W)



###################
#Set outcomes and 
###################

# Expand out factors into indicators 
dim(df)
df<-df[,-5] #remove duplicate "round" variable
head(df)

Y<-df$diar7d
table(is.na(df$individ))
table(is.na(df$hhid))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
id<-df$vilid  #SE's clusted on village id
#id<-df$hhid  #SE's clusted on house id
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

intdate<-df$intdate
Wfac<-subset(df,select= -c(hhid,intdate,bdate,mdiar7d, diar2d, diar7d, diar14d, hcgi7d, diardays, h2s ))
colnames(Wfac)

H2S<-df$h2s

#Check missingness of W
table(is.na(Wfac))

#Save variable names
Wvars <- colnames(Wfac)
Wvars <- Wvars[-c(1:2,4:5)] #drop ID variables from adjustment covariates
survey<-cbind(intdate,Y,id,H2S,Wfac)



save(df,names.svy,names.W,
     file="C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data/bl_covariates.Rdata")

save(survey, Wvars, 
     file="C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data/survey_dataset.Rdata")
