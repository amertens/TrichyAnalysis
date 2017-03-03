
rm(list=ls())

library(dplyr)
library(SuperLearner)
library(tmle)
library(tmle.npvi)
library(washb)
library(h2o)
library(h2oEnsemble)
library(cvAUC) 


setwd("C:/Users/andre/Documents/Trichy analysis")
load("prescreened_blcov.Rdata")

#subset to complete cases
dim(df)
df<-subset(df, select=c("diar7d","age","primwat","thatch","bfcur","anygrp"))

df<-df[complete.cases(df),]
dim(df)

#data<-cbind(df$diar7d,W)
df$diar7d<-as.numeric(df$diar7d)
df$thatch<-as.numeric(df$thatch)

df<-df[1:1000,] #subset for speed

W=subset(df, select=-c(thatch,diar7d))
W<-design_matrix(W)
A=df$thatch-1
Y=df$diar7d

  test<-tmle(Y=Y, A=A,W=W,Q.SL.library = c("SL.glm") ,g.SL.library = c("SL.glm"), family="binomial")
  test
  
  Q.data<-cbind(Y,A,W)
  Q.glm<-glm(Y~. ,data=Q.data, family="binomial")
  summary(Q.glm)
  g.data<-cbind(A,W)
  g.glm<-glm(A~. ,data=g.data, family="binomial")
  summary(g.glm)
  
  newdata0<-newdata1<-Q.data
  newdata1$A<-1
  newdata0$A<-0

  Q1<-predict(Q.glm,newdata=newdata1,type="response")
  Q0<-predict(Q.glm,newdata=newdata0,type="response")
  Q<-cbind(Q0,Q1)
  g<-predict(g.glm,type="response")
  
  test2<-tmle(Y=Y, A=A,W=W,Q=Q,g1W=g, family="binomial")
  test2

  #Consistent results between tmle package and glm predictions with 
