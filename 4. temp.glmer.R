

rm(list=ls())
if(!require("devtools")){
  install.packages("devtools", repos = "http://cran.us.r-project.org")
  devtools::install_github("hadley/devtools")
  library("devtools")}

if(!require("washb")){install_github("ben-arnold/washb")  ; library("washb")}
if(!require("dplyr")){install.packages("dplyr", repos = "http://cran.us.r-project.org"); library("dplyr")}
if(!require("lme4")){install.packages("lme4", repos = "http://cran.us.r-project.org"); library("lme4")}
if(!require("caret")){install.packages("caret", repos = "http://cran.us.r-project.org"); library("caret")}
if(!require("zoo")){install.packages("zoo", repos = "http://cran.us.r-project.org"); library("zoo")}


try(setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data"))

load("analysis_datasets.Rdata")


##########################
#Unadjusted analysis
##########################


#Mean diarrhea by quartile:
Mean_lag7<-d%>%
  group_by(quartile = ntile(temp.ave7.lag7, 4))%>%
  summarize(quart.n=n(),meantemp=mean(temp.ave7.lag7),meandiar=mean(Y), meanH2s=mean(H2S, na.rm=T), sumH2s=sum(H2S, na.rm=T))
Mean_lag14<-d%>%
  group_by(quartile = ntile(temp.ave7.lag14, 4))%>%
  summarize(quart.n=n(),meantemp=mean(temp.ave7.lag14),meandiar=mean(Y), meanH2s=mean(H2S, na.rm=T), sumH2s=sum(H2S, na.rm=T))
Mean_lag21<-d%>%
  group_by(quartile = ntile(temp.ave7.lag21, 4))%>%
  summarize(quart.n=n(),meantemp=mean(temp.ave7.lag21),meandiar=mean(Y), meanH2s=mean(H2S, na.rm=T), sumH2s=sum(H2S, na.rm=T))

#Save quartile cutoffs
temp_quartiles <- data.frame(tempQ=quantile(d$avetemp)[2:4],
                             tempQ_lag7=quantile(d$temp.ave7.lag7)[2:4],
                             tempQ_lag14=quantile(d$temp.ave7.lag14)[2:4],
                             tempQ_lag21=quantile(d$temp.ave7.lag21)[2:4])


# save(temp_quartiles, 
#      file="C:/Users/andre/Dropbox/Trichy analysis/Results/temp.quartiles.Rdata")


m1 = glmer(Y ~ tempQ7 + (1 | stdywk) + (1|vilid), data = d, family = "binomial")

res <- summary(m1)
res$coefficients

####################
#Quartile analysis function
####################


quart_glmner<-function(d,Wvars=NULL,weather="tempQ7",q=c(1,4),cutoff=4, library=NULL){
colnames(d)[which( colnames(d)==paste0(weather))]="weather"

# df<-d%>%
#   mutate(quartile = ntile(weather, cutoff))%>%
#   subset(quartile==q[1]|quartile==q[2])%>%
#   mutate(A=ifelse(quartile==q[2],1,0))

d$quartile<-cut(d$weather, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))

df <- d[d$quartile==levels(d$quartile)[q[1]]|
           d$quartile==levels(d$quartile)[q[2]],]
df$A=ifelse(df$quartile==levels(df$quartile)[q[2]],1,0)

print(table(df$Y, df$quartile))

if(!is.null(Wvars)){
  Wall<-subset(df, select=Wvars)
}


if(is.null(Wvars)){
  res<-as.numeric(washb_glm(Y=df$Y, tr=df$A, id=df$id, contrast=c(0,1), family=poisson(link='log'))$TR)
        
}else{
 
  Wscreen <- washb_prescreen(Y=df$Y,Ws=Wall,family="binomial", pval = 0.2 ,print=T)
  W<-subset(Wall, select=Wscreen)
  W<-design_matrix(W)
  #remove near zero variance columns
  preproc = caret::preProcess(W, method = c("zv", "nzv"))
  W = predict(preproc, W)
  W<-droplevels(W)
  d<-cbind(df$Y,df$A,W)
  d<-na.omit(d)


      res<-tmle(Y=d[,1],A=d[,2],W=d[,-c(1:2)], family="binomial",Q.SL.library=library,g.SL.library=library,id=df$id)
          }  

  return(res)
}



####################
#Unadjusted MH tests
####################

temp.unadj1v4<-temp.unadj2v4<-temp.unadj1v2<-temp.unadj1v3<-temp.unadj3v4<-matrix(NA, nrow=3, ncol=7)
for(j in 1:3){
  i<-j*7
  set.seed(12345)
  temp.unadj1v4[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,4),cutoff=4)
  temp.unadj2v4[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(2,4),cutoff=4)
  temp.unadj1v2[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,2),cutoff=4)
  temp.unadj1v3[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,3),cutoff=4)
  temp.unadj3v4[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(3,4),cutoff=4)
}
temp.unadj<-NULL
for(i in 1:3){
  temp.unadj<-rbind(temp.unadj,temp.unadj1v2[i,],temp.unadj1v3[i,],temp.unadj1v4[i,])
}
rownames(temp.unadj)<-c("lag7 1v2","lag7 1v3","lag7 1v4",
                        "lag14 1v2","lag14 1v3","lag14 1v4",
                        "lag21 1v2","lag21 1v3","lag21 1v4")
colnames(temp.unadj)<- c("PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
temp.unadj

save(Mean_lag7, Mean_lag14, Mean_lag21, temp.unadj,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/temp.Results.unadjusted.Rdata")

#Save dataset
save(d, 
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/temp.datasets.Rdata")







##############
#Effect of temperature on 
#hydrogen sulfide in the 
#stored drinking water
##############



d$Ydiar<-d$Y
d$Y<-d$H2S
H2S.unadj1v4<-H2S.unadj2v4<-H2S.unadj1v2<-H2S.unadj1v3<-H2S.unadj3v4<-matrix(NA, nrow=3, ncol=7)
for(j in 1:3){
  if(j==1){
    i<-1
  }else{
  i<-(j-1)*7
  }
  set.seed(6789)
  H2S.unadj1v4[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,4),cutoff=4)
  H2S.unadj2v4[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(2,4),cutoff=4)
  H2S.unadj1v2[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,2),cutoff=4)
  H2S.unadj1v3[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(1,3),cutoff=4)
  H2S.unadj3v4[j,]<-quart_glmner(d=d,Wvars=NULL,weather=paste0("temp.ave7.lag",i),q=c(3,4),cutoff=4)
}
H2S.unadj<-NULL
for(i in 1:3){
  H2S.unadj<-rbind(H2S.unadj,H2S.unadj1v2[i,],H2S.unadj1v3[i,],H2S.unadj1v4[i,])
}
rownames(H2S.unadj)<-c("lag7 1v2","lag7 1v3","lag7 1v4",
                        "lag14 1v2","lag14 1v3","lag14 1v4",
                        "lag21 1v2","lag21 1v3","lag21 1v4")
colnames(H2S.unadj)<- c("PR", "ci.lb", "ci.ub", "logPR", "se.logPR",  "Z", "p")
H2S.unadj
save(H2S.unadj, file="temp.h2s.Rdata")

d$Y<-d$Ydiar







####################
#Adjusted TMLE
####################

#Set library
library=c("SL.mean","SL.glm","SL.glmnet", "SL.gam", "SL.bayesglm")

temp.adj <- list()
for(j in 1:3){
  i<-j*7
  set.seed(438759)
  temp.adj1v4<-quart_glmner(d=d,Wvars=Wvars,weather=paste0("temp.ave7.lag",i),q=c(1,4),cutoff=4, library=library)
  temp.adj1v3<-quart_glmner(d=d,Wvars=Wvars,weather=paste0("temp.ave7.lag",i),q=c(1,3),cutoff=4, library=library)
  temp.adj1v2<-quart_glmner(d=d,Wvars=Wvars,weather=paste0("temp.ave7.lag",i),q=c(1,2),cutoff=4, library=library)
  temp.adj[[paste0("temp.ave7.lag",i,".1v4")]]<-temp.adj1v4
  temp.adj[[paste0("temp.ave7.lag",i,".1v3")]]<-temp.adj1v3
  temp.adj[[paste0("temp.ave7.lag",i,".1v2")]]<-temp.adj1v2
}





save(Mean_lag7, Mean_lag14, Mean_lag21, temp.unadj, temp.adj,
     file="C:/Users/andre/Dropbox/Trichy analysis/Results/temp.Results.adjusted.Rdata")



