
rm(list=ls())

library(dplyr)
library(SuperLearner)
library(tmle)
library(tmle.npvi)
library(washb)
library(h2o)
library(h2oEnsemble)
library(cvAUC) 
library(caret)
library(metafor)


setwd("C:/Users/andre/Documents/Trichy analysis")
load("prescreened_blcov.Rdata")
load("LaggedWeather.Rdata")



###################
#Weather data processing
###################
head(weather)

rainvars<-sapply(1:100, function(x) paste0("rain.",x))
rain<-subset(weather, select=c("year","month","day",rainvars))
#Change "trace" to 0.01 in rainfall variables
for(i in 1:100){
  levels(rain[,3+i])[which(levels(rain[,3+i])=="Trace")]<-"0.1"
}

#change rain vars to numeric:
rain<-apply(rain, 2, function(x) as.numeric(as.character(x))) 


ave.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) mean(x)))
  return(windowmean)
}

averain<-matrix(data=NA,nrow=nrow(rain),ncol=15)
ave.names<-rep(NA, 15)
for(i in 7:21){
  averain[,i-6]<-ave.window(i, 7, "rain", rain)
  ave.names[i-6]<- paste0("rain.ave7.lag",i)
}
averain<-cbind(rain[,1:3],averain)
colnames(averain)[4:ncol(averain)]<-ave.names
averain<-as.data.frame(averain)

#Create long term rain average
averain$LT8<-ave.window(15, 60, "rain", rain) #set i to 7 days after the associated 7-day lagged variable
averain$LT15<-ave.window(22, 60, "rain", rain)
averain$LT22<-ave.window(29, 60, "rain", rain)

averain$LT8<-ave.window(15, 60, "rain", rain) #set 3-level category
averain$LT15<-ave.window(22, 60, "rain", rain)
averain$LT22<-ave.window(29, 60, "rain", rain)

#Set to date format for merge
averain$intdate<-as.Date(paste(averain$month,averain$day,averain$year, sep="/"),"%m/%d/%Y")
averain<-subset(averain, select=-c(year,month,day))



###################
#Survey data processing
###################

# Expand out factors into indicators (ignore ordinality of several factors).
# Alternatively we could convert the ordinal factors to numerics.
dim(df)
df<-df[,-5] #remove duplicate "round" variable
head(df)

Y<-df$diar7d
id<-df$vilid
intdate<-df$intdate
Wfac<-subset(df,select= -c(vilid,hhid,individ,intdate,bdate,mdiar7d, diar2d, diar7d, diar14d, hcgi7d, diardays ))



#complete cases so design matrix works. XXX-Need to figure out where missingness is and fix
#Y<-Y[complete.cases(Wfac)]
#id<-id[complete.cases(Wfac)]
#intdate<-intdate[complete.cases(Wfac)]
#Wfac<-Wfac[complete.cases(Wfac),]
table(is.na(Wfac))

W<-design_matrix(Wfac)
dim(W)

# Remove zero variance (constant) and near-zero-variance columns.
# This can help reduce overfitting and also helps us use a basic glm().
# However, there is a slight risk that we are discarding helpful information.
preproc = caret::preProcess(W, method = c("zv", "nzv"))
W = predict(preproc, W)
rm(preproc)
dim(W) # Review our dimensions.

survey<-cbind(intdate,Y,id,W)



##########################
#Unadjusted analysis
##########################

#merge suvery and weather
dim(survey)
dim(averain)
d<-merge(survey,averain, by="intdate", all.x = F, all.y = F)  #Do I need to impute missing weather data 
#for the dates towards the end of the study?
dim(d)

colnames(d)


#Mean diarrhea by quartile:

d_rain<-subset(d, rain.ave7.lag7>0)
low_hi_cut<-quantile(d_rain$rain.ave7.lag7, prob = seq(0, 1, length = 4), type = 5, na.rm=TRUE)[2:3]

d$raincat.7<-d$raincat.14<-d$raincat.21<-"No rain"
d$raincat.7[d$rain.ave7.lag7>0]<-"Low rain"
d$raincat.14[d$rain.ave7.lag14>0]<-"Low rain"
d$raincat.21[d$rain.ave7.lag21>0]<-"Low rain"
d$raincat.7[d$rain.ave7.lag7>low_hi_cut[1]]<-"Med rain"
d$raincat.14[d$rain.ave7.lag14>low_hi_cut[1]]<-"Med rain"
d$raincat.21[d$rain.ave7.lag21>low_hi_cut[1]]<-"Med rain"
d$raincat.7[d$rain.ave7.lag7>low_hi_cut[2]]<-"High rain"
d$raincat.14[d$rain.ave7.lag14>low_hi_cut[2]]<-"High rain"
d$raincat.21[d$rain.ave7.lag21>low_hi_cut[2]]<-"High rain"
d$raincat.7<-factor(d$raincat.7, levels=c("No rain", "Low rain", "Med rain", "High rain"))
d$raincat.14<-factor(d$raincat.14, levels=c("No rain", "Low rain", "Med rain", "High rain"))
d$raincat.21<-factor(d$raincat.21, levels=c("No rain", "Low rain", "Med rain", "High rain"))



#########################
#Summary tables.
#########################

d%>%
  group_by(raincat.7)%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meandiar=mean(Y), Ndiar=sum(Y))
d%>%
  group_by(raincat.14)%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meandiar=mean(Y), Ndiar=sum(Y))
d%>%
  group_by(raincat.21)%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meandiar=mean(Y), Ndiar=sum(Y))


d%>%
  group_by(raincat.7, quartile = ntile(LT8, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meandiar=mean(Y), Ndiar=sum(Y))
d%>%
  group_by(raincat.14, quartile = ntile(LT15, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meandiar=mean(Y), Ndiar=sum(Y))
d%>%
  group_by(raincat.21, quartile = ntile(LT22, 3))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meandiar=mean(Y), Ndiar=sum(Y))



d%>%
  group_by(quartile = ntile(rain.ave7.lag7, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag7),meandiar=mean(Y))
d%>%
  group_by(quartile = ntile(rain.ave7.lag14, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag14),meandiar=mean(Y))
d%>%
  group_by(quartile = ntile(rain.ave7.lag21, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag21),meandiar=mean(Y))




####################
#Unadjusted MH tests
####################

rain_tmle<-function(d,W=NULL,weather="rain.ave7.lag7",q=c(1,4),cutoff=4, library=NULL){
colnames(d)[which( colnames(d)==paste0(weather))]="weather"
df<-d%>%
  mutate(quartile = ntile(weather, cutoff))%>%
  subset(quartile==q[1]|quartile==q[2])%>%
  mutate(A=ifelse(quartile==q[2],1,0))

if(is.null(W)){
  mhtab <- table(df$A, df$Y)
    mhtab <- mhtab[c(2:1), c(2:1)]
    mh.est <- rma.mh(ai = mhtab[1, 1], bi = mhtab[1, 2], ci = mhtab[2, 1], di = mhtab[2,2], measure = "RR", drop00 = TRUE, level = 95) 
        res <- c(exp(mh.est$b), exp(mh.est$ci.lb), exp(mh.est$ci.ub), 
            mh.est$b, mh.est$se, mh.est$zval, mh.est$pval)
        names(res) <- c("PR,", "ci.lb", "ci.ub", "logPR", "se.logPR", 
            "Z", "p")
}else{
      res<-tmle(Y=df$Y,A=df$A,W=W, family="binomial",Q.SL.library=library,g.SL.library=library,id=df$id)
  }  

  return(res)
}

rain.unadj1v4<-rain.unadj2v4<-rain.unadj1v2<-rain.unadj1v3<-rain.unadj3v4<-matrix(NA, nrow=3, ncol=7)
for(j in 1:3){
  i<-j*7
  set.seed(12345)
  rain.unadj1v4[j,]<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,4),cutoff=4)
  rain.unadj2v4[j,]<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(2,4),cutoff=4)
  rain.unadj1v2[j,]<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,2),cutoff=4)
  rain.unadj1v3[j,]<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,3),cutoff=4)
  rain.unadj3v4[j,]<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(3,4),cutoff=4)
}


rain.unadj1v2
rain.unadj1v3
rain.unadj1v4

rain.unadj1v4
rain.unadj2v4
rain.unadj3v4



################
#Heavy rain events
################
W<-matrix(1:nrow(d),nrow=nrow(d),ncol=1)
W[]<-1 #Make empty covariate matrix   
Hrain.unadj<-matrix(NA, nrow=28, ncol=3)
for(i in 1:28){
  set.seed(12345)
 colnames(d)[which( colnames(d)==paste0("rain.ave7.lag",i))]="weather"
df<-d%>%
  mutate(quartile = ntile(weather, 10))%>%
  mutate(A=ifelse(quartile==10|quartile==9,1,0)) # >80% as heavy rainfall event 

  fit<-tmle(Y=df$Y,A=df$A,W=W, family="binomial",Q.SL.library=c("SL.glm"),g.SL.library=c("SL.glm"),id=df$id)
  Hrain.unadj[i,1]<-fit$estimates$RR$psi
  Hrain.unadj[i,c(2,3)]<-fit$estimates$RR$CI
  
  colnames(d)[which( colnames(d)=="weather")]=paste0("rain.ave7.lag",i)
}
Hrain.unadj



##############
#Adjusted analysis
##############























##############
#Supplimentary figures
##############

#Sensitivity analysis to lag time

rain.unadj1v4<-rain.unadj2v4<-rain.unadj1v2<-rain.unadj1v3<-rain.unadj3v4<-matrix(NA, nrow=28, ncol=3)
for(i in 1:28){
  set.seed(12345)
  temp<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,4),cutoff=4)
  rain.unadj1v4[i,1]<-temp$estimates$RR$psi
  rain.unadj1v4[i,c(2,3)]<-temp$estimates$RR$CI
  
  temp<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(2,4),cutoff=4)
  rain.unadj2v4[i,1]<-temp$estimates$RR$psi
  rain.unadj2v4[i,c(2,3)]<-temp$estimates$RR$CI  
  
  temp<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,2),cutoff=4)
  rain.unadj1v2[i,1]<-temp$estimates$RR$psi
  rain.unadj1v2[i,c(2,3)]<-temp$estimates$RR$CI
  
  temp<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(1,3),cutoff=4)
  rain.unadj1v3[i,1]<-temp$estimates$RR$psi
  rain.unadj1v3[i,c(2,3)]<-temp$estimates$RR$CI
  
  temp<-quart_tmle(d=d,W=NULL,weather=paste0("rain.ave7.lag",i),q=c(3,4),cutoff=4)
  rain.unadj3v4[i,1]<-temp$estimates$RR$psi
  rain.unadj3v4[i,c(2,3)]<-temp$estimates$RR$CI
}

rain.unadj1v2
rain.unadj1v3
rain.unadj1v4

rain.unadj1v4
rain.unadj2v4
rain.unadj3v4

test<-glm(d$Y~d$rain.ave7.lag13, family=poisson(link=`log`))
summary(test)
test<-glm(d$Y~d$rain.ave7.lag14, family=poisson(link=`log`))
summary(test)
test<-glm(d$Y~d$rain.ave7.lag15, family=poisson(link=`log`))
summary(test)

plotdf<-data.frame(rep(1:nrow(rain.unadj1v4)),rain.unadj1v4)
colnames(plotdf)<-c("x","RD","lo","up")
ggplot(data = plotdf) + 
  geom_point(mapping = aes(x=x, y=RD), size=3) +
  geom_errorbar(mapping = aes(x=x, y=RD, ymin=lo, ymax=up), size=1) +
  geom_hline(yintercept=1) +
  ggtitle("Q4 v Q1 weekly rain average rainfall") +
  labs(x="Days lagged",y="Risk Ratio of Diarrheal Disease") + 
  theme(plot.title=element_text(""), 
        axis.text.x = element_text(color="#666666", face="bold",angle = 45, hjust = 1)) +
  theme(plot.title = element_text(color="#666666", face="bold", size=22, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=12)) +
  theme_bw()

plotdf<-data.frame(rep(1:nrow(rain.unadj2v4)),rain.unadj2v4)
colnames(plotdf)<-c("x","RD","lo","up")
ggplot(data = plotdf) + 
  geom_point(mapping = aes(x=x, y=RD), size=3) +
  geom_errorbar(mapping = aes(x=x, y=RD, ymin=lo, ymax=up), size=1) +
  geom_hline(yintercept=1) +
  ggtitle("Q4 v Q2 weekly rain average rainf") +
  labs(x="Days lagged",y="Risk Ratio of Diarrheal Disease") + 
  theme(plot.title=element_text(""), 
        axis.text.x = element_text(color="#666666", face="bold",angle = 45, hjust = 1)) +
  theme(plot.title = element_text(color="#666666", face="bold", size=22, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=12)) +
  theme_bw()

plotdf<-data.frame(rep(1:nrow(rain.unadj3v4)),rain.unadj1v4)
colnames(plotdf)<-c("x","RD","lo","up")
ggplot(data = plotdf) + 
  geom_point(mapping = aes(x=x, y=RD), size=3) +
  geom_errorbar(mapping = aes(x=x, y=RD, ymin=lo, ymax=up), size=1) +
  geom_hline(yintercept=1) +
  ggtitle("Q4 v Q3 weekly rain average rainf") +
  labs(x="Days lagged",y="Risk Ratio of Diarrheal Disease") + 
  theme(plot.title=element_text(""), 
        axis.text.x = element_text(color="#666666", face="bold",angle = 45, hjust = 1)) +
  theme(plot.title = element_text(color="#666666", face="bold", size=22, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=12)) +
  theme_bw()

plotdf<-data.frame(rep(1:nrow(rain.unadj2v4)),rain.unadj2v4)
colnames(plotdf)<-c("x","RD","lo","up")
ggplot(data = plotdf) + 
  geom_point(mapping = aes(x=x, y=RD), size=3) +
  geom_errorbar(mapping = aes(x=x, y=RD, ymin=lo, ymax=up), size=1) +
  geom_hline(yintercept=1) +
  ggtitle("Q4 v Q3 weekly rain average rainf") +
  labs(x="Days lagged",y="Risk Ratio of Diarrheal Disease") + 
  theme(plot.title=element_text(""), 
        axis.text.x = element_text(color="#666666", face="bold",angle = 45, hjust = 1)) +
  theme(plot.title = element_text(color="#666666", face="bold", size=22, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=12)) +
  theme_bw()



################
#Heavy rain events
################
W<-matrix(1:nrow(d),nrow=nrow(d),ncol=1)
W[]<-1 #Make empty covariate matrix   
Hrain.unadj<-matrix(NA, nrow=28, ncol=3)
for(i in 1:28){
  set.seed(12345)
 colnames(d)[which( colnames(d)==paste0("rain.ave7.lag",i))]="weather"
df<-d%>%
  mutate(quartile = ntile(weather, 10))%>%
  mutate(A=ifelse(quartile==10|quartile==9,1,0)) # >80% as heavy rainfall event 

  fit<-tmle(Y=df$Y,A=df$A,W=W, family="binomial",Q.SL.library=c("SL.glm"),g.SL.library=c("SL.glm"),id=df$id)
  Hrain.unadj[i,1]<-fit$estimates$RR$psi
  Hrain.unadj[i,c(2,3)]<-fit$estimates$RR$CI
  
  colnames(d)[which( colnames(d)=="weather")]=paste0("rain.ave7.lag",i)
}
Hrain.unadj
