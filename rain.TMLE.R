
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
  levels(rain[,3+i])[which(levels(rain[,3+i])=="Trace")]<-"0.01"
}

#change rain vars to numeric:
rain<-apply(rain, 2, function(x) as.numeric(as.character(x))) 

#

ave.window<-function(i, winsize, var, data){
  col<-which(colnames(data) %in% paste0(var,".",i))
  windowmean<-NA
  try(windowmean<-apply(data[,col:(col-1+winsize)],1,function(x) mean(x)))
  return(windowmean)
}

averain<-matrix(data=NA,nrow=nrow(rain),ncol=28)
ave.names<-rep(NA, 28)
for(i in 1:28){
  averain[,i]<-ave.window(i, 7, "rain", rain)
  ave.names[i]<- paste0("rain.ave7.lag",i)
}
averain<-cbind(rain[,1:3],averain)
colnames(averain)[4:ncol(averain)]<-ave.names
averain<-as.data.frame(averain)

#Create long term rain average
averain$LT8<-ave.window(15, 60, "rain", rain) #set i to 7 days after the associated 7-day lagged variable
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
Y<-Y[complete.cases(Wfac)]
id<-id[complete.cases(Wfac)]
intdate<-intdate[complete.cases(Wfac)]
Wfac<-Wfac[complete.cases(Wfac),]


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
#H20 analysis
##########################

#merge suvery and weather
dim(survey)
dim(averain)
d<-merge(survey,averain, by="intdate", all.x = F, all.y = F)  #Do I need to impute missing weather data 
#for the dates towards the end of the study?
dim(d)

colnames(d)


#Mean diarrhea by quartile:
d%>%
  group_by(quartile = ntile(rain.ave7.lag1, 4))%>%
  summarize(quart.n=n(),meanrain=mean(rain.ave7.lag1),meandiar=mean(Y))
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

quart_tmle<-function(d,W=NULL,weather="rain.ave7.lag1",q=c(1,4),cutoff=4){
colnames(d)[which( colnames(d)==paste0(weather))]="weather"
df<-d%>%
  mutate(quartile = ntile(weather, cutoff))%>%
  subset(quartile==q[1]|quartile==q[2])%>%
  mutate(A=ifelse(quartile==q[2],1,0))

if(is.null(W)){
W<-matrix(1:nrow(df),nrow=nrow(df),ncol=1)
W[]<-1 #Make empty covariate matrix    
  }  

  res<-tmle(Y=df$Y,A=df$A,W=W, family="binomial",Q.SL.library=c("SL.glm"),g.SL.library=c("SL.glm"),id=df$id)
  return(res)
}


rain.unadj1v4<-rain.unadj2v4<-rain.unadj1v2<-rain.unadj1v3<-rain.unadj3v4<-matrix(NA, nrow=28, ncol=3)
for(i in 1:28){
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

### NOTE:
#time lag 1-7 have temporality issues due to the 7 day diarrheal disease recall. Use lag=8, 15, and 22 for
#1,2, and 3 week lag for adjusted analysis

################
#Heavy rain events
################
W<-matrix(1:nrow(d),nrow=nrow(d),ncol=1)
W[]<-1 #Make empty covariate matrix   
Hrain.unadj<-matrix(NA, nrow=28, ncol=3)
for(i in 1:28){
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


#create a tmle function that fits Q and A with H20
h2o_tmle<-function(d,A="rain.ave7.lag8",q=c(1,4),cutoff=4,learners, metalearner, longtermlag=NULL){
  
j=1
if(!is.null(longtermlag)){
    col<-which(colnames(d) %in% paste0("LT",longtermlag))
    d$LTtertile = ntile(d[,col], 3)
  j=3
  d.full<-d
  Astore<-A
}

for(i in 1:j){
  
  if(!is.null(longtermlag)){
     d<-d.full
     A<-Astore
     d<-subset(d, LTtertile==i)
      d<-subset(d, select= -LTtertile)
      cat("Tertile:",i," Y(0,1):",table(d$Y),"\n")
    }
  
d$Acont<-d[,which( colnames(d)==paste0(A))]
df<-d%>%
  mutate(quartile = ntile(Acont, cutoff))%>%
  subset(quartile==q[1]|quartile==q[2])%>%
  mutate(A=ifelse(quartile==q[2],1,0))%>%
  select(-(rain.ave7.lag1:rain.ave7.lag28))%>%
  select(-c(Acont,quartile,id,intdate))

cat("\nA(0,1):",table(d$A),"\n")

#remove additional variables with no variance in the subset contrast 
pre<-ncol(df)
preproc = caret::preProcess(df, method = c("zv"))
df = predict(preproc, df)
rm(preproc)
cat("\nNum variables removed due to zero var:",pre-ncol(df),"\n")


#Data into H20 from R data.frame
df$Y<-as.factor(df$Y) #binary Y needs to be a factor for h20 ensemble
df$A<-as.factor(df$A)
data<-as.h2o(df)
y<-"Y"
A<-"A"

x<-setdiff(names(data), y)

#check that outcome is binary, and table count of outcome
h2o.unique(data[y])
  h2o.table(data[y])
  
  n <- nrow(data)  # Total number of training samples
  h2o.table(data[y])['Count']/n
  
#################
#superlearner
#################

#Run superlearner model:
  fit_Q = h2o.ensemble(x = x, y = y,  training_frame = data,
                    family = "binomial",  learner = learners,
                    metalearner = metalearner, cvControl = list(V = 5),
                    parallel="multicore")
  fit_Q$runtime 

 
  #save model
  #h2o.save_ensemble()   #comment out until full model is fit.
  
  
  #Predict from model (test case with A=primwatPrivateWell):
  newdata1<-data[c(y,x)]
  newdata1["A"]<-1
  Q1<-predict(fit_Q, newdata1)
  Q1
  
  newdata0<-data[c(y,x)]
  newdata0["A"]<-0
  Q0<-predict(fit_Q, newdata0)
  Q0
  
  Q<-cbind(as.data.frame(Q0$pred)[,3],as.data.frame(Q1$pred)[,3])
  
  
  #Fit g(A|W)
w<-setdiff(names(data), c(y,A))
  
  #Run superlearner model:
  fit_g = h2o.ensemble(x = w, y = A,  training_frame = data,
                    family = "binomial",  learner = learners,
                    metalearner = metalearner, cvControl = list(V = 5),
                    parallel="multicore")
  fit_g$runtime 
  
  g<-predict(fit_g, newdata=data)
  g<-as.data.frame(g$pred)[,3]
  
  
  W<-df[,w]
  df$Y<-as.numeric(df$Y)-1 #convert back to binary dummy variables from factors for H20
  df$A<-as.numeric(df$A)-1
  tmlefit<-tmle(Y=df$Y, A=df$A,W=W,Q=Q,g1W=g, family="binomial", df$id)
    if(!is.null(longtermlag)){
      if(i==1) tmlefit1<-tmlefit
      if(i==2) tmlefit2<-tmlefit
      if(i==3) tmlefit3<-tmlefit
    }
}

    if(is.null(longtermlag)){return(list("tmlefit"=tmlefit, "predQ"=Q,"predg"=g,"fit_Q"=fit_Q,"fit_g"=fit_g))}
    if(!is.null(longtermlag)){return(list("tmle.tert1"=tmlefit1, "tmle.tert2"=tmlefit2, "tmle.tert3"=tmlefit3))}
}


#set up learner library and metalearners
  #install.packages("https://h2o-release.s3.amazonaws.com/h2o-ensemble/R/h2oEnsemble_0.1.8.tar.gz", repos = NULL)
      # Used to calculate test set AUC (requires version >=1.0.1 of cvAUC)
  
  #set learners
  learners = c("h2o.glm.wrapper", "h2o.randomForest.wrapper",
             "h2o.gbm.wrapper")
  
  #custom learners
    h2o.glm.1 <- function(..., alpha = 0.0) h2o.glm.wrapper(..., alpha = alpha)
    h2o.glm.2 <- function(..., alpha = 0.5) h2o.glm.wrapper(..., alpha = alpha)
    h2o.glm.3 <- function(..., alpha = 1.0) h2o.glm.wrapper(..., alpha = alpha)
      #learners = c("h2o.glm.wrapper","h2o.glm.1","h2o.glm.2","h2o.glm.3")
    
h2o.glm.1 <- function(..., alpha = 0.0) h2o.glm.wrapper(..., alpha = alpha)
h2o.glm.2 <- function(..., alpha = 0.5) h2o.glm.wrapper(..., alpha = alpha)
h2o.glm.3 <- function(..., alpha = 1.0) h2o.glm.wrapper(..., alpha = alpha)
h2o.randomForest.1 <- function(..., ntrees = 200, nbins = 50, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
h2o.randomForest.2 <- function(..., ntrees = 200, sample_rate = 0.75, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
h2o.randomForest.3 <- function(..., ntrees = 200, sample_rate = 0.85, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
h2o.randomForest.4 <- function(..., ntrees = 200, nbins = 50, balance_classes = TRUE, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, balance_classes = balance_classes, seed = seed)
h2o.gbm.1 <- function(..., ntrees = 100, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, seed = seed)
h2o.gbm.2 <- function(..., ntrees = 100, nbins = 50, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
h2o.gbm.3 <- function(..., ntrees = 100, max_depth = 10, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, max_depth = max_depth, seed = seed)
h2o.gbm.4 <- function(..., ntrees = 100, col_sample_rate = 0.8, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
h2o.gbm.5 <- function(..., ntrees = 100, col_sample_rate = 0.7, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
h2o.gbm.6 <- function(..., ntrees = 100, col_sample_rate = 0.6, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
h2o.gbm.7 <- function(..., ntrees = 100, balance_classes = TRUE, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, balance_classes = balance_classes, seed = seed)
h2o.gbm.8 <- function(..., ntrees = 100, max_depth = 3, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, max_depth = max_depth, seed = seed)
h2o.deeplearning.1 <- function(..., hidden = c(500,500), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
h2o.deeplearning.2 <- function(..., hidden = c(200,200,200), activation = "Tanh", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
h2o.deeplearning.3 <- function(..., hidden = c(500,500), activation = "RectifierWithDropout", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
h2o.deeplearning.4 <- function(..., hidden = c(500,500), activation = "Rectifier", epochs = 50, balance_classes = TRUE, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, balance_classes = balance_classes, seed = seed)
h2o.deeplearning.5 <- function(..., hidden = c(100,100,100), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
h2o.deeplearning.6 <- function(..., hidden = c(50,50), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
h2o.deeplearning.7 <- function(..., hidden = c(100,100), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning

learners <- c("h2o.glm.wrapper",
             "h2o.randomForest.1", "h2o.randomForest.2",
             "h2o.gbm.1", "h2o.gbm.6", "h2o.gbm.8",
             "h2o.deeplearning.1", "h2o.deeplearning.6", "h2o.deeplearning.7")

#Metalearning
  #Historically, methods such as GLM or non-negative least squares (NNLS) have been used to find the optimal 
  #weighted combination of the base learners, however any supervised learning algorithm can be used as a metalearner. 
  #To use a GLM with non-negative weights, you simply pass  non_negative = TRUE  to the generic,  h2o.glm.wrapper  
  #function as follows:
  h2o.glm_nn <- function(..., non_negative = TRUE) {
    h2o.glm.wrapper(..., non_negative = non_negative)
  }
  metalearner <- "h2o.glm_nn"

  
  

  #set up h20 connection
# Clean slate - just in case the cluster was already running.
h2o.removeAll()

# Start an H2O cluster with nthreads = num cores on your machine.
# -1 means to use all cores.
h2o.init(nthreads = -1)
  
###Run adjusted TMLE functions

set.seed(84753)
test<-h2o_tmle(d=d,A="rain.ave7.lag8",q=c(1,4),cutoff=4,learners=learners, metalearner=metalearner)
test

set.seed(84753)
test<-h2o_tmle(d=d,A="rain.ave7.lag15",q=c(1,4),cutoff=4,learners=learners, metalearner=metalearner)
test

set.seed(84753)
test<-h2o_tmle(d=d,A="rain.ave7.lag22",q=c(1,4),cutoff=4,learners=learners, metalearner=metalearner)
test

set.seed(84753)
test<-h2o_tmle(d=d,A="rain.ave7.lag8",q=c(2,4),cutoff=4,learners=learners, metalearner=metalearner)
test


#test out long term stratification:
test<-h2o_tmle(d=d,A="rain.ave7.lag8",q=c(1,4),cutoff=4,learners=learners, metalearner=metalearner,longtermlag=8)
test

test<-h2o_tmle(d=d,A="rain.ave7.lag15",q=c(1,4),cutoff=4,learners=learners, metalearner=metalearner,longtermlag=15)
test

  h2o.shutdown(prompt = F) #shut down cluster
