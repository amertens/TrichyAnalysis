###Rain curves




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
library(ggplot2)
library(gridExtra)


setwd("C:/Users/andre/Documents/Trichy analysis")
load("prescreened_blcov.Rdata")
load("LaggedWeather.Rdata")



###################
#Weather data processing
###################
head(weather)

#select weather vars
rainvars<-sapply(1:100, function(x) paste0("rain.",x))
maxtempvars<-sapply(1:100, function(x) paste0("maxt.",x))
mintempvars<-sapply(1:100, function(x) paste0("mint.",x))
humvars<-sapply(1:100, function(x) paste0("hum.",x))

weather<-subset(weather, select=c("year","month","day","maxtemp","mintemp","rain","hum0830","hum1730",rainvars,maxtempvars,mintempvars,humvars))


#Change "trace" to 0.01 in rainfall variables
for(i in 1:100){
  levels(weather[,8+i])[which(levels(weather[,8+i])=="Trace")]<-"0.01"
}

#change rain vars to numeric:
weather<-apply(weather, 2, function(x) as.numeric(as.character(x))) 


#Set to date format for merge
weather<-as.data.frame(weather)
weather$intdate<-as.Date(paste(weather$month,weather$day,weather$year, sep="/"),"%m/%d/%Y")
weather<-subset(weather, select=-c(year,month,day))



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

#Bind baseline data together
survey<-cbind(intdate,Y,id,W)


#merge suvery and weather
dim(survey)
dim(W)
d<-merge(survey,weather, by="intdate", all.x = F, all.y = F)  #Do I need to impute missing weather data?
#for the dates towards the end of the study?
dim(d)



#set up h20 clusters
h2o.removeAll()
h2o.init(nthreads = -1)

#Data into H20 from R data.frame
d$Y<-as.factor(d$Y) #binary Y needs to be a factor for h20 ensemble




#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate"))


#Metalearning
  h2o.glm_nn <- function(..., non_negative = TRUE) {
    h2o.glm.wrapper(..., non_negative = non_negative)
  }
  metalearner <- "h2o.glm_nn"


  train_frame<-as.h2o(d)

  #check outcome distribution
  print(h2o.table(train_frame[y]))

  
 

#Random Grid Search
#Demo of random grid search, all from Erin LeDell's code. This will be pretty slow. First we fit the base models. 

# Random Grid Search (e.g. 120 second maximum)
# This is set to run fairly quickly, increase max_runtime_secs 
# or max_models to cover more of the hyperparameter space.
# Also, you can expand the hyperparameter space of each of the 
# algorithms by modifying the hyper param code below.

search_criteria <- list(strategy = "RandomDiscrete", 
                        max_runtime_secs = 120)
V<-10

# GBM Hyperparamters
learn_rate_opt <- c(0.01, 0.03) 
max_depth_opt <- c(3, 4, 5, 6, 9)
sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
col_sample_rate_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
hyper_params <- list(learn_rate = learn_rate_opt,
                     max_depth = max_depth_opt, 
                     sample_rate = sample_rate_opt,
                     col_sample_rate = col_sample_rate_opt)

gbm_grid <- h2o.grid("gbm", x = x, y = y,
                     training_frame = train_frame,
                     ntrees = 100,
                     seed = 1,
                     nfolds = V,
                     fold_assignment = "Modulo",
                     keep_cross_validation_predictions = TRUE,
                     hyper_params = hyper_params,
                     search_criteria = search_criteria)
gbm_models <- lapply(gbm_grid@model_ids, function(model_id) h2o.getModel(model_id))



# RF Hyperparamters
mtries_opt <- 8:20 
max_depth_opt <- c(5, 10, 15, 20, 25)
sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
col_sample_rate_per_tree_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
hyper_params <- list(mtries = mtries_opt,
                     max_depth = max_depth_opt,
                     sample_rate = sample_rate_opt,
                     col_sample_rate_per_tree = col_sample_rate_per_tree_opt)

rf_grid <- h2o.grid("randomForest", x = x, y = y,
                    training_frame = train_frame,
                    ntrees = 200,
                    seed = 1,
                    nfolds = V,
                    fold_assignment = "Modulo",
                    keep_cross_validation_predictions = TRUE,                    
                    hyper_params = hyper_params,
                    search_criteria = search_criteria)
rf_models <- lapply(rf_grid@model_ids, function(model_id) h2o.getModel(model_id))

# GLM Hyperparamters
alpha_opt <- seq(0,1,0.1)
lambda_opt <- c(0,1e-7,1e-5,1e-3,1e-1)
hyper_params <- list(alpha = alpha_opt,
                     lambda = lambda_opt)

glm_grid <- h2o.grid("glm", x = x, y = y,
                     training_frame = train_frame,
                     family = "binomial",
                     nfolds = V,
                     fold_assignment = "Modulo",
                     keep_cross_validation_predictions = TRUE,                    
                     hyper_params = hyper_params,
                     search_criteria = search_criteria)
glm_models <- lapply(glm_grid@model_ids, function(model_id) h2o.getModel(model_id))

# Create a list of all the base models
models <- c(gbm_models, rf_models, glm_models)



#Combine the base models using nnls metalearner.
# Let's stack!
SLfit <- h2o.stack(models = models, 
                   response_frame = train_frame[, y],
                   metalearner = metalearner, seed=248357)



  
  
  
  

# Generate predictions on the test set
pred0 <- predict.h2o.ensemble(SLfit, train_frame)
predictions <- as.data.frame(pred0$pred)[,c("p1")]  #p1 is P(Y==1)
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.df<-c(0,me-ci, me, me+ci)

pred_frame<-train_frame
index<-which(colnames(pred_frame) %in% "rain.1")

for(i in 1:20){
pred_frame[,index:(index+99)]<-pred_frame[,index:(index+99)]+0.5
  
pred <- predict.h2o.ensemble(SLfit, pred_frame)
predictions <- as.data.frame(pred$pred)[,c("p1")] 
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.i<-c(0.5*i, me-ci, me, me+ci)

pred.df<-rbind(pred.df, pred.i)
}
rownames(pred.df)<-NULL
pred.df<-as.data.frame(pred.df)
colnames(pred.df)<-c("addrain", "l.ci", "mean", "u.ci")
pred.df



#Getting adjusted  curves
p1<-ggplot(data=pred.df,aes(addrain, mean)) +
  geom_jitter(mapping=aes(x=addrain, y=mean), shape=1, size=3, height=0.01, width=0.1) + 
  #How to add confidence intervals to the predictions and curve?
  geom_smooth(method='loess', span=2)+ coord_cartesian(ylim=c(0,0.05), xlim=c(0,10))  +xlab("Added Rainfall (mm)")+ ylab("Mean Predicted Diarrheal disease risk")+
  ggtitle("Predicted Diarrheal Risk \n Under Simulated Added Rain")+ theme(axis.text=element_text(size=12), title=element_text(size=12))

p1


########################
#Temperature predictions
########################

predictions <- as.data.frame(pred0$pred)[,c("p1")]  
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.df<-c(0,me-ci, me, me+ci)

pred_frame<-train_frame
index1<-which(colnames(pred_frame) %in% "maxt.1")
index2<-which(colnames(pred_frame) %in% "mint.1")

for(i in 1:20){
pred_frame[,index1:(index1+99)]<-pred_frame[,index1:(index1+99)]+0.25
pred_frame[,index2:(index2+99)]<-pred_frame[,index2:(index2+99)]+0.25
  
pred <- predict.h2o.ensemble(SLfit, pred_frame)
predictions <- as.data.frame(pred$pred)[,c("p1")] 
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.i<-c(0.25*i, me-ci, me, me+ci)

pred.df<-rbind(pred.df, pred.i)
}
rownames(pred.df)<-NULL
pred.df<-as.data.frame(pred.df)
colnames(pred.df)<-c("addT", "l.ci", "mean", "u.ci")
pred.df

#Getting adjusted  curves
p2<-ggplot(data=pred.df,aes(addT, mean)) +
  geom_jitter(mapping=aes(x=addT, y=mean), shape=1, size=3, height=0.01, width=0.1) + 
  #How to add confidence intervals to the predictions and curve?
  geom_smooth(method='loess', span=2)+ coord_cartesian(ylim=c(0,0.1), xlim=c(0,5))  +xlab("Added Temp (C)")+ ylab("Mean Predicted Diarrheal disease risk")+
  ggtitle("Predicted Diarrheal Risk \n Under Simulated Added Temp")+ theme(axis.text=element_text(size=12), title=element_text(size=12))

p2





########################
#Heavy rain predictions
########################
predictions <- as.data.frame(pred0$pred)[,c("p1")]  #p1 is P(Y==1)
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.df<-c(0,me-ci, me, me+ci)

pred_frame<-train_frame
index<-which(colnames(pred_frame) %in% "rain.1")

for(i in 1:20){
pred_frame[,index+6+i]<-pred_frame[,index+6+i]+25
  
pred <- predict.h2o.ensemble(SLfit, pred_frame)
predictions <- as.data.frame(pred$pred)[,c("p1")] 
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.i<-c(0.5*i, me-ci, me, me+ci)

pred.df<-rbind(pred.df, pred.i)
}
rownames(pred.df)<-NULL
pred.df<-as.data.frame(pred.df)
colnames(pred.df)<-c("addHrainday", "l.ci", "mean", "u.ci")
pred.df



#Getting adjusted  curves
p3<-ggplot(data=pred.df,aes(addHrainday, mean)) +
  geom_jitter(mapping=aes(x=addHrainday, y=mean), shape=1, size=3, height=0.01, width=0.1) + 
  #How to add confidence intervals to the predictions and curve?
  geom_smooth(method='loess', span=2)+ coord_cartesian(ylim=c(0,0.05), xlim=c(0,10))  +xlab("Added Rainfall (mm)")+ ylab("Mean Predicted Diarrheal disease risk")+
  ggtitle("Predicted Diarrheal Risk Under \n Under Simulated Added Heavy Rain Days")+ theme(axis.text=element_text(size=12), title=element_text(size=12))

p3


#Examining predicted diarrhea under no rain:
predictions <- as.data.frame(pred0$pred)[,c("p1")]  #p1 is P(Y==1)
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.df<-c(0,me-ci, me, me+ci)

pred_frame<-train_frame
index<-which(colnames(pred_frame) %in% "rain.1")

pred_frame[,index:(index+99)]<-0

pred <- predict.h2o.ensemble(SLfit, pred_frame)
predictions <- as.data.frame(pred$pred)[,c("p1")] 
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.i<-c(0, me-ci, me, me+ci)

pred.df
pred.i


########################
#Humidity predictions
########################
predictions <- as.data.frame(pred0$pred)[,c("p1")]  #p1 is P(Y==1)
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.df<-c(0,me-ci, me, me+ci)

pred_frame<-train_frame
index<-which(colnames(pred_frame) %in% "hum.1")

for(i in 1:20){
pred_frame[,index:(index+99)]<-pred_frame[,index:(index+99)]+1
  
pred <- predict.h2o.ensemble(SLfit, pred_frame)
predictions <- as.data.frame(pred$pred)[,c("p1")] 
me<-mean(predictions)
ci <- qnorm(.975)*(sd(predictions)/sqrt(length(predictions))) 
pred.i<-c(1*i, me-ci, me, me+ci)

pred.df<-rbind(pred.df, pred.i)
}
rownames(pred.df)<-NULL
pred.df<-as.data.frame(pred.df)
colnames(pred.df)<-c("addhum", "l.ci", "mean", "u.ci")
pred.df



#Getting adjusted  curves
p4<-ggplot(data=pred.df,aes(addhum, mean)) +
  geom_jitter(mapping=aes(x=addhum, y=mean), shape=1, size=3, height=0.01, width=0.1) + 
  #How to add confidence intervals to the predictions and curve?
  geom_smooth(method='loess', span=2)+ coord_cartesian(ylim=c(0,0.05), xlim=c(0,10))  +xlab("Added Humidity (%)")+ ylab("Mean Predicted Diarrheal disease risk")+
  ggtitle("Predicted Diarrheal Risk Under \n Simulated Added Humididty")+ theme(axis.text=element_text(size=12), title=element_text(size=12))

p4

#Shut down cluster
h2o.shutdown(prompt = F)

















