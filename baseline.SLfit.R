
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
#h2o.init(nthreads = -1)
h2o.init(nthreads = 2)

#Set up library
#h2o.glm.1 <- function(..., alpha = 0.0) h2o.glm.wrapper(..., alpha = alpha)
#h2o.glm.2 <- function(..., alpha = 0.5) h2o.glm.wrapper(..., alpha = alpha)
#h2o.glm.3 <- function(..., alpha = 1.0) h2o.glm.wrapper(..., alpha = alpha)
#h2o.glm.4 <- function(..., lambda = 0.0) h2o.glm.wrapper(..., lambda = lambda)
#h2o.glm.5 <- function(..., lambda = 1e-5) h2o.glm.wrapper(..., lambda = lambda)
#h2o.glm.6 <- function(..., lambda = 1e-1) h2o.glm.wrapper(..., lambda = lambda)
#h2o.randomForest.1 <- function(..., ntrees = 200, nbins = 50, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
#h2o.randomForest.2 <- function(..., ntrees = 200, sample_rate = 0.75, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
#h2o.randomForest.3 <- function(..., ntrees = 200, sample_rate = 0.85, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
#h2o.randomForest.4 <- function(..., ntrees = 200, nbins = 50, balance_classes = TRUE, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, balance_classes = balance_classes, seed = seed)
#h2o.gbm.1 <- function(..., ntrees = 100, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, seed = seed)
#h2o.gbm.2 <- function(..., ntrees = 100, nbins = 50, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
#h2o.gbm.3 <- function(..., ntrees = 100, max_depth = 10, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, max_depth = max_depth, seed = seed)
#h2o.gbm.4 <- function(..., ntrees = 100, col_sample_rate = 0.8, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
#h2o.gbm.5 <- function(..., ntrees = 100, col_sample_rate = 0.7, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
#h2o.gbm.6 <- function(..., ntrees = 100, col_sample_rate = 0.6, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
#h2o.gbm.7 <- function(..., ntrees = 100, balance_classes = TRUE, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, balance_classes = balance_classes, seed = seed)
#h2o.gbm.8 <- function(..., ntrees = 100, max_depth = 3, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, max_depth = max_depth, seed = seed)
#h2o.deeplearning.1 <- function(..., hidden = c(500,500), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.2 <- function(..., hidden = c(200,200,200), activation = "Tanh", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.3 <- function(..., hidden = c(500,500), activation = "RectifierWithDropout", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.4 <- function(..., hidden = c(500,500), activation = "Rectifier", epochs = 50, balance_classes = TRUE, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, balance_classes = balance_classes, seed = seed)
#h2o.deeplearning.5 <- function(..., hidden = c(100,100,100), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.6 <- function(..., hidden = c(50,50), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.7 <- function(..., hidden = c(100,100), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning

#learners <- c("h2o.glm.wrapper",
#             "h2o.randomForest.1", "h2o.randomForest.2",
#             "h2o.gbm.1", "h2o.gbm.6", "h2o.gbm.8",
#             "h2o.deeplearning.1", "h2o.deeplearning.6", "h2o.deeplearning.7")


#Metalearning
  h2o.glm_nn <- function(..., non_negative = TRUE) {
    h2o.glm.wrapper(..., non_negative = non_negative)
  }
  metalearner <- "h2o.glm_nn"



#Set up Random Grid Search hyperparameters
#Demo of random grid search, all from Erin LeDell's code. This will be pretty slow. First we fit the base models. 

# Random Grid Search (e.g. 120 second maximum)
# This is set to run fairly quickly, increase max_runtime_secs 
# or max_models to cover more of the hyperparameter space.
# Also, you can expand the hyperparameter space of each of the 
# algorithms by modifying the hyper param code below.

search_criteria <- list(strategy = "RandomDiscrete", 
                        max_runtime_secs = 300)

# GBM Hyperparamters
learn_rate_opt <- c(0.01, 0.03) 
max_depth_opt <- c(3, 4, 5, 6, 9)
sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
col_sample_rate_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
hyper_params.GBM <- list(learn_rate = learn_rate_opt,
                     max_depth = max_depth_opt, 
                     sample_rate = sample_rate_opt,
                     col_sample_rate = col_sample_rate_opt)

# RF Hyperparameters for grid search 
mtries_opt <- 8:20 
max_depth_opt <- c(5, 10, 15, 20, 25)
sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
col_sample_rate_per_tree_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
hyper_params.RF <- list(mtries = mtries_opt,
                     max_depth = max_depth_opt,
                     sample_rate = sample_rate_opt,
                     col_sample_rate_per_tree = col_sample_rate_per_tree_opt)

# GLM Hyperparamters
alpha_opt <- seq(0,1,0.1)
lambda_opt <- c(0,1e-7,1e-5,1e-3,1e-1)
hyper_params.GLM <- list(alpha = alpha_opt,
                     lambda = lambda_opt)



#Run CV superlearner model:

  # Divide into 10 folds
V=10 
seed<-23189 #seed for h20 functions
set.seed(14623)
flds <- createFolds(d$Y, k = V, list = TRUE, returnTrain = FALSE)

cv.h2oEnsemble<-function(d, flds, V, y, x, learners, metalearner, cvControl=5){
for(i in 1:V){
  
  cat("\n\nXXXXXXXXXXXXXXXXXX\nCross-validated fold left out: ",i,"\nXXXXXXXXXXXXXXXXXX\n\n")
  
  #create empty dfs to hold results
  if(i==1){
    cv.AUC<-rep(NA,V)
    ensemble.AUC<-NULL
  }
  
  #split v-fold into test and train
  test<-d[flds[[i]],]
  train<-d[!(1:nrow(d) %in% flds[[i]]),]
  
  holdout_frame<-as.h2o(test)
  train_frame<-as.h2o(train)

  #check outcome distribution
  print(h2o.table(train_frame[y]))
  print(h2o.table(holdout_frame[y]))
  
  
  #fit_bl = h2o.ensemble(x = x, y = y,  training_frame = train_frame,
  #                  family = "binomial",  learner = learners,
  #                  metalearner = metalearner, cvControl = list(V = cvControl),
  #                  parallel="multicore")
  
  
gbm_grid <- h2o.grid("gbm", x = x, y = y,
                     training_frame = train_frame,
                     ntrees = 100,
                     seed = 1,
                     nfolds = 10,
                     fold_assignment = "Modulo",
                     keep_cross_validation_predictions = TRUE,
                     hyper_params = hyper_params.GBM,
                     search_criteria = search_criteria)
gbm_models <- lapply(gbm_grid@model_ids, function(model_id) h2o.getModel(model_id))


rf_grid <- h2o.grid("randomForest", x = x, y = y,
                    training_frame = train_frame,
                    ntrees = 200,
                    seed = 1,
                    nfolds = 10,
                    fold_assignment = "Modulo",
                    keep_cross_validation_predictions = TRUE,                    
                    hyper_params = hyper_params.RF,
                    search_criteria = search_criteria)
rf_models <- lapply(rf_grid@model_ids, function(model_id) h2o.getModel(model_id))


glm_grid <- h2o.grid("glm", x = x, y = y,
                     training_frame = train_frame,
                     family = "binomial",
                     nfolds = 10,
                     fold_assignment = "Modulo",
                     keep_cross_validation_predictions = TRUE,                    
                     hyper_params = hyper_params.GLM,
                     search_criteria = search_criteria)
glm_models <- lapply(glm_grid@model_ids, function(model_id) h2o.getModel(model_id))


#Deep1<-h2o.deeplearning(x=x, y=y, training_frame = train_frame, activation ="Rectifier",
#                        hidden = c(50, 50), epochs = 50, seed = 35489, nfolds=V, fold_assignment = "Modulo")
#deep_model <-  h2o.getModel(Deep1@model_id)

#Deep2<-h2o.deeplearning(x=x, y=y, training_frame = train_frame, activation ="Tanh", nfolds=V, fold_assignment = "Modulo")

# Create a list of all the base models
models <- c(gbm_models, rf_models, glm_models)



#Combine the base models using nnls metalearner.
fit_bl <- h2o.stack(models = models, 
                   response_frame = train_frame[, y],
                   metalearner = metalearner, seed=248357)
  
  
  #cat("Fold ",i,"runtime: \n") 
  #print(fit_bl$runtime)

# Generate predictions on the test set
pred <- predict.h2o.ensemble(fit_bl, holdout_frame)
predictions <- as.data.frame(pred$pred)[,c("p1")]  #p1 is P(Y==1)
labels <- as.data.frame(holdout_frame[,c(y)])[,1]



# Ensemble test AUC 
cv.AUC[i]<-cvAUC::AUC(predictions = predictions , labels = labels)



# Base learner test AUC (for comparison)
L <- length(learners)
auc <- sapply(seq(L), function(l) AUC(predictions = as.data.frame(pred$basepred)[,l], labels = labels)) 
ensemble.AUC<-rbind(ensemble.AUC,data.frame(learners, auc,i))
}
  
#Calc normal curve based 95% CI
me <- qnorm(.975)*(sd(cv.AUC)/sqrt(V)) 
mean(cv.AUC)-me
mean(cv.AUC)+me

sl.AUC=c(mean(cv.AUC),mean(cv.AUC)-me,mean(cv.AUC)+me)
names(sl.AUC)<-c("Mean.AUC","95.lb","95.ub")

return(list("sl.AUC"=sl.AUC,"auc.V"=auc, "ensemble.AUC"=ensemble.AUC))
}



################################
#Run model with only HH survey data
################################

#Data into H20 from R data.frame
d$Y<-as.factor(d$Y) #binary Y needs to be a factor for h20 ensemble

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate","stdywk","round1","round2","round3","round4","round5","round6","round7","round8","round9","round10","round11","maxtemp","mintemp","rain","hum0830","hum1730",rainvars,maxtempvars,mintempvars,humvars))   #remove time variables

BL.slfit<-cv.h2oEnsemble(d=d, flds=flds, V=V, y=y, x=x, learners = learners, metalearner = metalearner)
BL.slfit$sl.AUC
BL.slfit$ensemble.AUC



################################
#Load and fit full model with all weather data
################################

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate")) 

Weather.slfit<-cv.h2oEnsemble(d=d, flds=flds, V=V, y=y, x=x, learners = learners, metalearner = metalearner)
Weather.slfit$sl.AUC
Weather.slfit$ensemble.AUC




#Shut down cluster
h2o.shutdown(prompt = F)




