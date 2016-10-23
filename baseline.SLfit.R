

library(dplyr)
library(h2o)
library(SuperLearner)
library(tmle)
library(tmle.npvi)
library(washb)
library(h2oEnsemble)
library(cvAUC) 

load("prescreened_blcov.Rdata")

#subset to complete cases
dim(df)
df<-df[complete.cases(W),]
dim(df)

#transform W to all binary indicators for SL model:

W<- df[,colnames(df) %in% names.W]
W<-design_matrix(W)

data<-cbind(df$diar7d,W)
colnames(data)[1]<-"Y"


  #Set 10% of data as validation
  data$split<-rbinom(n=nrow(data),size=1,prob=0.1)


#set up h20 connection
localH20<-h2o.init(nthreads=-1,max_mem_size="8G")

#Use this code to check connection if already connected:
h2o.init()
#If not connected elsewhere, will set up defaults on your machine


#Data into H20 from R data.frame
data<-as.h2o(data)
y<-"Y"
x<-setdiff(names(data),y)

fit<-h2o.gbm(x=x,y=y,training_frame=data)

pred<-h20.predict(fit, test)

plot(fit)




#check that outcome is binary
h2o.unique(data[y])


#We can query the categorical "levels" as well 
  h2o.levels(data[y])

  
  h2o.nacnt(data[y])
  
  h2o.nacnt(data)
  
  h2o.table(data[y])
  
  
  n <- nrow(data)  # Total number of training samples
  h2o.table(data[y])['Count']/n
  
  
  #Split data into train, validate, and test set
  #train <- data[data['split']=="train",]
  #nrow(train)
  
  valid <- data[data['split']==1,]
  nrow(valid)
  
  #test <- data[data['split']=="test",]
  #nrow(test)
  

  x <- setdiff(names(data), c("Y","split"))  #Remove the split and Y columns
  x

#Now that we have specified x and y, we can train the GBM model using a few 
#non-default model parameters. Since we are predicting a binary response, 
#we set distribution = "bernoulli".

  
  
  #Cross-validated model
  cvmodel <- h2o.gbm(x = x, y = y,
                     training_frame = data,
                     validation_frame = valid,
                     #distribution = "bernoulli",
                     ntrees = 100,
                     max_depth = 4,
                     learn_rate = 0.1,
                     nfolds = 5)
    print(cvmodel)
  
  perf <- h2o.performance(model = cvmodel, newdata = valid)
  class(perf)
  
  h2o.r2(perf)
  h2o.auc(perf)
  h2o.mse(perf)
  
  
  print(h2o.auc(cvmodel, train = TRUE))
  print(h2o.auc(cvmodel, xval = TRUE))
  
  #Grid search over different parameters
  ntrees_opt <- c(5,50,100)
  max_depth_opt <- c(2,3,5)
  learn_rate_opt <- c(0.1,0.2)
  
  hyper_params = list('ntrees' = ntrees_opt,
                      'max_depth' = max_depth_opt,
                      'learn_rate' = learn_rate_opt)
  
  gs <- h2o.grid(algorithm = "gbm", 
                 grid_id = "eeg_demo_gbm_grid",
                 hyper_params = hyper_params,
                 x = x, y = y, 
                 training_frame = data, 
                 validation_frame = valid,
                 nfolds=5)
  
  
  print(gs)
  
  
  # print out the auc for all of the models
  auc_table <- h2o.getGrid(grid_id = "eeg_demo_gbm_grid", sort_by = "auc", decreasing = TRUE)
  print(auc_table)
  
  best_model <- h2o.getModel(auc_table@model_ids[[1]])
  h2o.auc(best_model, valid = TRUE)  #Validation AUC for best model
  
  best_perf <- h2o.performance(model = best_model, newdata = test)
  h2o.auc(best_perf)
  
  
  #Try out superlearner
  #install.packages("https://h2o-release.s3.amazonaws.com/h2o-ensemble/R/h2oEnsemble_0.1.8.tar.gz", repos = NULL)
# Used to calculate test set AUC (requires version >=1.0.1 of cvAUC)
  
#Metalearning
  #Historically, methods such as GLM or non-negative least squares (NNLS) have been used to find the optimal 
  #weighted combination of the base learners, however any supervised learning algorithm can be used as a metalearner. 
  #To use a GLM with non-negative weights, you simply pass  non_negative = TRUE  to the generic,  h2o.glm.wrapper  
  #function as follows:
  h2o.glm_nn <- function(..., non_negative = TRUE) {
    h2o.glm.wrapper(..., non_negative = non_negative)
  }
  metalearner <- "h2o.glm_nn"
  
  
#Run superlearner model:
  super<-h2o.ensemble(x = x, y = y,
               training_frame = train,
               validation_frame = valid,
               family = "binomial")
  
  pp <- predict(super, test) 
  predictions <- as.data.frame(pp$pred)[,3] #third column, p1 is P(Y==1) labels <- as.data.frame(test[,y])[,1]
  labels <- as.data.frame(test[,y])[,1]
  cvAUC::AUC(predictions = predictions, labels = labels) 
  
#Using trichy data
  load("TrichyAnalysis_Clean.Rdata")
  set.seed(123456)
  villages <- levels(as.factor(data[,1]))
  trn <- sample(villages, length(villages)*0.5)
  
  df.one <- subset(data, data[,1] %in% trn)
  df.two <- subset(data, !(data[,1] %in% trn))
  
  #check good split of diarrhea cases
  table(df.one[,31])
  table(df.two[,31])
  #Looks ok
  
  #Run exposure definition searching with df. one
  #Using "data" as dataframe variable to preserve previous code
  #save unsplit dataset:
  data.whole<-data
  data<-df.one
  data <- h2o.importFile(data)
  data = h2o.importFile(path = normalizePath("../data/covtype.full.csv"))
  
  #Save full covariates matrix
  covariates<-X
  
  #Create covariate matrix for cv.superlearner
  X<-data[,c(8,16,21,38,52,53,56,60,128:131)]
  #Create missing data indicators for binary variables:
  X[,13:18]=0
  X[,13][is.na(X[,2])]=1
  X[,14][is.na(X[,3])]=1
  X[,15][is.na(X[,4])]=1
  X[,16][is.na(X[,5])]=1
  X[,17][is.na(X[,6])]=1
  X[,18][is.na(X[,7])]=1
  X[,c(2:7)][is.na(X[,c(2:7)])]=0
  colnames(X)[13:18]<-c("momworkNA","hwsoapNA","latsoapNA","fecessmNA","chobs_handNA","chobs_clothNA")
  
  #Select outcome data
  Y<-as.matrix(data$diar7d)
  
  #superlearner
  super<-h2o.ensemble(x = X, y = Y, training_frame = X,
                      cvControl = list(V = 10, shuffle = TRUE),
                      family = "binomial")