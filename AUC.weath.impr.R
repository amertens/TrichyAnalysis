
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



#Functions:
getUnivariateR2 <- function(Y, psiHat.Pnv0, return.IC){
    n <- length(Y[,1])
    J <- ncol(Y)
    
    # compute Y bar for each outcome
    Ybar <- colMeans(Y)
    
    # MSE for each outcome
    MSE <- apply(matrix(1:J), 1, function(j){
        mean((Y[,j] - psiHat.Pnv0[,j])^2)
    })
    
    # ic for MSE for each outcome 
    IC.MSE <- apply(matrix(1:J), 1, function(j){
        (Y[,j] - psiHat.Pnv0[,j])^2 - MSE[j]
    })

    # Var for each otucome
    Var <- apply(matrix(1:J), 1, function(j){
        mean((Y[,j] - Ybar[j])^2)
    })
    
    # ic for Variance -- denominator
    IC.Var <- apply(matrix(1:J), 1, function(j){
        (Y[,j] - Ybar[j])^2 - Var[j]
    })
    
    # ic for log(numerator) - log(denominator)
    se.logR2 <- apply(matrix(1:J),1,function(j){
        grad <- matrix(c(1/MSE[j],-1/Var[j]),nrow=1)
        IC <- cbind(IC.MSE[,j],IC.Var[,j])
        sqrt(grad%*%t(IC)%*%IC%*%t(grad))/n
    })
    
    # point est + 95% CI
    out <- lapply(split(1:J,1:J),function(j){
        est <- 1 - MSE[j]/Var[j]
        ci.low <- 1 - exp(
            log(MSE[j]/Var[j]) + 1.96*se.logR2[j]
        )
        ci.high <- 1 - exp(
            log(MSE[j]/Var[j]) - 1.96*se.logR2[j]
        )
        return(c(est, ci.low, ci.high))
    })
    names(out) <- colnames(Y)
    out$MSE <- MSE
    out$Var <- Var
    out$IC <- NULL
    if(return.IC){
        out$IC <- list(IC.MSE = IC.MSE, IC.Var = IC.Var)
    }
    return(out)
}












###################
#Weather data processing
###################
head(weather)

#select weather vars
rainvars<-sapply(1:30, function(x) paste0("rain.",x))
maxtempvars<-sapply(1:30, function(x) paste0("maxt.",x))
mintempvars<-sapply(1:30, function(x) paste0("mint.",x))
humvars<-sapply(1:30, function(x) paste0("hum.",x))

#other to drop
Orainvars<-sapply(31:100, function(x) paste0("rain.",x))
Omaxtempvars<-sapply(31:100, function(x) paste0("maxt.",x))
Omintempvars<-sapply(31:100, function(x) paste0("mint.",x))
Ohumvars<-sapply(31:100, function(x) paste0("hum.",x))

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
#h2o.init(nthreads = 3)

#Set up library
h2o.glm.1 <- function(..., alpha = 0.0) h2o.glm.wrapper(..., alpha = alpha)
h2o.glm.2 <- function(..., alpha = 0.5) h2o.glm.wrapper(..., alpha = alpha)
h2o.glm.3 <- function(..., alpha = 1.0) h2o.glm.wrapper(..., alpha = alpha)
h2o.glm.4 <- function(..., lambda = 0.0) h2o.glm.wrapper(..., lambda = lambda)
h2o.glm.5 <- function(..., lambda = 1e-5) h2o.glm.wrapper(..., lambda = lambda)
h2o.glm.6 <- function(..., lambda = 1e-1) h2o.glm.wrapper(..., lambda = lambda)
h2o.randomForest.1 <- function(..., ntrees = 200, nbins = 50, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
#h2o.randomForest.2 <- function(..., ntrees = 200, sample_rate = 0.75, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
h2o.randomForest.3 <- function(..., ntrees = 200, sample_rate = 0.85, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, sample_rate = sample_rate, seed = seed)
#h2o.randomForest.4 <- function(..., ntrees = 200, nbins = 50, balance_classes = TRUE, seed = 1) h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, balance_classes = balance_classes, seed = seed)
h2o.gbm.1 <- function(..., ntrees = 100, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, seed = seed)
#h2o.gbm.2 <- function(..., ntrees = 100, nbins = 50, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
#h2o.gbm.3 <- function(..., ntrees = 100, max_depth = 10, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, max_depth = max_depth, seed = seed)
#h2o.gbm.4 <- function(..., ntrees = 100, col_sample_rate = 0.8, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
#h2o.gbm.5 <- function(..., ntrees = 100, col_sample_rate = 0.7, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
h2o.gbm.6 <- function(..., ntrees = 100, col_sample_rate = 0.6, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, col_sample_rate = col_sample_rate, seed = seed)
#h2o.gbm.7 <- function(..., ntrees = 100, balance_classes = TRUE, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, balance_classes = balance_classes, seed = seed)
h2o.gbm.8 <- function(..., ntrees = 100, max_depth = 3, seed = 1) h2o.gbm.wrapper(..., ntrees = ntrees, max_depth = max_depth, seed = seed)
#h2o.deeplearning.1 <- function(..., hidden = c(500,500), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.2 <- function(..., hidden = c(200,200,200), activation = "Tanh", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.3 <- function(..., hidden = c(500,500), activation = "RectifierWithDropout", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
h2o.deeplearning.4 <- function(..., hidden = c(500,500), activation = "Rectifier", epochs = 50, balance_classes = TRUE, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, balance_classes = balance_classes, seed = seed)
#h2o.deeplearning.5 <- function(..., hidden = c(100,100,100), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.6 <- function(..., hidden = c(50,50), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
#h2o.deeplearning.7 <- function(..., hidden = c(100,100), activation = "Rectifier", epochs = 50, seed = 1)  h2o.deeplearning

learners <- c("h2o.glm.wrapper","h2o.glm.1","h2o.glm.2","h2o.glm.3","h2o.glm.4","h2o.glm.5","h2o.glm.6",
             "h2o.randomForest.1", "h2o.randomForest.3",
             "h2o.gbm.1", "h2o.gbm.6", "h2o.gbm.8","h2o.deeplearning.4")

learners <- c("h2o.glm.wrapper","h2o.glm.1","h2o.glm.2")
learners <- c("h2o.glm.1","h2o.glm.2","h2o.glm.4","h2o.glm.6",
              "h2o.randomForest.3",
              "h2o.gbm.6")

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
  
  
  fit_bl = h2o.ensemble(x = x, y = y,  training_frame = train_frame,
                    family = "binomial",  learner = learners,
                    metalearner = metalearner, cvControl = list(V = cvControl),
                    parallel="multicore")
  

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

return(list("sl.AUC"=sl.AUC,"auc.V"=auc, "ensemble.AUC"=ensemble.AUC,"fit"=fit_bl))
}



################################
#Run model with only HH survey data
################################

#Data into H20 from R data.frame
d$Y<-as.factor(d$Y) #binary Y needs to be a factor for h20 ensemble

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate","stdywk","round1","round2","round3","round4","round5","round6","round7","round8","round9","round10","round11","maxtemp","mintemp","rain","hum0830","hum1730",rainvars,maxtempvars,mintempvars,humvars,Orainvars,Omaxtempvars,Omintempvars,Ohumvars))   #remove time variables

set.seed(12345)
BL.slfit<-cv.h2oEnsemble(d=d, flds=flds, V=V, y=y, x=x, learners = learners, metalearner = metalearner)
BL.slfit$sl.AUC
BL.slfit$ensemble.AUC



    
    
################################
#Load and fit full model with all rain data
################################

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate", humvars,maxtempvars,mintempvars,Orainvars,Omaxtempvars,Omintempvars,Ohumvars)) 


set.seed(12345)
rain.slfit<-cv.h2oEnsemble(d=d, flds=flds, V=V, y=y, x=x, learners = learners, metalearner = metalearner)
rain.slfit$sl.AUC


################################
#Load and fit full model with all temp  data
################################

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate", humvars,rainvars,Orainvars,Omaxtempvars,Omintempvars,Ohumvars)) 
X<-subset(d, select=x)
Y<-subset(d, select=y)
id<-subset(d, select=id)

temp.slfit<-cv.h2oEnsemble(d=d, flds=flds, V=V, y=y, x=x, learners = learners, metalearner = metalearner)
temp.slfit$sl.AUC
    
################################
#Load and fit full model with all temp and rain data
################################

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate", humvars,Orainvars,Omaxtempvars,Omintempvars,Ohumvars)) 
full.slfit<-cv.h2oEnsemble(d=d, flds=flds, V=V, y=y, x=x, learners = learners, metalearner = metalearner)
full.slfit$sl.AUC

################################
#Load and fit full model with only temp and rain data
################################

#Set up x and y variables for h2o model
y<-"Y"
DWeath<-subset(d, select=c("Y",rainvars,maxtempvars,mintempvars))
x<-setdiff(names(DWeath), c(y)) 


onlyWeath.slfit<-cv.h2oEnsemble(d=DWeath, flds=flds, V=V, y=y, x=x, learners = learners, metalearner = metalearner)
onlyWeath.slfit$sl.AUC


#Shut down cluster
h2o.shutdown(prompt = F)



    #plotting function
    library(scales)

    
uplot <- function(x,varlab="",svy,cols,yaxis=F, n) {

  # empty plot
  plot(1:8,1:8,type="n",
       xaxt="n",xlab="",xlim=c(0,6),
       yaxt="n",ylab="",bty="n",ylim=c(0.5,1)
  )

  # Y-axis
  if(yaxis==TRUE) mtext(seq(.5,1,by=.1),side=2,line=0,at=seq(0.5,1,by=0.1),las=1,col="gray20")
  segments(x0=0,x1=8,y0=seq(.5,1,by=0.1),col="gray80",lty=2)

  # X-axis labels
  mtext(c("Rain&Temp Only","Survey","+Rain","+Temp","+Rain&Temp"),side=1,line=0,at=seq(1,5,1),col=cols,cex=0.8,las=1)
  mtext(svy,side=1,line=2.5,col="gray20",cex=1)

  # plot points
  arrows(x0=1:n, y0= as.numeric(x[2,]), y1=as.numeric(x[3,]), col=cols,lwd=1,length=0.085,angle=90,code=3)
  points(1:n,x[1,],pch=21,cex=1.5,lwd=1,col=cols,bg="white")
  points(1:n,x[1,],pch=21,cex=1.5,lwd=0,col=cols,bg=alpha(cols,alpha=0.5))

}

# general label plot
ulabplot <- function(title) {
  plot(1,1,type="n",
       xaxt="n",xlab="",xlim=c(0,1),
       yaxt="n",ylab="",bty="n",ylim=c(0,1)
  )
  text(1,0.5,title,adj=1,cex=1.5)
}



#pdf("~/Dropbox/WBK-primary-analysis/results/figures/kenya-uptake.pdf",width=10.5,height=14)
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# cols <- c("gray30",cbPalette[c(2:4,5,6:8)])

black = "#000004FF"
blue = "#3366AA"
teal = "#11AA99"
green = "#66AA55"
chartr = "#CCCC55"
magent = "#992288"
red = "#EE3333"
orange = "#EEA722"
yellow = "#FFEE33"
grey = "#777777"
#cols=c(black,chartr,blue,teal,green,orange,red,magent)
cols=c(magent,black,blue,red,green)


BL.slfit$sl.AUC
rain.slfit$sl.AUC
temp.slfit$sl.AUC
full.slfit$sl.AUC


    auc.mat<-cbind(c(.62,.595,.645), BL.slfit$sl.AUC, rain.slfit$sl.AUC, temp.slfit$sl.AUC, full.slfit$sl.AUC)



uplot(x=as.data.frame(auc.mat[,1]),cols=cols,svy="AUC by predictor variable group",yaxis=T,n=1)
uplot(x=as.data.frame(auc.mat[,c(1:2)]),cols=cols,svy="AUC by predictor variable group",yaxis=T,n=2)
uplot(x=as.data.frame(auc.mat[,1:3]),cols=cols,svy="AUC by predictor variable group",yaxis=T,n=3)
uplot(x=as.data.frame(auc.mat[,1:4]),cols=cols,svy="AUC by predictor variable group",yaxis=T,n=4)
uplot(x=auc.mat,cols=cols,svy="AUC by predictor variable group",yaxis=T,n=5)

