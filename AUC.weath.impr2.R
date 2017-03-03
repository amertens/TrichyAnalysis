
rm(list=ls())

library(dplyr)
library(SuperLearner)
library(tmle)
library(washb)
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




getPredictionsOnValidation <- function(out, whichAlgorithm){
    Reduce("cbind",lapply(out, FUN=function(s){
        if(whichAlgorithm=="SuperLearner"){
            s$SL.predict
        }else if(whichAlgorithm=="discreteSuperLearner"){
            s$discreteSL.predict
        }else{
            s$library.predict[,which(
                s$libraryNames==whichAlgorithm
            )]
        }
    }))
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

#replace 99 with NAs
#for(i in 1:ncol(df)){
#  df[which(df[,i]==99),i]<-NA
#}

Y<-df$diar7d
id<-df$vilid
intdate<-df$intdate
Wfac<-subset(df,select= -c(vilid,hhid,individ,intdate,bdate,mdiar7d, diar2d, diar7d, diar14d, hcgi7d, diardays,  hw1, hw2, hw3, hw4, hw5, hw6, hw7, hw8, hw9, hw10, hw11, hw12, hws1, hws2, hws3, hws4, hws5, hws6, hws7, hws8, hws9, hws10, hws11, hws12))





#complete cases so design matrix works. XXX-Need to figure out where missingness is and fix
length(Y)
Y<-Y[complete.cases(Wfac)]
id<-id[complete.cases(Wfac)]
intdate<-intdate[complete.cases(Wfac)]
Wfac<-Wfac[complete.cases(Wfac),]
length(Y)

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






################################
#Run model with only HH survey data
################################

#remove other packages causing trouble:
try(detach(package:caret))


screen.corRank10<-function (Y, X, family, method = "pearson", rank = 10, ...) 
{
    listp <- apply(X, 2, function(x, Y, method) {
        ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
    }, Y = Y, method = method)
    whichVariable <- (rank(listp) <= rank)
    return(whichVariable)
}
screen.corRank4<-function (Y, X, family, method = "pearson", rank = 4, ...) 
{
    listp <- apply(X, 2, function(x, Y, method) {
        ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
    }, Y = Y, method = method)
    whichVariable <- (rank(listp) <= rank)
    return(whichVariable)
}


library<-c("SL.gam","SL.glmnet", "SL.glm", "SL.mean", "SL.bayesglm")

library = list(c("SL.glm", "screen.corP"),c("SL.glm", "screen.corRank"),c("SL.gam","screen.corP"),"SL.glmnet")
library = list(c("SL.glm", "screen.corP"),c("SL.glm", "screen.corRank10"),c("SL.glm.interaction","screen.corRank10"),"SL.glmnet")
                  
library = list(c("SL.glm", "screen.corP"),c("SL.glm", "screen.corRank10"),"SL.glmnet","SL.mean")
library = list(c("SL.glm", "screen.corP"),c("SL.glm", "screen.corRank10"),"SL.mean")


colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate","stdywk","round1","round2","round3","round4","round5","round6","round7","round8","round9","round10","round11","maxtemp","mintemp","rain","hum0830","hum1730",rainvars,maxtempvars,mintempvars,humvars))   #remove time variables
X<-subset(d, select=x)
Y<-subset(d, select=y)
id<-subset(d, select=id)

        #need to incorporate ID
        set.seed(12345)
        fit <- SuperLearner::CV.SuperLearner(
            Y = d$Y, X = X, SL.library = library, family = "binomial", V = 10, method = "method.AUC",
            verbose =T,parallel = 'multicore'
        )
   
    AUC<-mean(summary.CV.SuperLearner(fit)$Risk.SL)
      se<-summary.CV.SuperLearner(fit)$Table[1,3]
    AUC.ci<-c(AUC-1.96*se, AUC+1.96*se)
      
    blAUC<-list(auc=AUC, ci=AUC.ci, se=se)

    
    
################################
#Load and fit full model with all rain data
################################

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate", humvars,maxtempvars,mintempvars)) 
X<-subset(d, select=x)


#drop variables with >30% missing data
dim(X)
toomissing<-apply(X, 2, function(x) sum(is.na(x))/nrow(X)>0.3)
X<- subset(X, select=c(!toomissing)) 
dim(X)


#Add dummy variables for missing X:
        X.cont.dist = X
        xp = ncol(X.cont.dist)
        n.cont = nrow(X.cont.dist)
        sum_nas = apply(X.cont.dist, 2, function(x) sum(is.na(x)))
        nmesX = colnames(X.cont.dist)
        miss.cont = NULL
        nmesm = NULL
        for (k in 1:ncol(X)) {
            if (sum_nas[k] > 0) {
                ix = as.numeric(is.na(X.cont.dist[, k]))
                miss.cont = cbind(miss.cont, ix)
                nmesm = c(nmesm, paste("miss_", nmesX[k], sep = ""))
            }
        }
        colnames(miss.cont)=nmesm
        X<-cbind(X,miss.cont)

#KNN imputation for missing data, with dummy indicator variable added for missing data
impute_info = caret::preProcess(X, method = "medianImpute")
X = predict(impute_info, X)
rm(impute_info)


# Remove zero variance (constant) and near-zero-variance columns.
# This can help reduce overfitting and also helps us use a basic glm().
# However, there is a slight risk that we are discarding helpful information.
dim(X)
preproc = caret::preProcess(X, method = c("zv", "nzv"))
W = predict(preproc, X)
rm(preproc)
dim(X) # Review our dimensions.



        set.seed(12345)
        fitrain <- SuperLearner::CV.SuperLearner(
            Y = d$Y, X = X, SL.library = library, family = "binomial", V = 10, method = "method.AUC",
            verbose =T,parallel = 'multicore'
        )
   
    AUC<-mean(summary.CV.SuperLearner(fitrain)$Risk.SL)
      se<-summary.CV.SuperLearner(fitrain)$Table[1,3]
    AUC.ci<-c(AUC-1.96*se, AUC+1.96*se)
      
    rainAUC<-list(auc=AUC, ci=AUC.ci, se=se)
    rainAUC
    
################################
#Load and fit full model with all temp  data
################################

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate", humvars,rainvars)) 
X<-subset(d, select=x)
Y<-subset(d, select=y)
id<-subset(d, select=id)

#drop variables with >30% missing data
dim(X)
toomissing<-apply(X, 2, function(x) sum(is.na(x))/nrow(X)>0.3)
X<- subset(X, select=c(!toomissing)) 
dim(X)


#Add dummy variables for missing X:
        X.cont.dist = X
        xp = ncol(X.cont.dist)
        n.cont = nrow(X.cont.dist)
        sum_nas = apply(X.cont.dist, 2, function(x) sum(is.na(x)))
        nmesX = colnames(X.cont.dist)
        miss.cont = NULL
        nmesm = NULL
        for (k in 1:ncol(X)) {
            if (sum_nas[k] > 0) {
                ix = as.numeric(is.na(X.cont.dist[, k]))
                miss.cont = cbind(miss.cont, ix)
                nmesm = c(nmesm, paste("miss_", nmesX[k], sep = ""))
            }
        }
        colnames(miss.cont)=nmesm
        X<-cbind(X,miss.cont)

#KNN imputation for missing data, with dummy indicator variable added for missing data
impute_info = caret::preProcess(X, method = "medianImpute")
X = predict(impute_info, X)
rm(impute_info)


# Remove zero variance (constant) and near-zero-variance columns.
# This can help reduce overfitting and also helps us use a basic glm().
# However, there is a slight risk that we are discarding helpful information.
dim(X)
preproc = caret::preProcess(X, method = c("zv", "nzv"))
W = predict(preproc, X)
rm(preproc)
dim(X) # Review our dimensions.



        set.seed(12345)
        fittemp <- SuperLearner::CV.SuperLearner(
            Y = d$Y, X = X, SL.library = library, family = "binomial", V = 10, method = "method.AUC",
            verbose =T,parallel = 'multicore'
        )
   
    AUC<-mean(summary.CV.SuperLearner(fittemp)$Risk.SL)
      se<-summary.CV.SuperLearner(fittemp)$Table[1,3]
    AUC.ci<-c(AUC-1.96*se, AUC+1.96*se)
      
    tempAUC<-list(auc=AUC, ci=AUC.ci, se=se)
    
    
################################
#Load and fit full model with all temp and rain data
################################

#Set up x and y variables for h2o model
colnames(d)
y<-"Y"
x<-setdiff(names(d), c(y, "intdate", humvars)) 
X<-subset(d, select=x)
Y<-subset(d, select=y)
id<-subset(d, select=id)

#drop variables with >30% missing data
dim(X)
toomissing<-apply(X, 2, function(x) sum(is.na(x))/nrow(X)>0.3)
X<- subset(X, select=c(!toomissing)) 
dim(X)


#Add dummy variables for missing X:
        X.cont.dist = X
        xp = ncol(X.cont.dist)
        n.cont = nrow(X.cont.dist)
        sum_nas = apply(X.cont.dist, 2, function(x) sum(is.na(x)))
        nmesX = colnames(X.cont.dist)
        miss.cont = NULL
        nmesm = NULL
        for (k in 1:ncol(X)) {
            if (sum_nas[k] > 0) {
                ix = as.numeric(is.na(X.cont.dist[, k]))
                miss.cont = cbind(miss.cont, ix)
                nmesm = c(nmesm, paste("miss_", nmesX[k], sep = ""))
            }
        }
        colnames(miss.cont)=nmesm
        X<-cbind(X,miss.cont)

#KNN imputation for missing data, with dummy indicator variable added for missing data
impute_info = caret::preProcess(X, method = "medianImpute")
X = predict(impute_info, X)
rm(impute_info)


# Remove zero variance (constant) and near-zero-variance columns.
# This can help reduce overfitting and also helps us use a basic glm().
# However, there is a slight risk that we are discarding helpful information.
dim(X)
preproc = caret::preProcess(X, method = c("zv", "nzv"))
W = predict(preproc, X)
rm(preproc)
dim(X) # Review our dimensions.



        set.seed(12345)
        fitfull <- SuperLearner::CV.SuperLearner(
            Y = d$Y, X = X, SL.library = library, family = "binomial", V = 10, method = "method.AUC",
            verbose =T,parallel = 'multicore'
        )
   
    AUC<-mean(summary.CV.SuperLearner(fitfull)$Risk.SL)
      se<-summary.CV.SuperLearner(fitfull)$Table[1,3]
    AUC.ci<-c(AUC-1.96*se, AUC+1.96*se)
      
    fullAUC<-list(auc=AUC, ci=AUC.ci, se=se)
    fullAUC
    
    
    
    
################################
#Load and fit full model with only temp and rain data
################################

#Set up x and y variables for h2o model
y<-"Y"
X<-subset(d, select=c(rainvars,maxtempvars,mintempvars))
Y<-subset(d, select=y)
id<-subset(d, select=id)

#drop variables with >30% missing data
dim(X)
toomissing<-apply(X, 2, function(x) sum(is.na(x))/nrow(X)>0.3)
X<- subset(X, select=c(!toomissing)) 
dim(X)


#Add dummy variables for missing X:
        X.cont.dist = X
        xp = ncol(X.cont.dist)
        n.cont = nrow(X.cont.dist)
        sum_nas = apply(X.cont.dist, 2, function(x) sum(is.na(x)))
        nmesX = colnames(X.cont.dist)
        miss.cont = NULL
        nmesm = NULL
        for (k in 1:ncol(X)) {
            if (sum_nas[k] > 0) {
                ix = as.numeric(is.na(X.cont.dist[, k]))
                miss.cont = cbind(miss.cont, ix)
                nmesm = c(nmesm, paste("miss_", nmesX[k], sep = ""))
            }
        }
        colnames(miss.cont)=nmesm
        X<-cbind(X,miss.cont)

#KNN imputation for missing data, with dummy indicator variable added for missing data
impute_info = caret::preProcess(X, method = "medianImpute")
X = predict(impute_info, X)
rm(impute_info)


# Remove zero variance (constant) and near-zero-variance columns.
# This can help reduce overfitting and also helps us use a basic glm().
# However, there is a slight risk that we are discarding helpful information.
dim(X)
preproc = caret::preProcess(X, method = c("zv", "nzv"))
W = predict(preproc, X)
rm(preproc)
dim(X) # Review our dimensions.



        set.seed(12345)
        fitOnlyWeath <- SuperLearner::CV.SuperLearner(
            Y = d$Y, X = X, SL.library = library, family = "binomial", V = 10, method = "method.AUC",
            verbose =T,parallel = 'multicore'
        )
   
    AUC<-mean(summary.CV.SuperLearner(fitOnlyWeath)$Risk.SL)
      se<-summary.CV.SuperLearner(fitOnlyWeath)$Table[1,3]
    AUC.ci<-c(AUC-1.96*se, AUC+1.96*se)
      
    onlyWeathAUC<-list(auc=AUC, ci=AUC.ci, se=se)
    onlyWeathAUC    
    
    
    
    
    
    
    

    #Assess if AUC confidence intervals are significantly diff:
    jointci.lb<-(blAUC$auc-fullAUC$auc)-1.96*sqrt(blAUC$se^2-fullAUC$se^2)
    jointci.ub<-(blAUC$auc-fullAUC$auc)+1.96*sqrt(blAUC$se^2-fullAUC$se^2)
    jointci.lb
    jointci.ub
    
    
    
    summary.CV.SuperLearner(fit)
    summary.CV.SuperLearner(fitrain)
    summary.CV.SuperLearner(fittemp)
    summary.CV.SuperLearner(fitfull)
    
     

    #Calc normal curve based 95% CI
    me <- qnorm(.975)*(sd(summary(fit)$Risk.SL)/sqrt(10)) 
    lb<-mean(summary(fit)$Risk.SL)-me
    ub<-mean(summary(fit)$Risk.SL)+me
    bl.auc<-c(mean(summary(fit)$Risk.SL),lb,ub)
    
    me <- qnorm(.975)*(sd(summary(fitrain)$Risk.SL)/sqrt(10)) 
    lb<-mean(summary(fitrain)$Risk.SL)-me
    ub<-mean(summary(fitrain)$Risk.SL)+me
    rain.auc<-c(mean(summary(fitrain)$Risk.SL),lb,ub)
    
    me <- qnorm(.975)*(sd(summary(fittemp)$Risk.SL)/sqrt(10)) 
    lb<-mean(summary(fittemp)$Risk.SL)-me
    ub<-mean(summary(fittemp)$Risk.SL)+me
    temp.auc<-c(mean(summary(fittemp)$Risk.SL),lb,ub)
    
    me <- qnorm(.975)*(sd(summary(fitfull)$Risk.SL)/sqrt(10)) 
    lb<-mean(summary(fitfull)$Risk.SL)-me
    ub<-mean(summary(fitfull)$Risk.SL)+me
    full.auc<-c(mean(summary(fitfull)$Risk.SL),lb,ub)
    
    me <- qnorm(.975)*(sd(summary(fitOnlyWeath)$Risk.SL)/sqrt(10)) 
    lb<-mean(summary(fitOnlyWeath)$Risk.SL)-me
    ub<-mean(summary(fitOnlyWeath)$Risk.SL)+me
    Wonly.auc<-c(mean(summary(fitOnlyWeath)$Risk.SL),lb,ub)
    
    bl.auc
    rain.auc
    temp.auc
    full.auc
    Wonly.auc
    
    auc.mat<-cbind(Wonly.auc, bl.auc, rain.auc, temp.auc, full.auc)
    
    
    
    
    
    
    
    
    #plotting function
    library(scales)

    
uplot <- function(x,varlab="",svy,cols,yaxis=F, n) {

  # empty plot
  plot(1:8,1:8,type="n",
       xaxt="n",xlab="",xlim=c(0,6),
       yaxt="n",ylab="",bty="n",ylim=c(0.5,1)
  )

  # Y-axis
  if(yaxis==TRUE) mtext(seq(0,100,by=20),side=2,line=0,at=seq(0.5,1,by=0.1),las=1,col="gray20")
  segments(x0=0,x1=8,y0=seq(0,1,by=0.2),col="gray80",lty=2)

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


uplot(x=as.data.frame(auc.mat[,1]),cols=cols,svy="AUC by predictor variable group",yaxis=T,n=1)
uplot(x=as.data.frame(auc.mat[,c(1:2)]),cols=cols,svy="AUC by predictor variable group",yaxis=T,n=2)
uplot(x=as.data.frame(auc.mat[,1:3]),cols=cols,svy="AUC by predictor variable group",yaxis=T,n=3)
uplot(x=as.data.frame(auc.mat[,1:4]),cols=cols,svy="AUC by predictor variable group",yaxis=T,n=4)
uplot(x=auc.mat,cols=cols,svy="AUC by predictor variable group",yaxis=T,n=5)
