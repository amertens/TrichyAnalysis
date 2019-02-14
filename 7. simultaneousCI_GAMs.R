
rm(list=ls())
library(devtools)
#install_github("ben-arnold/tmleAb")
library(tmleAb)
library(washb)
library(dplyr)
library(scales)
library(SuperLearner)
    try(detach(package:mgcv))
    library(gam)



rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
}



GAM_simulCI<-function (Y, Age, W = NULL, id = NULL, SL.library = c( "SL.gam"), cvControl = list(V=5), 
    gamdf = NULL, imputeX=F){
  
    require(SuperLearner)
    if(is.null(id)) 
        id <- 1:length(Y)
    if(is.null(W)){
        nullW <- TRUE
        fulld <- data.frame(id, Y, Age)
    }else{
        nullW <- FALSE
        Wdesign <- design_matrix(W)
        fulld <- data.frame(id, Y, Age, Wdesign)
    }
    if(imputeX==T){
    As <- seq(0, max(fulld$Age), by=0.1)
    }else{
    As <- unique(fulld$Age)
    }
    pY <- rep(NA, length(As))
    fitd <- fulld[complete.cases(fulld), ]
    n.orig <- dim(fulld)[1]
    n.fit <- dim(fitd)[1]
    

    
    if (n.orig > n.fit) 
        warning(paste("\n\n", n.orig - n.fit, "observations were dropped due to missing values\n in the outcome, age, or adjustement covariates. \n The original dataset contained", 
            n.orig, "observations,\n but GAM_simulCI is fitting the curve using", 
            n.fit, "observations."))
    X <- subset(fitd, select = -c(1:2))
    if (length(grep("SL.gam", SL.library)) > 0) {
      set.seed(123456)
        cvGAM <- ab_cvGAM(Y = fitd$Y, X = X, id = fitd$id, SL.library = SL.library, 
            cvControl = cvControl, df = gamdf)
        SL.library <- cvGAM$SL.library
    }
    
    print(cvGAM$df_opt)

    try(detach(package:gam))
    require(mgcv)
        m <- gam(Y ~ s(Age, k = cvGAM$df_opt), data = fitd,  method = "REML")
        pval<- unlist(summary(m))$s.pv

      
Vb <- vcov(m)
#https://stats.stackexchange.com/questions/110091/how-to-calculate-the-robust-standard-error-of-predicted-y-from-a-linear-regressi
#Vb <- sandwichSE(fitd, m, fitd$id)
newd <- seq(min(Age), max(Age), length = nrow(fitd))
pred <- predict(m, data.frame(Age = newd),  se.fit = TRUE)
se.fit <- pred$se.fit

sandwichSE(fitd, m, fitd$id)

set.seed(123456)
N <- 10000

BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
Cg <- predict(m, data.frame(Age = newd), type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)


absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))


masd <- apply(absDev, 2L, max)
crit <- quantile(masd, prob = 0.95, type = 8)
pred <- transform(cbind(data.frame(pred), newd),
                  uprP = fit + (2 * se.fit),
                  lwrP = fit - (2 * se.fit),
                  uprS = fit + (crit * se.fit),
                  lwrS = fit - (crit * se.fit))  
        
pred<-data.frame(Y=fitd$Y, X=fitd$Age, pred, Pval=rep(pval, nrow(fitd)), degrees.freedom=cvGAM$df_opt)  


    return(pred)    
}



######################################
#GAM temperature and rainfall fits
######################################

# setwd("C:/Users/andre/Dropbox/Trichy analysis/Results/")
# load("temp.datasets.Rdata")

try(setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data"))

load("analysis_datasets.Rdata")
head(d)



set.seed(12345)
fit.temp7.unadj <- GAM_simulCI(Y=d$Y,Age=d$temp.ave7.lag7,id=d$id,  gamdf = c(1:10))
fit.temp14.unadj <- GAM_simulCI(Y=d$Y,Age=d$temp.ave7.lag14,id=d$id,  gamdf = c(1:10))
fit.temp21.unadj <- GAM_simulCI(Y=d$Y,Age=d$temp.ave7.lag21,id=d$id,  gamdf = c(1:10))

save(fit.temp7.unadj,fit.temp14.unadj,fit.temp21.unadj, file="temp.unadjusted.GAMfits.Rdata")



# 
# load("C:/Users/andre/Dropbox/Trichy analysis/Results/rain.datasets.Rdata")
# head(d)
set.seed(12345)
 fit.rain7.unadj <- GAM_simulCI(Y=d$Y,Age=log(d$rain.ave7.lag8+.1),id=d$id, imputeX=T, gamdf = c(1:10))
 fit.rain14.unadj <- GAM_simulCI(Y=d$Y,Age=log(d$rain.ave7.lag15+.1),id=d$id, imputeX=T, gamdf = c(1:10))
 fit.rain21.unadj <- GAM_simulCI(Y=d$Y,Age=log(d$rain.ave7.lag22+.1),id=d$id, imputeX=T, gamdf = c(1:10))
save(fit.rain7.unadj,fit.rain14.unadj,fit.rain21.unadj, file="rain.unadjusted.GAMfits.Rdata")


fit.temp7.unadj$degrees.freedom[1]
fit.temp14.unadj$degrees.freedom[1]
fit.temp21.unadj$degrees.freedom[1]
fit.rain7.unadj$degrees.freedom[1]
fit.rain14.unadj$degrees.freedom[1]
fit.rain21.unadj$degrees.freedom[1]




######################################
#GAM long-term rainfall subgroup fits
######################################

head(d)
df<-subset(d, select=c(Y, id, rain.ave7.lag7, rain.ave7.lag14, rain.ave7.lag21,LT8, LT15, LT22))

  LT<-df$LT8
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  df$LTQ8<-2
  df$LTQ8[LT < ntiles[1]]<-1
  df$LTQ8[LT >= ntiles[2]]<-3
  
  LT<-df$LT15
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  df$LTQ15<-2
  df$LTQ15[LT < ntiles[1]]<-1
  df$LTQ15[LT >= ntiles[2]]<-3
  
  LT<-df$LT22
  ntiles<-as.vector(quantile(LT, probs = seq(0, 1, 0.3333), na.rm=T)[2:3])
  df$LTQ22<-2
  df$LTQ22[LT < ntiles[1]]<-1
  df$LTQ22[LT >= ntiles[2]]<-3
  
   
df <- subset(df, select=c(Y, id, rain.ave7.lag7, rain.ave7.lag14, rain.ave7.lag21,LTQ8, LTQ15, LTQ22))

#Diarrhea outcome fits
set.seed(1111)
fit.rain7.unadjLT1 <- GAM_simulCI(Y=df$Y[df$LTQ8==1],Age=log(df$rain.ave7.lag7[df$LTQ8==1]+.1),id=df$id[df$LTQ8==1], imputeX=T, gamdf = c(1:10))
fit.rain7.unadjLT2 <- GAM_simulCI(Y=df$Y[df$LTQ8==2],Age=log(df$rain.ave7.lag7[df$LTQ8==2]+.1),id=df$id[df$LTQ8==2], imputeX=T, gamdf = c(1:10))
fit.rain7.unadjLT3 <- GAM_simulCI(Y=df$Y[df$LTQ8==3],Age=log(df$rain.ave7.lag7[df$LTQ8==3]+.1),id=df$id[df$LTQ8==3], imputeX=T, gamdf = c(1:10))

fit.rain14.unadjLT1 <- GAM_simulCI(Y=df$Y[df$LTQ15==1],Age=log(df$rain.ave7.lag14[df$LTQ15==1]+.1),id=df$id[df$LTQ15==1], imputeX=T, gamdf = c(1:10))
fit.rain14.unadjLT2 <- GAM_simulCI(Y=df$Y[df$LTQ15==2],Age=log(df$rain.ave7.lag14[df$LTQ15==2]+.1),id=df$id[df$LTQ15==2], imputeX=T, gamdf = c(1:10))
fit.rain14.unadjLT3 <- GAM_simulCI(Y=df$Y[df$LTQ15==3],Age=log(df$rain.ave7.lag14[df$LTQ15==3]+.1),id=df$id[df$LTQ15==3], imputeX=T, gamdf = c(1:10))

fit.rain21.unadjLT1 <- GAM_simulCI(Y=df$Y[df$LTQ22==1],Age=log(df$rain.ave7.lag21[df$LTQ22==1]+.1),id=df$id[df$LTQ22==1], imputeX=T, gamdf = c(1:10))
fit.rain21.unadjLT2 <- GAM_simulCI(Y=df$Y[df$LTQ22==2],Age=log(df$rain.ave7.lag21[df$LTQ22==2]+.1),id=df$id[df$LTQ22==2], imputeX=T, gamdf = c(1:10))
fit.rain21.unadjLT3 <- GAM_simulCI(Y=df$Y[df$LTQ22==3],Age=log(df$rain.ave7.lag21[df$LTQ22==3]+.1),id=df$id[df$LTQ22==3], imputeX=T, gamdf = c(1:10))

save(fit.rain7.unadjLT1,fit.rain7.unadjLT2,fit.rain7.unadjLT3,
     fit.rain14.unadjLT1,fit.rain14.unadjLT2,fit.rain14.unadjLT3,
     fit.rain21.unadjLT1,fit.rain21.unadjLT2,fit.rain21.unadjLT3,
     file="rain.unadjusted.GAMfits.stratified.Rdata")


  



