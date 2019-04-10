

####################
# GAMM AR1 functions
####################



#Make GAMM analysis function
gammAR1 <- function(data, Avar, Yvar="Diarrhea", strat=NULL){
  set.seed(12345)
  df <-data %>% subset(., select=c("Y","stdywk", "individ", "vilid",Avar))
  colnames(df)[colnames(df)==Avar] <- "A"
  df <-df[!is.na(df$A),]
  nlevels <- length(unique(df$A))
  
  m<-NULL
  try(m <- gamm(Y ~ A , data=df, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log')))
  if(is.null(m)){m <- gamm(Y ~ A , data=df, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))}
 
  res<-summary(m$gam)
  #print(res)
  resdf <- data.frame(outcome=rep(Yvar,nlevels-1), exposure=rep(Avar,nlevels-1), 
                      PR=exp(res$p.coeff)[2:nlevels], 
                      ci.lb=exp(res$p.coeff[2:nlevels] - res$se[2:nlevels] * 1.96), 
                      ci.ub=exp(res$p.coeff[2:nlevels] + res$se[2:nlevels] * 1.96))
  if(is.null(strat)){
    resdf$strat="Unstratified"
  }else{
    resdf$strat=strat
  }
  return(resdf)
}


#Make GAMM analysis function
gammAR1_adj <- function(data, Avar, Yvar="Diarrhea", strat=NULL, weathervar, Wvars){
  set.seed(12345)
  df <-data %>% subset(., select=c("Y","stdywk", "individ", "vilid",Avar, weathervar, Wvars))
  colnames(df)[colnames(df)==Avar] <- "A"
  colnames(df)[colnames(df)==weathervar] <- "Weather"

  df <-df[!is.na(df$A),]
  nlevels <- length(unique(df$A))
  
  Wall<-subset(df, select=c("Weather", Wvars))
  Wall<-droplevels(Wall)

  Wscreen <- hbgdki_prescreen(Y=df$Y,Ws=Wall,family="binomial", pval = 0.2 ,print=T, ncases=sum(df$Y, na.rm=T))
  

  d.W<-subset(Wall, select=Wscreen)

 Ws <- colnames(d.W)
 dm <- data.frame(Y=df$Y, A=df$A, stdywk=df$stdywk, vilid=df$vilid, individ=df$individ, d.W)

  
  frm <- paste0("Y ~ A + ", paste0(colnames(dm)[-c(1:5)], collapse = " + "))
  
  m<-NULL
  try(m <- gamm(formula = as.formula(frm), data=dm, random =  list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log')))
  if(is.null(m)){try(m <- gamm(formula = as.formula(frm), data=dm, random =  list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log')))}
  if(is.null(m)){
    #If error from "Singularity in backsolve", try dropping zero-variance covariates:
    d.W<-design_matrix(d.W)
    #remove near zero variance columns
    preproc = caret::preProcess(d.W, method = c("zv", "nzv"))
    d.W = predict(preproc, d.W)
    d.W<-droplevels(d.W)

    Ws <- colnames(d.W)
    dm <- data.frame(Y=df$Y, A=df$A, stdywk=df$stdywk, vilid=df$vilid, individ=df$individ, d.W)
    frm <- paste0("Y ~ A + ", paste0(colnames(dm)[-c(1:5)], collapse = " + "))
    
    m <- gamm(formula = as.formula(frm), data=dm, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
    }
  
  
  res <- summary(m$gam)
  
  resdf <- data.frame(outcome=rep(Yvar,nlevels-1), exposure=rep(Avar,nlevels-1), 
                      PR=exp(res$p.coeff)[2:nlevels], 
                      ci.lb=exp(res$p.coeff[2:nlevels] - res$se[2:nlevels] * 1.96), 
                      ci.ub=exp(res$p.coeff[2:nlevels] + res$se[2:nlevels] * 1.96))
  print(resdf)
  if(is.null(strat)){
    resdf$strat="Unstratified"
  }else{
    resdf$strat=strat
  }
  return(list(resdf=resdf, cov=Wscreen, Avar=Avar, Yvar=Yvar, strat=strat))
}



trichy_gamm <- function(d, A, Y="Diarrhea", strat=NULL, Wvars=NULL, weathervar=NULL){
  set.seed(12345)
  if(is.null(strat)){
    if(is.null(Wvars)){
      res<-gammAR1(data=d, Avar=A, Yvar=Y, strat=NULL)
    }else{
      res<-gammAR1_adj(data=d, Avar=A, Yvar=Y, strat=NULL, weathervar=weathervar, Wvars=Wvars)
    }
  }else{
    V= d[,strat]
    d<-d[!is.na(V),]
    V<-factor(V[!is.na(V)])
    res<-NULL
    for(i in 1:length(unique(V))){
      dsub= d[V==levels(V)[i],]
          if(is.null(Wvars)){
            res_strat <- gammAR1(data=dsub, Avar=A, Yvar=Y, strat=levels(V)[i])
            res<-rbind(res, res_strat)
          }else{
            res_strat<-gammAR1_adj(data=dsub, Avar=A, Yvar=Y, strat=levels(V)[i], weathervar=weathervar, Wvars=Wvars)
            res<-c(res, res_strat)
            names(res)[i] <- levels(V)[i]
          }
    }
    
  }
  if(is.null(Wvars)){
    res$adjusted="N" 
  }else{
    res$resdf$adjusted="Y" 
  }
  return(res)
}






#Complete case sensitivity analysis function

trichy_gammCC <- function(d, A, Y="Diarrhea", strat=NULL, adj_set, Wvars, weathervar=NULL){
  set.seed(12345)
  if(is.null(strat)){
      res<-gammAR1_adjCC(data=d, Avar=A, Yvar=Y, strat=NULL, adj_set=adj_set, weathervar=weathervar, Wvars=Wvars)
  }else{
    V= d[,strat]
    d<-d[!is.na(V),]
    V<-factor(V[!is.na(V)])
    res<-NULL
    for(i in 1:length(unique(V))){
      dsub= d[V==levels(V)[i],]
        res_strat<-gammAR1_adjCC(data=dsub, Avar=A, Yvar=Y, strat=levels(V)[i], adj_set=adj_set[[i]], weathervar=weathervar, Wvars=Wvars)
        res<-c(res, res_strat)
        names(res)[i] <- levels(V)[i]
    }
    
  }
  return(res)
}



#Make GAMM analysis function for complete case sensitivity analysis function
gammAR1_adjCC <- function(data, Avar, Yvar="Diarrhea", strat=NULL, adj_set, weathervar, Wvars){
  set.seed(12345)
  df <-data %>% subset(., select=c("Y","stdywk", "individ", "vilid",Avar, weathervar, Wvars))
  colnames(df)[colnames(df)==Avar] <- "A"
  colnames(df)[colnames(df)==weathervar] <- "Weather"
  
  df <-df[!is.na(df$A),]
  nlevels <- length(unique(df$A))
  
  Wall<-subset(df, select=c("Weather", Wvars))
  Wall<-droplevels(Wall)

  
  d.W<-subset(Wall, select=adj_set)
  
  Ws <- colnames(d.W)
  dm <- data.frame(Y=df$Y, A=df$A, stdywk=df$stdywk, vilid=df$vilid, individ=df$individ, d.W)
  Nfull <- nrow(dm)
  dm <- dm[complete.cases(dm),]
  N <- nrow(dm)
  cat(Nfull-N, " observations dropped\n")
  
  
  frm <- paste0("Y ~ A + ", paste0(colnames(dm)[-c(1:5)], collapse = " + "))
  
  m<-NULL
  try(m <- gamm(formula = as.formula(frm), data=dm, random =  list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=binomial(link='log')))
  if(is.null(m)){try(m <- gamm(formula = as.formula(frm), data=dm, random =  list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log')))}
  if(is.null(m)){
    #If error from "Singularity in backsolve", try dropping zero-variance covariates:
    d.W<-design_matrix(d.W)
    #remove near zero variance columns
    preproc = caret::preProcess(d.W, method = c("zv", "nzv"))
    d.W = predict(preproc, d.W)
    d.W<-droplevels(d.W)
    
    Ws <- colnames(d.W)
    dm <- data.frame(Y=df$Y, A=df$A, stdywk=df$stdywk, vilid=df$vilid, individ=df$individ, d.W)
    frm <- paste0("Y ~ A + ", paste0(colnames(dm)[-c(1:5)], collapse = " + "))
    
    m <- gamm(formula = as.formula(frm), data=dm, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  }
  
  
  res <- summary(m$gam)
  
  resdf <- data.frame(outcome=rep(Yvar,nlevels-1), exposure=rep(Avar,nlevels-1), 
                      PR=exp(res$p.coeff)[2:nlevels], 
                      ci.lb=exp(res$p.coeff[2:nlevels] - res$se[2:nlevels] * 1.96), 
                      ci.ub=exp(res$p.coeff[2:nlevels] + res$se[2:nlevels] * 1.96))
  print(resdf)
  if(is.null(strat)){
    resdf$strat="Unstratified"
  }else{
    resdf$strat=strat
  }
  return(list(resdf=resdf, cov=adj_set, Avar=Avar, Yvar=Yvar, strat=strat, N=N, Ndropped=Nfull-N))
}




####################
#Plotting functions
####################





#------------------------------
# Function to clean results 
# dataframes into plot formats
#------------------------------


plotdf_format <- function(
                  temp = temp.unadj,
                  Hrain = Hrain.unadj,
                  Hrain1 = Hrain.unadjLT1,
                  Hrain2 = Hrain.unadjLT2,
                  Hrain3 = Hrain.unadjLT3,
                  H2S = H2S.Hrain.unadj,
                  H2S1 = H2s.unadjLT1,
                  H2S2 = H2s.unadjLT2,
                  H2S3 = H2s.unadjLT3){

    raindf <- data.frame(exposure=rep("Heavy rain",12),
                 strata=c(rep("Unstratified",3),
                 rep("Low",3),
                 rep("Medium",3),
                 rep("High",3)),
                 x=as.character(c(1,1,1,2,2,2,3,3,3,4,4,4)),
                 lag=rep(c("One week lag", "Two week lag", "Three week lag"),4),
                  rbind(
                  Hrain,
                  Hrain1,
                  Hrain2,
                  Hrain3))

    raindf$strata <- factor(raindf$strata, levels=unique(raindf$strata))
    raindf$lag <- factor(raindf$lag, levels=unique(raindf$lag))
    
    raindf<- subset(raindf, select=-c(a,b,c,d, logPR,se.logPR,Z,p))
    
    #Add blank rows for reference levels
    temp <- rbind(NA,temp[1:3,], NA,temp[4:6,], NA,temp[7:9,])
    rownames(temp) <- NULL
    tempdf <- data.frame(exposure=rep("Temperature",12),
                 strata=rep(c("(ref.)","1v2","1v3","1v4"),3),
                              x=as.character(c(1,2,3,4,1,2,3,4,1,2,3,4)),
                 lag=c(rep("One week lag",4), rep("Two week lag",4), rep("Three week lag",4)),
                             temp) %>% 
      subset(., select=-c(logPR,se.logPR,Z,p))
    tempdf$strata <- factor(tempdf$strata, levels=unique(tempdf$strata))
    tempdf$lag <- factor(tempdf$lag, levels=unique(tempdf$lag))
    
    
    H2Sdf <- data.frame(exposure=rep("Heavy rain",12),
                 strata=c(rep("Unstratified",3),
                 rep("Low",3),
                 rep("Medium",3),
                 rep("High",3)),
                 x=as.character(c(1,1,1,2,2,2,3,3,3,4,4,4)),
                 lag=rep(c("One week lag", "Two week lag", "Three week lag"),4),
                  rbind(
                  H2S,
                  H2S1,
                  H2S2,
                  H2S3))

    H2Sdf$strata <- factor(H2Sdf$strata, levels=unique(H2Sdf$strata))
    H2Sdf$lag <- factor(H2Sdf$lag, levels=unique(H2Sdf$lag))
    
    H2Sdf<- subset(H2Sdf, select=-c(a,b,c,d, logPR,se.logPR,Z,p))

    df <- list(raindf=raindf, H2Sdf=H2Sdf, tempdf=tempdf)
    
  return(df)
}
  









# --------------------------------------
# function to prescreen covariates and 
# select 1 covariate per 10 cases
# --------------------------------------

hbgdki_prescreen <- function (Y, Ws, ncases, family = "binomial", pval = 0.2,  print = TRUE){
  
   n<-nrow(Ws)
   if(n-ncases < ncases){ncases<-n-ncases}  
  
    require(lmtest)
    if(family[[1]]=="neg.binom"){
       require(MASS)
    }
    if(pval > 0.99 | pval < 0){
        stop("P-value threshold not set between 0 and 1.")
    }
    Ws <- as.data.frame(Ws)
    dat <- data.frame(Ws, Y)
    dat <- dat[complete.cases(dat), ]
    nW <- ncol(Ws)
    LRp <- matrix(rep(NA, nW), nrow = nW, ncol = 1)
    rownames(LRp) <- names(Ws)
    colnames(LRp) <- "P-value"
    if(family[[1]] != "neg.binom"){
        for(i in 1:nW) {
            dat$W <- dat[, i]
            if(class(dat$W) == "factor" & dim(table(dat$W)) == 
                1) {
                fit1 <- fit0 <- glm(Y ~ 1, data = dat, family = family)
            }
            else{
                fit1 <- glm(Y ~ W, data = dat, family = family)
                fit0 <- glm(Y ~ 1, data = dat, family = family)
            }
            LRp[i] <- lrtest(fit1, fit0)[2, 5]
        }
    }
    else{
        if(!requireNamespace("MASS", quietly = TRUE)){
            stop("Pkg needed forthis function to work. Please install it.", 
                call. = FALSE)
        }
        else{
            for(i in 1:nW){
                dat$W <- dat[, i]
                if(class(dat$W) == "factor" & dim(table(dat$W)) == 
                  1) {
                  fit1 <- fit0 <- glm(Y ~ 1, data = dat, family = family)
                }
                else{
                  fit1 <- glm.nb(Y ~ W, data = dat, family = family)
                  fit0 <- glm.nb(Y ~ 1, data = dat, family = family)
                }
                LRp[i] <- lrtest(fit1, fit0)[2, 5]
            }
        }
    }
    p20 <- ifelse(LRp < pval, 1, 0)
    if(print == TRUE) {
        cat("\nLikelihood Ratio Test P-values:\n")
        print(round(LRp, 5))
        if(sum(p20) > 0) {
            LRps <- matrix(LRp[p20 == 1, ], ncol = 1)
            rownames(LRps) <- names(Ws)[p20 == 1]
            colnames(LRps) <- "P-value"
            cat(paste("\n\nCovariates selected (P<", pval, "):\n", 
                sep = ""))
            print(LRps)
        }
        else{
            cat(paste("\nNo covariates were associated with the outcome at P<", 
                pval))
        }
    }
    
    W <- data.frame(wvar=names(Ws), p=as.numeric(LRp), pthres=as.numeric(p20))
    if(floor(ncases/10) > 0){
    W <- W %>% arrange(p) %>% slice(1:floor(ncases/10))
                cat(paste("\n\nTop ",floor(ncases/10), " Covariates selected:\n", 
                sep = ""))
     print(W$wvar)

    }else{
      W$pthres = 0
      cat("\nNot enough cases for adjusted analysis")
    }
    return(as.character(W$wvar[W$pthres == 1]))
}

