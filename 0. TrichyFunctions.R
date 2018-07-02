




####################
#Quartile analysis function
####################


quart_tmle<-function(d,Wvars=NULL,weather="tempQ7",q=c(1,4),cutoff=4, library=NULL){
colnames(d)[which( colnames(d)==paste0(weather))]="quartile"

# df<-d%>%
#   mutate(quartile = ntile(weather, cutoff))%>%
#   subset(quartile==q[1]|quartile==q[2])%>%
#   mutate(A=ifelse(quartile==q[2],1,0))

#d$quartile<-cut(d$weather, breaks = c(0,tempQ,999), labels=c("Q1","Q2","Q3","Q4"))

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
#Heavy rain TMLE function
####################

HeavyRainTMLE<-function(dat=d,i,Wvars=NULL,sl.library=library){
  W<-NULL
  res<-matrix(NA, nrow=1, ncol=8)
  set.seed(12345)
 colnames(dat)[which(colnames(dat)==paste0("HeavyRain.lag",i))]="HeavyRain"

 df<-dat

  a=round(sum(df$HeavyRain==1 & df$Y==1, na.rm=T),0)
  b=round(sum(df$HeavyRain==1 & df$Y==0, na.rm=T),0)
  c=round(sum(df$HeavyRain==0 & df$Y==1, na.rm=T),0)
  D=round(sum(df$HeavyRain==0 & df$Y==0, na.rm=T),0)
  
if(!is.null(Wvars)){
  Wall<-subset(df, select=Wvars)
}
if(is.null(Wvars)){
  
      res<-as.numeric(washb_glm(Y=df$Y, tr=df$HeavyRain, id=df$id, contrast=c(0,1), family=poisson(link='log'))$TR)
    res.vec<-c(a, b, c, D, res)
    names(res.vec) <- c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR", "Z", "p")

}else{
  Keep<-rep(T, ncol(Wall))
  for(i in 1:ncol(Wall)){
    if(length(levels(Wall[,i]))==1){
      print(table(Wall$wqsource))
      Keep[i]<-F
    }
  }
  
  Wall<-droplevels(Wall)
  Wscreen <- washb_prescreen(Y=df$Y,Ws=Wall,family="binomial", pval = 0.2 ,print=T)
  d.W<-subset(Wall, select=Wscreen)
  d.W<-design_matrix(d.W)
  #remove near zero variance columns
  preproc = caret::preProcess(d.W, method = c("zv", "nzv"))
  d.W = predict(preproc, d.W)
  fit<-tmle(Y=df$Y,A=df$HeavyRain,W=d.W, family="binomial",Q.SL.library=sl.library,g.SL.library=sl.library,id=df$id)
  res[1,c(1:4)]<-c(a,b,c,D)
  res[1,5]<-fit$estimates$RR$psi
  res[1,c(6,7)]<-fit$estimates$RR$CI
  res[1,8]<-fit$estimates$RR$pvalue
  res.vec<-as.vector(res)
}
  return(res.vec)
}
  



####################
# H2S Heavy rain TMLE function
####################


H2STMLE<-function(dat=d,i,Wvars=NULL,sl.library=library, seed=12345){
  W<-NULL
  res<-matrix(NA, nrow=1, ncol=8)
  set.seed(seed)
 colnames(dat)[which(colnames(dat)==paste0("HeavyRain.lag",i))]="HeavyRain"
 df<-dat
 df$h2s <- df$H2S

  a=round(sum(df$HeavyRain==1 & df$h2s==1, na.rm=T),0)
  b=round(sum(df$HeavyRain==1 & df$h2s==0, na.rm=T),0)
  c=round(sum(df$HeavyRain==0 & df$h2s==1, na.rm=T),0)
  D=round(sum(df$HeavyRain==0 & df$h2s==0, na.rm=T),0)
  cat(a,"\n",b,"\n",c,"\n",D,"\n")
  
if(!is.null(Wvars)){
  Wall<-subset(df, select=Wvars)
}
if(is.null(Wvars)){
 
    res<-as.numeric(washb_glm(Y=df$h2s, tr=df$HeavyRain, id=df$id, contrast=c(0,1), family=poisson(link='log'))$TR)
    res.vec<-c(a, b, c, D, res)
    names(res.vec) <- c("a","b","c","d","PR", "ci.lb", "ci.ub", "logPR", "se.logPR", "Z", "p")

}else{
  Wall<-Wall[!is.na(df$h2s),]
  Wall<-droplevels(Wall)

  Wscreen <- washb_prescreen(Y=df$h2s[!is.na(df$h2s)],Ws=Wall,family="binomial", pval = 0.2 ,print=T)
  d.W<-subset(Wall, select=Wscreen)
  d.W<-design_matrix(d.W)
  #remove near zero variance columns
  preproc = caret::preProcess(d.W, method = c("zv", "nzv"))
  d.W = predict(preproc, d.W)
  d.W<-droplevels(d.W)
  fit<-tmle(Y=df$h2s[!is.na(df$h2s)],A=df$HeavyRain[!is.na(df$h2s)],W=d.W, family="binomial",Q.SL.library=sl.library,g.SL.library=sl.library,id=df$id[!is.na(df$h2s)])
  res[1,c(1:4)]<-c(a,b,c,D)
  res[1,5]<-fit$estimates$RR$psi
  res[1,c(6,7)]<-fit$estimates$RR$CI
  res[1,8]<-fit$estimates$RR$pvalue
  res.vec<-as.vector(res)
}
  return(res.vec)
}



####################
# GAMM AR1 functions
####################



#Make GAMM analysis function
gammAR1 <- function(data, Avar, Yvar="Diarrhea", strat=NULL){
  df <-data %>% subset(., select=c("Y","stdywk", "individ", "vilid",Avar))
  colnames(df)[colnames(df)==Avar] <- "A"
  df <-df[!is.na(df$A),]
  nlevels <- length(unique(df$A))
  
  m <- gamm(Y ~ A , data=df, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
 
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
 
  df <-data %>% subset(., select=c("Y","stdywk", "individ", "vilid",Avar, weathervar, Wvars))
  colnames(df)[colnames(df)==Avar] <- "A"
  colnames(df)[colnames(df)==weathervar] <- "Weather"

  df <-df[!is.na(df$A),]
  nlevels <- length(unique(df$A))
  
  #Wall<-subset(df, select=c("Weather", Wvars[-c(1:34,75:86)]))
  Wall<-subset(df, select=c("Weather", Wvars))
  Wall<-droplevels(Wall)

  Wscreen <- washb_prescreen(Y=df$Y,Ws=Wall,family="binomial", pval = 0.2 ,print=T)
  

  d.W<-subset(Wall, select=Wscreen)
  # d.W<-design_matrix(d.W)
  # #remove near zero variance columns
  # preproc = caret::preProcess(d.W, method = c("zv", "nzv"))
  # d.W = predict(preproc, d.W)
  # d.W<-droplevels(d.W)
  
 Ws <- colnames(d.W)
 dm <- data.frame(Y=df$Y, A=df$A, stdywk=df$stdywk, vilid=df$vilid, individ=df$individ, d.W)

  
  frm <- paste0("Y ~ A + ", paste0(colnames(dm)[-c(1:5)], collapse = " + "))
  m <- gamm(formula = as.formula(frm), data=dm, random = list(vilid=~1), correlation = corAR1(form = ~ stdywk|individ), family=poisson(link='log'))
  res <- summary(m$gam)
  
  
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

trichy_gamm <- function(d, A, Y="Diarrhea", strat=NULL, Wvars=NULL, weathervar=NULL){
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
          }else{
            res_strat<-gammAR1_adj(data=dsub, Avar=A, Yvar=Y, strat=levels(V)[i], weathervar=weathervar, Wvars=Wvars)
          }
      
      res<-rbind(res, res_strat)
    }
    
  }
  if(is.null(Wvars)){
    res$adjusted="N" 
  }else{
    res$adjusted="Y" 
  }
  return(res)
}





screen.lasso <- function(Y, X, family, alpha = 1, minscreen = 2, nfolds = 5, nlambda = 20){
    require("glmnet")

      if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
    }
  
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = "deviance", 
        nfolds = nfolds, family = "binomial", alpha = alpha, 
        nlambda = nlambda)
    whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 
        0)
    if (sum(whichVariable) < minscreen) {
        warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
        sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, 
            function(x) sum((x != 0)))
        newCut <- which.max(sumCoef >= minscreen)
        whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, 
            newCut] != 0)
    }
    return(X[,whichVariable])
}




####################
#Plotting functions
####################

#axis formatting
scaleFUN <- function(x) sprintf("%.2f", x)



# 
# RR_plotfun <- function(df, xlab="Quartile contrast", title="", yticks=c(0.125,0.25,0.5,1,2,4,8), cols= tempcol){
#   
#   p<-ggplot(df, aes(x=strata)) + 
#     geom_point(aes(y=PR, fill=strata, color=strata), size = 4) +
#     geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strata),
#                    alpha=0.5, size = 3) +
#     labs(x = xlab, y = "Prevalence Ratio") +
#     geom_hline(yintercept = 1) +
#     coord_cartesian(ylim=range(yticks)) +
#     scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
#     scale_fill_manual(values=rep(cols,4)) +
#     scale_colour_manual(values=rep(cols,4)) +
#     theme(strip.background = element_blank(),
#       legend.position="none",
#       strip.text.x = element_text(size=12),
#       axis.text.x = element_text(size=12)) +
#     facet_wrap(~lag,  scales = "fixed") +
#     ggtitle(title)
#   
#   return(p)
# }
# 
# 
# 
# 
# 
# 
# 
# #------------------------------
# # Function to combine plots 
# # into one figure
# #------------------------------
# 
# trichy_panel_plot <- function(raind=raindf, H2Sd=H2Sdf, tempd=tempdf){
#   
# prain <- RR_plotfun(df=raind, xlab="Long-term rainfall strata", title="Heavy rainfall - diarrhea association",cols= raincol)
# pH2S <- RR_plotfun(df=H2Sd, xlab="Long-term rainfall strata", title="Heavy rainfall - drinking water H2S association", cols= raincol)
# ptemp <- RR_plotfun(df=tempd, xlab="Quartile contrast", title="Temperature - diarrhea association", cols= tempcol)
#   
#   p <- plot_grid(ptemp, prain, pH2S, labels = c("A", "B","C"), ncol = 1)
#   p
#   
#  return(p) 
# }



#------------------------------
# Function to combine plots 
# into one figure
#------------------------------
PR_plotfun <- function(df, xlab="",title="", yticks=c(0.125,0.25,0.5,1,2,4,8)){

  df$strat <- factor(df$strat, levels=unique(df$strat))
  df$lag <- factor(df$lag, levels=unique(df$lag))

  p<-ggplot(df, aes(x=strat)) + 
      geom_point(aes(y=PR, fill=strat, color=strat), size = 4) +
      geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strat),
                     alpha=0.5, size = 3) +
      labs(x = xlab, y = "Prevalence Ratio") +
      geom_hline(yintercept = 1) +
      coord_cartesian(ylim=range(yticks)) +
      scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
      scale_fill_manual(name = "strat", values=cols) +
      scale_colour_manual(name = "strat",values=cols) +
      theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1)) +
      facet_wrap(~lag,  scales = "fixed") +
      ggtitle(title)
  return(p)
}


trichy_panel_plot <- function(d){
  
    #plot parameters
    theme_set(theme_bw())
    cols=c("(ref.)"="#000000",
           "Unstratified"="#000000",
           "1v2"="#ED8021", 
           "1v3"="#E15324", 
           "1v4"="#D62728", 
           "Low"="#56B4E9", 
           "Medium"="#4C8FE5", 
           "High"="#426AE3")
           
    
    #Clean data frames 
     d$Xvar <- d$exposure
     d$lag <- str_split_fixed(d$exposure, "Q|lag",2)[,2]
     d$exposure <- str_split_fixed(d$exposure, "Q|lag",2)[,1]
    
    d$quartile<-str_sub( str_split_fixed(rownames(d), "Q",2)[,2],1,1)
    d$strat[d$strat=="Unstratified" & d$quartile!=""] <- d$quartile[d$strat=="Unstratified" & d$quartile!=""]
    
    #Add blank rows for reference levels
    if(any(grepl("temp",d$exposure))){
      temp_blank <- d[1:3,]
      temp_blank$quartile <- 1
      temp_blank$strat <- "(ref.)"
      temp_blank$lag <- c("7","14","21") 
      temp_blank$PR <- 1
      temp_blank$ci.lb <- temp_blank$ci.ub <-temp_blank$Xvar <- NA
      d <- rbind(temp_blank,d)
    }
    
      #Clean up labels
    d$lag[d$lag=="7"] <- "One week lag"
    d$lag[d$lag=="14"] <- "Two week lag"
    d$lag[d$lag=="21"] <- "Three week lag"
    
    d$strat[d$strat=="2"] <- "1v2"
    d$strat[d$strat=="3"] <- "1v3"
    d$strat[d$strat=="4"] <- "1v4"
    
    d$strat[d$strat=="T1"] <- "Low"
    d$strat[d$strat=="T2"] <- "Medium"
    d$strat[d$strat=="T3"] <- "High"
    
    d$strat[d$strat=="temp"] <- "Temperature"
    d$strat[d$strat=="HeavyRain."] <- "Heavy Rain"
    
    
    #Split into plot dataframes
    d1<-d[grepl("temp",d$exposure) & grepl("Diarrhea",d$outcome) ,]
    d2<-d[grepl("HeavyRain",d$exposure) & grepl("Diarrhea",d$outcome) ,]
    d3<-d[grepl("temp",d$exposure) & grepl("H2S",d$outcome) ,]
    d4<-d[grepl("HeavyRain",d$exposure) & grepl("H2S",d$outcome) ,]  
      
    ptemp <- prain <- ptempH2s <- pH2S <- NULL
      if(nrow(d1)>0){
        ptemp <- PR_plotfun(df=d1, xlab="Quartile contrast", title="Temperature - diarrhea association")
      }
    
      if(nrow(d2)>0){
        prain <- PR_plotfun(df=d2, xlab="Long-term rainfall strata", title="Heavy rainfall - diarrhea association")
      }
    
      if(nrow(d3)>0){
        ptempH2s <- PR_plotfun(df=d3, xlab="Quartile contrast", title="Temperature - drinking water H2S association")
      }
    
      if(nrow(d4)>0){
        pH2S <- PR_plotfun(df=d4, xlab="Long-term rainfall strata", title="Heavy rainfall - drinking water H2S association")
      }
    
    plist <- list(ptemp, prain, ptempH2s,pH2S)
    names(plist) <- seq_along(plist)
    plist[sapply(plist, is.null)] <- NULL
    
    nplots <- length(plist)
    if(nplots==1){
      p <- plist[[1]]
    }
    if(nplots==2){
      p <- plot_grid(plist[[1]], plist[[2]],  labels = c("A", "B"), ncol = 1)
    }
    if(nplots==3){
      p <- plot_grid(plist[[1]], plist[[2]], plist[[3]], labels = c("A", "B","C"), ncol = 1)
    }
  
 return(p) 
}



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
  



#Print tables to two decimals places

tab_format <- function(d){
  d$PR <- sprintf("%1.2f", d$PR)
    d$ci.lb <- sprintf("%1.2f", d$ci.lb)
  d$ci.ub <- sprintf("%1.2f", d$ci.ub)
  
  d$tab_format= paste0(d$PR," (", d$ci.lb,", ",d$ci.ub,")")

  return(d)
}

