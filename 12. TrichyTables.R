

#Trichy analysis table printing
rm(list=ls())
library(xtable)



# Table functions
cleantable <- function(x,digits) {
 print( xtable(x,digits=digits),
        sanitize.text.function=function(y){y},
        floating=FALSE,
        include.rownames=FALSE,
        include.colnames=FALSE,
        only.contents=TRUE,
        hline.after=NULL
 )
}




#Results tables
setwd("C:/Users/andre/Dropbox/Trichy analysis/Results/")
load("temp.Results.unadjusted.Rdata")
Mean_lag7
Mean_lag14
Mean_lag21
temp.unadj
load("temp.Results.adjusted.Rdata")



TempTab<-data.frame(lag=rep(NA,12),
                    Q=rep(NA,12),
                    N=rep(NA,12),
                    ndiar=rep(NA,12),
                    meandiar=rep(NA,12),
                    meantemp=rep(NA,12),
                    PR=rep(NA,12),
                    P=rep(NA,12),
                    adjPR=rep(NA,12),
                    adjP=rep(NA,12))
TempTab$lag<-c(rep("One week lag",4),rep("Two week lag",4),rep("Three week lag",4))             
TempTab$Q<-rep(c("~~Quartile 1","~~Quartile 2","~~Quartile 3","~~Quartile 4"),3)
TempTab$N[1:4]<-Mean_lag7$quart.n
TempTab$N[5:8]<-Mean_lag14$quart.n
TempTab$N[9:12]<-Mean_lag21$quart.n
TempTab$ndiar[1:4]<-sprintf("%1.0f", Mean_lag7$quart.n*Mean_lag7$meandiar)
TempTab$ndiar[5:8]<-sprintf("%1.0f", Mean_lag14$quart.n*Mean_lag14$meandiar)
TempTab$ndiar[9:12]<-sprintf("%1.0f", Mean_lag21$quart.n*Mean_lag21$meandiar)
TempTab$meandiar[1:4]<-Mean_lag7$meandiar
TempTab$meandiar[5:8]<-Mean_lag14$meandiar
TempTab$meandiar[9:12]<-Mean_lag21$meandiar
TempTab$meantemp[1:4]<-Mean_lag7$meantemp
TempTab$meantemp[5:8]<-Mean_lag14$meantemp
TempTab$meantemp[9:12]<-Mean_lag21$meantemp
TempTab$PR[c(2:4,6:8,10:12)]<-paste(sprintf("%1.2f",temp.unadj[,1])," (",sprintf("%1.2f",temp.unadj[,2]),", ",sprintf("%1.2f",temp.unadj[,3]),")", sep="")
TempTab$P[c(2:4,6:8,10:12)]<-sprintf("%1.3f",temp.unadj[,7])
TempTab$adjPR[2]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v2$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v2$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v2$estimates$RR$CI)[2]),")", sep="")
TempTab$adjPR[3]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v3$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v3$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v3$estimates$RR$CI)[2]),")", sep="")
TempTab$adjPR[4]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v4$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v4$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag7.1v4$estimates$RR$CI)[2]),")", sep="")
TempTab$adjPR[6]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v2$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v2$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v2$estimates$RR$CI)[2]),")", sep="")
TempTab$adjPR[7]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v3$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v3$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v3$estimates$RR$CI)[2]),")", sep="")
TempTab$adjPR[8]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v4$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v4$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag14.1v4$estimates$RR$CI)[2]),")", sep="")
TempTab$adjPR[10]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v2$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v2$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v2$estimates$RR$CI)[2]),")", sep="")
TempTab$adjPR[11]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v3$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v3$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v3$estimates$RR$CI)[2]),")", sep="")
TempTab$adjPR[12]<-paste(sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v4$estimates$RR$psi))," (",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v4$estimates$RR$CI)[1]),", ",sprintf("%1.2f",unlist(temp.adj$temp.ave7.lag21.1v4$estimates$RR$CI)[2]),")", sep="")
TempTab$adjP[2]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag7.1v2$estimates$RR$pvalue))
TempTab$adjP[3]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag7.1v3$estimates$RR$pvalue))
TempTab$adjP[4]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag7.1v4$estimates$RR$pvalue))
TempTab$adjP[6]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag14.1v2$estimates$RR$pvalue))
TempTab$adjP[7]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag14.1v3$estimates$RR$pvalue))
TempTab$adjP[8]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag14.1v4$estimates$RR$pvalue))
TempTab$adjP[10]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag21.1v2$estimates$RR$pvalue))
TempTab$adjP[11]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag21.1v3$estimates$RR$pvalue))
TempTab$adjP[12]<-sprintf("%1.3f",unlist(temp.adj$temp.ave7.lag21.1v4$estimates$RR$pvalue))


TempTab[c(1,5,9),c(7,9)]<-"Referent"
TempTab[,5]<-sprintf("%1.1f",TempTab[,5]*100)

#Drop p-values 
TempTab<-TempTab[,-c(8,10)]
cleantable(TempTab[,-1],2)



#----------------------------------------
# Heavy rain- unstratified
#----------------------------------------
rm(list=ls())
library(xtable)



# Table functions
cleantable <- function(x,digits) {
 print( xtable(x,digits=digits),
        sanitize.text.function=function(y){y},
        floating=FALSE,
        include.rownames=FALSE,
        include.colnames=FALSE,
        only.contents=TRUE,
        hline.after=NULL
 )
}

load("rain.Results.unadjusted.Rdata")
load("rain.Results.adjusted.Rdata")
load("HRmean.Results.Rdata")
load("rain.Results.adjusted.Rdata")
load("H2S.Results.adjusted.Rdata")



HRTab<-data.frame(lag=rep(NA,6),
                    Q=rep(NA,6),
                    N=rep(NA,6),
                    ndiar=rep(NA,6),
                    meandiar=rep(NA,6),
                    meanHR=rep(NA,6),
                    PR=rep(NA,6),
                    P=rep(NA,6),
                    adjPR=rep(NA,6),
                    adjP=rep(NA,6))
HRTab$lag<-c(rep("One week lag",2),rep("Two week lag",2),rep("Three week lag",2))             
HRTab$Q<-rep(c("~~No heavy rainfall","~~Heavy rainfall"),3)
HRTab$N[c(1,3,5)]<-Hrain.unadj[,3]+Hrain.unadj[,4]
HRTab$N[c(2,4,6)]<-Hrain.unadj[,1]+Hrain.unadj[,2]
HRTab$ndiar[c(1,3,5)]<-Hrain.unadj[,3]
HRTab$ndiar[c(2,4,6)]<-Hrain.unadj[,1]
HRTab$meandiar <- HRTab$ndiar / HRTab$N * 100
HRTab$meanHR[c(1:2)]<-HR.rainmean7[,2]
HRTab$meanHR[c(3:4)]<-HR.rainmean14[,2]
HRTab$meanHR[c(5:6)]<-HR.rainmean21[,2]
HRTab$PR[c(1,3,5)] <- "Referent"
HRTab$PR[c(2,4,6)] <- paste(sprintf("%1.2f",Hrain.unadj[,5])," (",sprintf("%1.2f",Hrain.unadj[,6]),", ",sprintf("%1.2f",Hrain.unadj[,7]),")", sep="")
HRTab$P[c(2,4,6)]<-sprintf("%1.3f",Hrain.unadj[,11])
HRTab$adjPR[c(1,3,5)] <- "Referent"
HRTab$adjPR[c(2,4,6)] <- paste(sprintf("%1.2f",Hrain.adj[,5])," (",sprintf("%1.2f",Hrain.adj[,6]),", ",sprintf("%1.2f",Hrain.adj[,7]),")", sep="")
HRTab$adjP[c(2,4,6)]<-sprintf("%1.3f",Hrain.adj[,8])


 HRTab[,3]<-paste(sprintf("%1.0f", HRTab[,3]))
 HRTab[,4]<-paste(sprintf("%1.0f", HRTab[,4]))
 HRTab[,5]<-paste(sprintf("%1.1f", HRTab[,5]))
 
 #drop p-values 
 HRTab<-HRTab[,-c(8,10)]
cleantable(HRTab[,-1],2)


#----------------------------------------
# Heavy rain- stratified
#----------------------------------------

load("C:/Users/andre/Dropbox/Trichy analysis/Results/HRmean-stratified.Results.Rdata")


HRstrat<-data.frame(lag=rep(NA,18),
                    Q=rep(NA,18),
                    N=rep(NA,18),
                    ndiar=rep(NA,18),
                    meandiar=rep(NA,18),
                    meanLongRain=rep(NA,18),
                    meanHR=rep(NA,18),
                    PR=rep(NA,18),
                    P=rep(NA,18),
                    adjPR=rep(NA,18),
                    adjP=rep(NA,18))
HRstrat$lag<-rep(c(rep("~~One week lag",2),rep("~~Two week lag",2),rep("~~Three week lag",2)),3)             
HRstrat$Q<-rep(c("No heavy rainfall","Heavy rainfall"),9)
HRstrat$N[c(1,3,5)]<-Hrain.unadjLT1[,3]+Hrain.unadjLT1[,4]
HRstrat$N[c(2,4,6)]<-Hrain.unadjLT1[,1]+Hrain.unadjLT1[,2]
HRstrat$N[c(7,9,11)]<-Hrain.unadjLT2[,3]+Hrain.unadjLT2[,4]
HRstrat$N[c(8,10,12)]<-Hrain.unadjLT2[,1]+Hrain.unadjLT2[,2]
HRstrat$N[c(13,15,17)]<-Hrain.unadjLT3[,3]+Hrain.unadjLT3[,4]
HRstrat$N[c(14,16,18)]<-Hrain.unadjLT3[,1]+Hrain.unadjLT3[,2]
HRstrat$ndiar[c(1,3,5)]<-Hrain.unadjLT1[,3]
HRstrat$ndiar[c(2,4,6)]<-Hrain.unadjLT1[,1]
HRstrat$ndiar[c(7,9,11)]<-Hrain.unadjLT2[,3]
HRstrat$ndiar[c(8,10,12)]<-Hrain.unadjLT2[,1]
HRstrat$ndiar[c(13,15,17)]<-Hrain.unadjLT3[,3]
HRstrat$ndiar[c(14,16,18)]<-Hrain.unadjLT3[,1]
HRstrat$meandiar <- HRstrat$ndiar / HRstrat$N * 100
HRstrat$meanLongRain<-
    c(0.6746329,NA,0.7135688,NA,0.7690993,NA,
    2.7613850,NA, 2.7933568,NA, 2.8433577,NA,
    6.5649825,NA, 6.6059397,NA, 6.6059397,NA)
HRstrat$meanHR[c(1:2,7:8,13:14)]<-HR.strat.rainmean7[,3]
HRstrat$meanHR[c(3:4,9:10,15:16)]<-HR.strat.rainmean14[,3]
HRstrat$meanHR[c(5:6,11:12,17:18)]<-HR.strat.rainmean21[,3]
HRstrat$PR[c(1,3,5)] <- "Referent"
HRstrat$PR[c(2,4,6)] <- paste(sprintf("%1.2f",Hrain.unadjLT1[,5])," (",sprintf("%1.2f",Hrain.unadjLT1[,6]),", ",sprintf("%1.2f",Hrain.unadjLT1[,7]),")", sep="")
HRstrat$P[c(2,4,6)]<-sprintf("%1.3f",Hrain.unadjLT1[,11])
HRstrat$adjPR[c(2,4,6)] <- paste(sprintf("%1.2f",Hrain.adjLT1[,5])," (",sprintf("%1.2f",Hrain.adjLT1[,6]),", ",sprintf("%1.2f",Hrain.adjLT1[,7]),")", sep="")
HRstrat$adjP[c(2,4,6)]<-sprintf("%1.3f",Hrain.adjLT1[,8])

HRstrat$PR[c(8,10,12)] <- paste(sprintf("%1.2f",Hrain.unadjLT2[,5])," (",sprintf("%1.2f",Hrain.unadjLT2[,6]),", ",sprintf("%1.2f",Hrain.unadjLT2[,7]),")", sep="")
HRstrat$P[c(8,10,12)]<-sprintf("%1.3f",Hrain.unadjLT2[,11])
HRstrat$adjPR[c(8,10,12)] <- paste(sprintf("%1.2f",Hrain.adjLT2[,5])," (",sprintf("%1.2f",Hrain.adjLT2[,6]),", ",sprintf("%1.2f",Hrain.adjLT2[,7]),")", sep="")
HRstrat$adjP[c(8,10,12)]<-sprintf("%1.3f",Hrain.adjLT2[,8])

HRstrat$PR[c(14,16,18)] <- paste(sprintf("%1.2f",Hrain.unadjLT3[,5])," (",sprintf("%1.2f",Hrain.unadjLT3[,6]),", ",sprintf("%1.2f",Hrain.unadjLT3[,7]),")", sep="")
HRstrat$P[c(14,16,18)]<-sprintf("%1.3f",Hrain.unadjLT3[,11])
HRstrat$adjPR[c(14,16,18)] <- paste(sprintf("%1.2f",Hrain.adjLT3[,5])," (",sprintf("%1.2f",Hrain.adjLT3[,6]),", ",sprintf("%1.2f",Hrain.adjLT3[,7]),")", sep="")
HRstrat$adjP[c(14,16,18)]<-sprintf("%1.3f",Hrain.adjLT3[,8])
HRstrat[c(2,4,6,8,10,12,14,16,18),1]<-""


 HRstrat[,3]<-paste(sprintf("%1.0f", HRstrat[,3]))
 HRstrat[,4]<-paste(sprintf("%1.0f", HRstrat[,4]))
 HRstrat[,5]<-paste(sprintf("%1.1f", HRstrat[,5]))

 #drop p-value 
 HRstrat<-HRstrat[,-c(9,11)]
 
cleantable(HRstrat,2)

      
     
     

#----------------------------------------
# Hydrogen sulfide- unstratified
#----------------------------------------

h2sTab<-data.frame(lag=rep(NA,24),
                    Q=rep(NA,24),
                    N=rep(NA,24),
                    nh2s=rep(NA,24),
                    meanh2s=rep(NA,24),
                    meanLongRain=rep(NA,24),
                    meanHR=rep(NA,24),
                    PR=rep(NA,24),
                    P=rep(NA,24),
                    adjPR=rep(NA,24),
                    adjP=rep(NA,24))

H2S.Hrain.unadj<-rbind(
  H2S.Hrain.unadj,
  H2s.unadjLT1,
  H2s.unadjLT2,
  H2s.unadjLT3)

H2S.Hrain.adj<-rbind(
  H2S.Hrain.adj,
  H2S.adjLT1,
  H2S.adjLT2,
  H2S.adjLT3)

h2sTab$lag<-rep(c(rep("~~One week lag",2),rep("~~Two week lag",2),rep("~~Three week lag",2)),4)             
h2sTab$Q<-rep(c("No heavy rainfall","Heavy rainfall"),12)
h2sTab$N[c(1,3,5,7,9,11,13,15,17,19,21,23)]<-H2S.Hrain.unadj[,3]+H2S.Hrain.unadj[,4]
h2sTab$N[c(2,4,6,8,10,12,14,16,18,20,22,24)]<-H2S.Hrain.unadj[,1]+H2S.Hrain.unadj[,2]
h2sTab$nh2s[c(1,3,5,7,9,11,13,15,17,19,21,23)]<-H2S.Hrain.unadj[,3]
h2sTab$nh2s[c(2,4,6,8,10,12,14,16,18,20,22,24)]<-H2S.Hrain.unadj[,1]
h2sTab$meanh2s <- h2sTab$nh2s / h2sTab$N * 100
h2sTab$meanLongRain<-c(NA,NA,NA,NA,NA,NA,
      0.5883034,NA,0.6746329,NA, 0.7135688, NA ,
      2.6517274,NA,2.7613850,NA, 2.7933568,NA,
      6.5230230,NA,6.5649825,NA, 6.6059397,NA)
h2sTab$meanHR[c(1:2)]<-HR.rainmean7[,2]
h2sTab$meanHR[c(3:4)]<-HR.rainmean14[,2]
h2sTab$meanHR[c(5:6)]<-HR.rainmean21[,2]
h2sTab$meanHR[c(7:8,13:14,19:20)]<-HR.strat.rainmean7[,3]
h2sTab$meanHR[c(9:10,15:16,21:22)]<-HR.strat.rainmean14[,3]
h2sTab$meanHR[c(11:12,17:18,23:24)]<-HR.strat.rainmean21[,3]
h2sTab$PR[c(1,3,5,7,9,11,13,15,17,19,21,23)] <- "Referent"
h2sTab$PR[c(2,4,6,8,10,12,14,16,18,20,22,24)] <- paste(sprintf("%1.2f",H2S.Hrain.unadj[,5])," (",sprintf("%1.2f",H2S.Hrain.unadj[,6]),", ",sprintf("%1.2f",H2S.Hrain.unadj[,7]),")", sep="")
h2sTab$P[c(2,4,6,8,10,12,14,16,18,20,22,24)]<-sprintf("%1.3f",H2S.Hrain.unadj[,11])
h2sTab$adjPR[c(1,3,5,7,9,11,13,15,17,19,21,23)] <- "Referent"
h2sTab$adjPR[c(2,4,6,8,10,12,14,16,18,20,22,24)] <- paste(sprintf("%1.2f",H2S.Hrain.adj[,5])," (",sprintf("%1.2f",H2S.Hrain.adj[,6]),", ",sprintf("%1.2f",H2S.Hrain.adj[,7]),")", sep="")
h2sTab$adjP[c(2,4,6,8,10,12,14,16,18,20,22,24)]<-sprintf("%1.3f",H2S.Hrain.adj[,8])


 h2sTab[,3]<-paste(sprintf("%1.0f", h2sTab[,3]))
 h2sTab[,4]<-paste(sprintf("%1.0f", h2sTab[,4]))
 h2sTab[,5]<-paste(sprintf("%1.1f", h2sTab[,5]))
 #remove p-value 
 h2sTab<-h2sTab[,-c(9,11)]

cleantable(h2sTab,2)







#Supplimentary table 2: temperature quartiles and h2s prevalence
load("temp.h2s.Rdata")




TempTab<-data.frame(lag=rep(NA,12),
                    Q=rep(NA,12),
                    N=c( 3029, 2976, 3016, 2965, 2932, 3020, 2990, 3007, 2852, 2996, 3069, 2992),
                    ndiar=c( 558, 421, 743, 1223, 504, 392, 777, 1272, 436, 370, 854, 1288),
                    meanH2s=c( 85.18, 83.69, 84.80, 83.40, 82.46, 84.93, 84.08, 84.59, 80.21, 84.26, 86.30, 84.16),
                    meantemp=c( 25.37, 27.46, 29.81, 31.85, 25.31, 27.35, 29.78, 32.05, 25.21, 27.19, 29.64, 32.06),
                    PR=rep(NA,12))


TempTab$lag<-c(rep("One week lag",4),rep("Two week lag",4),rep("Three week lag",4))             
TempTab$Q<-rep(c("~~Quartile 1","~~Quartile 2","~~Quartile 3","~~Quartile 4"),3)
TempTab$N[1:4]<-sprintf("%1.0f", Mean_lag7$quart.n*Mean_lag7$meanH2s)
TempTab$N[5:8]<-sprintf("%1.0f", Mean_lag14$quart.n*Mean_lag14$meanH2s)
TempTab$N[9:12]<-sprintf("%1.0f", Mean_lag21$quart.n*Mean_lag21$meanH2s)
TempTab$ndiar[1:4]<-sprintf("%1.0f", Mean_lag7$sumH2s*Mean_lag7$meanH2s)
TempTab$ndiar[5:8]<-sprintf("%1.0f", Mean_lag14$sumH2s*Mean_lag14$meanH2s)
TempTab$ndiar[9:12]<-sprintf("%1.0f", Mean_lag21$sumH2s*Mean_lag21$meanH2s)
TempTab$meanH2s[1:4]<-Mean_lag7$meanH2s
TempTab$meanH2s[5:8]<-Mean_lag14$meanH2s
TempTab$meanH2s[9:12]<-Mean_lag21$meanH2s
TempTab$meantemp[1:4]<-Mean_lag7$meantemp
TempTab$meantemp[5:8]<-Mean_lag14$meantemp
TempTab$meantemp[9:12]<-Mean_lag21$meantemp
TempTab$PR[c(2:4,6:8,10:12)]<-paste(sprintf("%1.2f",H2S.unadj[,1])," (",sprintf("%1.2f",H2S.unadj[,2]),", ",sprintf("%1.2f",H2S.unadj[,3]),")", sep="")



TempTab[c(1,5,9),c(7)]<-"Referent"
TempTab[,5]<-sprintf("%1.1f",TempTab[,5])

TempTab$N <- as.character(TempTab$N)
TempTab$ndiar <- as.character(TempTab$ndiar)
TempTab$meanH2s <- as.character(TempTab$meanH2s)

cleantable(TempTab[,-1],2)




#Supplimentary table 3: h2s in village sources
load("C:/Users/andre/Dropbox/Trichy analysis/Results/villageH2S.Results.unadjusted.Rdata")
load("C:/Users/andre/Dropbox/Trichy analysis/Results/villageh2s.HRmean.Results.Rdata")
load("C:/Users/andre/Dropbox/Trichy analysis/Results/villageh2s.HRmean-stratified.Results.Rdata")


vil.h2sTab<-data.frame(lag=rep(NA,24),
                    Q=rep(NA,24),
                    N=rep(NA,24),
                    nh2s=rep(NA,24),
                    meanh2s=rep(NA,24),
                    meanLongRain=rep(NA,24),
                    meanHR=rep(NA,24),
                    PR=rep(NA,24),
                    P=rep(NA,24))

H2S.Hrain.unadj<-rbind(
  H2S.Hrain.unadj,
  H2s.unadjLT1,
  H2s.unadjLT2,
  H2s.unadjLT3)


vil.h2sTab$lag<-rep(c(rep("~~One week lag",2),rep("~~Two week lag",2),rep("~~Three week lag",2)),4)             
vil.h2sTab$Q<-rep(c("No heavy rainfall","Heavy rainfall"),12)
vil.h2sTab$N[c(1,3,5,7,9,11,13,15,17,19,21,23)]<-H2S.Hrain.unadj[,3]+H2S.Hrain.unadj[,4]
vil.h2sTab$N[c(2,4,6,8,10,12,14,16,18,20,22,24)]<-H2S.Hrain.unadj[,1]+H2S.Hrain.unadj[,2]
vil.h2sTab$nh2s[c(1,3,5,7,9,11,13,15,17,19,21,23)]<-H2S.Hrain.unadj[,3]
vil.h2sTab$nh2s[c(2,4,6,8,10,12,14,16,18,20,22,24)]<-H2S.Hrain.unadj[,1]
vil.h2sTab$meanh2s <- vil.h2sTab$nh2s / vil.h2sTab$N * 100
vil.h2sTab$meanLongRain<-c(NA,NA,NA,NA,NA,NA,
      0.5883034,NA,0.6746329,NA, 0.7135688, NA ,
      2.6517274,NA,2.7613850,NA, 2.7933568,NA,
      6.5230230,NA,6.5649825,NA, 6.6059397,NA)
vil.h2sTab$meanHR[c(1:2)]<-HR.rainmean7[1:2,2]
vil.h2sTab$meanHR[c(3:4)]<-HR.rainmean14[1:2,2]
vil.h2sTab$meanHR[c(5:6)]<-HR.rainmean21[1:2,2]
vil.h2sTab$meanHR[c(7:8,13:14,19:20)]<-HR.strat.rainmean7[,3]
vil.h2sTab$meanHR[c(9:10,15:16,21:22)]<-HR.strat.rainmean14[,3]
vil.h2sTab$meanHR[c(11:12,17:18,23:24)]<-HR.strat.rainmean21[,3]
vil.h2sTab$PR[c(1,3,5,7,9,11,13,15,17,19,21,23)] <- "Referent"
vil.h2sTab$PR[c(2,4,6,8,10,12,14,16,18,20,22,24)] <- paste(sprintf("%1.2f",H2S.Hrain.unadj[,5])," (",sprintf("%1.2f",H2S.Hrain.unadj[,6]),", ",sprintf("%1.2f",H2S.Hrain.unadj[,7]),")", sep="")
vil.h2sTab$P[c(2,4,6,8,10,12,14,16,18,20,22,24)]<-sprintf("%1.3f",H2S.Hrain.unadj[,11])

 vil.h2sTab[,3]<-paste(sprintf("%1.0f", vil.h2sTab[,3]))
 vil.h2sTab[,4]<-paste(sprintf("%1.0f", vil.h2sTab[,4]))
 vil.h2sTab[,5]<-paste(sprintf("%1.1f", vil.h2sTab[,5]))

  #remove p-value for publication table
 vil.h2sTab<-vil.h2sTab[,-c(9,11)]

  #remove n-prevalent column 
 vil.h2sTab<-vil.h2sTab[,-4]
 
cleantable(vil.h2sTab,2)


