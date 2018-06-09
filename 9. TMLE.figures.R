



#Trichy analysis table printing
rm(list=ls())
library(xtable)
library(tidyverse)
library(colourpicker)
library(cowplot)

theme_set(theme_bw())

#Forest plot Color palette
cbPalette <- c( "#56B4E9" , rep("#999999",40))

#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
  "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cols <- cbPalette[c(1,3,7)]

#Temp pallete
tempcol <- c("#000000","#ED8021", "#E15324", "#D62728")
  
#Rain pallete
raincol <- c("#000000",  "#56B4E9", "#4C8FE5", "#426AE3")
  

# colfunc <- colorRampPalette(c("#ED8021","#D62728"))
# colfunc(3)

#colourPicker()


scaleFUN <- function(x) sprintf("%.2f", x)






#Load results
setwd("C:/Users/andre/Dropbox/Trichy analysis/Results/")
load("temp.Results.unadjusted.Rdata")
load("temp.Results.adjusted.Rdata")
load("rain.Results.unadjusted.Rdata")
load("rain.Results.adjusted.Rdata")


temp.unadj
Hrain.unadj
Hrain.unadjLT1
Hrain.unadjLT2
Hrain.unadjLT3
H2S.Hrain.unadj
H2s.unadjLT1
H2s.unadjLT2
H2s.unadjLT3


raindf <- data.frame(exposure=rep("Heavy rain",12),
             strata=c(rep("Unstratified",3),
             rep("Low",3),
             rep("Medium",3),
             rep("High",3)),
             x=as.character(c(1,1,1,2,2,2,3,3,3,4,4,4)),
             lag=rep(c("One week lag", "Two week lag", "Three week lag"),4),
              rbind(
              Hrain.unadj,
              Hrain.unadjLT1,
              Hrain.unadjLT2,
              Hrain.unadjLT3))
raindf<- raindf %>% rename(PR=PR.)
raindf$strata <- factor(raindf$strata, levels=unique(raindf$strata))
raindf$lag <- factor(raindf$lag, levels=unique(raindf$lag))

raindf<- subset(raindf, select=-c(a,b,c,d, logPR,se.logPR,Z,p))

#Add blank rows for reference levels
temp.unadj <- rbind(NA,temp.unadj[1:3,], NA,temp.unadj[4:6,], NA,temp.unadj[7:9,])
rownames(temp.unadj) <- NULL
tempdf <- data.frame(exposure=rep("Temperature",12),
             strata=rep(c("(ref.)","1v2","1v3","1v4"),3),
                          x=as.character(c(1,2,3,4,1,2,3,4,1,2,3,4)),
             lag=c(rep("One week lag",4), rep("Two week lag",4), rep("Three week lag",4)),
                         temp.unadj) %>% 
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
              H2S.Hrain.unadj,
              H2s.unadjLT1,
              H2s.unadjLT2,
              H2s.unadjLT3))
H2Sdf<- H2Sdf %>% rename(PR=PR.)
H2Sdf$strata <- factor(H2Sdf$strata, levels=unique(H2Sdf$strata))
H2Sdf$lag <- factor(H2Sdf$lag, levels=unique(H2Sdf$lag))

H2Sdf<- subset(H2Sdf, select=-c(a,b,c,d, logPR,se.logPR,Z,p))



#--------------------------------------
#plot rain
#--------------------------------------

  # Rainplot<-ggplot(data=raindf) + 
  #   labs( x = "Treatment arm", y = "Prevalence ratio") +
  #   geom_hline(yintercept = 1) +
  #   scale_y_continuous(breaks=c(0.125,0.25,0.5,1,2,4,8), trans='log10') +
  #   #coord_flip(ylim = c(0.25, 4,8)) +
  #   geom_pointrange( mapping=aes(x=strata, y= PR., ymin=ci.lb, ymax=ci.ub , colour=strata)) +
  #   facet_wrap(~lag, scales = "free_x") +
  #   theme(panel.border = element_blank(), 
  #         strip.background = element_blank(), legend.position = "none") + 
  #   ggtitle("Effect of heavy rain on diarrhea prevalence,\nunstratified and stratified by long-term rain trend") +theme(legend.position="none")
  # Rainplot

  
    yticks=c(0.125,0.25,0.5,1,2,4,8)
  
  
prain<-ggplot(raindf, aes(x=strata)) + 
  geom_point(aes(y=PR, fill=strata, color=strata), size = 4) +
  geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strata),
                 alpha=0.5, size = 3) +
  labs(x = "Long-term rainfall strata", y = "Prevalence Ratio") +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim=range(yticks)) +
  scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
  scale_fill_manual(values=rep(raincol,4)) +
  scale_colour_manual(values=rep(raincol,4)) +
  scale_size_continuous(range = c(0.5, 1))+
  theme(strip.background = element_blank(),
    legend.position="none",
    strip.text.x = element_text(size=12),
    axis.text.x = element_text(size=12)) +
  facet_wrap(~lag,  scales = "fixed") +
  ggtitle("Rain")
prain



#--------------------------------------
#plot H2S-rain
#--------------------------------------
  
    yticks=c(0.125,0.25,0.5,1,2,4,8)
  
  
pH2S<-ggplot(H2Sdf, aes(x=strata)) + 
  geom_point(aes(y=PR, fill=strata, color=strata), size = 4) +
  geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strata),
                 alpha=0.5, size = 3) +
  labs(x = "Long-term rainfall strata", y = "Prevalence Ratio") +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim=range(yticks)) +
  scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
  scale_fill_manual(values=rep(raincol,4)) +
  scale_colour_manual(values=rep(raincol,4)) +
  scale_size_continuous(range = c(0.5, 1))+
  theme(strip.background = element_blank(),
    legend.position="none",
    strip.text.x = element_text(size=12),
    axis.text.x = element_text(size=12)) +
  facet_wrap(~lag,  scales = "fixed") +
  ggtitle("H2S")
pH2S

#--------------------------------------
#plot temperature
#--------------------------------------


    yticks=c(0.125,0.25,0.5,1,2,4,8)
  
  
ptemp<-ggplot(tempdf, aes(x=strata)) + 
  geom_point(aes(y=PR, fill=strata, color=strata), size = 4) +
  geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strata),
                 alpha=0.5, size = 3) +
  labs(x = "Quartile contrast", y = "Prevalence Ratio") +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim=range(yticks)) +
  scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
  scale_fill_manual(values=rep(tempcol,4)) +
  scale_colour_manual(values=rep(tempcol,4)) +
  scale_size_continuous(range = c(0.5, 1))+
  theme(strip.background = element_blank(),
    legend.position="none",
    strip.text.x = element_text(size=12),
    axis.text.x = element_text(size=12)) +
  facet_wrap(~lag,  scales = "fixed") +
  ggtitle("Temperature")
ptemp

#Change the colors to be increasingly red


#--------------------------------------
#plot rain + temperature
#--------------------------------------

df <- rbind(raindf, tempdf)

  yticks=c(0.125,0.25,0.5,1,2,4,8)
  

  
# p<-ggplot(df, aes(x=x)) + 
#   geom_point(aes(y=PR, fill=strata, color=strata), size = 4) +
#   geom_linerange(aes( ymin=ci.lb, ymax=ci.ub, color=strata),
#                  alpha=0.5, size = 3) +
#   labs(x = "Risk factor level", y = "Prevalence Ratio") +
#   geom_hline(yintercept = 1) +
#   coord_cartesian(ylim=range(yticks)) +
#   scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
#   scale_x_discrete(labels=df$strata) +
#   scale_fill_manual(values=rep(tableau10,4)) +
#   scale_colour_manual(values=rep(tableau10,4)) +
#   scale_size_continuous(range = c(0.5, 1))+
#   theme(strip.background = element_blank(),
#     legend.position="none",
#     strip.text.x = element_text(size=12),
#     axis.text.x = element_text(size=12)) +
#   #facet_wrap(exposure~lag,  scales = "free_x", strip.position="right") +
#   facet_grid(exposure~lag,  scales = "free_x") +
#   ggtitle("")
# p
#   

#make two plots then cowplot.


#Maybe make the unstratified grey, with low, med, and high strata becoming more blue

#Then the Temp can be faceted below, with blank reference and increasingly orange for 1v2, 1v3, 1v4

  
  
  
#------------------------------
# Combine plots into one figure
#------------------------------
  
  
  p <- plot_grid(ptemp, prain, pH2S, labels = c("A", "B","C"), ncol = 1)
  p
  
  ggsave(p, file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Figure3.pdf")
  
  
      