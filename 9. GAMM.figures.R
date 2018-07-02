
rm(list=ls())
library(tidyverse)
library(colourpicker)
library(cowplot)
library(stringr)

#source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")
    #plot parameters
    theme_set(theme_bw())
    cols=c("(ref.)"="#919191",
           "Unstratified"="#000000",
           "1v2"="#ED8021",
           "1v3"="#E15324",
           "1v4"="#D62728",
           "Low"="#56B4E9",
           "Medium"="#4C8FE5",
           "High"="#426AE3")


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

PR_plotfun_strat <- function(df, xlab="",title="", yticks=c(0.0625,0.125,0.25,0.5,1,2,4,8)){

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
      facet_wrap(lag~int,  scales = "fixed", nrow = 1) +
      ggtitle(title)
  return(p)
}


trichy_multiplot <- function(d, titlelab=""){


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
        ptemp <- PR_plotfun(df=d1, xlab="Quartile contrast", title=paste0("Temperature - diarrhea association",titlelab))
      }

      if(nrow(d2)>0){
        prain <- PR_plotfun(df=d2, xlab="Long-term rainfall strata", title=paste0("Heavy rainfall - diarrhea association",titlelab))
      }

      if(nrow(d3)>0){
        ptempH2s <- PR_plotfun(df=d3, xlab="Quartile contrast", title=paste0("Temperature - drinking water H2S association",titlelab))
      }

      if(nrow(d4)>0){
        pH2S <- PR_plotfun(df=d4, xlab="Long-term rainfall strata", title=paste0("Heavy rainfall - drinking water H2S association",titlelab))
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



load("C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_results.Rdata")
#load("C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_plot_dfs.Rdata")



#primary unadjusted
res_prim <- rbind(T1, T2, T3, HR1, HR2, HR3, HR1_strat, HR2_strat, HR3_strat)

#primary adjusted
res_prim_adj <- rbind(T1_adj, T2_adj, T3_adj, 
                      HR1_adj, HR2_adj, HR3_adj,
                      HR1_strat_adj, HR2_strat_adj, HR3_strat_adj)

#h2s
res_H2S <- rbind( h2s.T1, h2s.T2, h2s.T3,
 h2s.HR1, h2s.HR2, h2s.HR3,
 h2s.HR1_strat, h2s.HR2_strat, h2s.HR3_strat)

res_H2S_adj <- rbind(h2s.HR1_adj, h2s.HR2_adj, h2s.HR3_adj,
h2s.HR1_strat_adj, h2s.HR2_strat_adj, h2s.HR3_strat_adj)

#subgroups



res_minmax_temp <- rbind(maxT1, maxT2, maxT3, minT1,  minT2, minT3)
res_wpi <- rbind(T1_wpi, T2_wpi, T3_wpi, HR1_wpi, HR2_wpi, HR3_wpi, HR1_wpi_strat, HR2_wpi_strat, HR3_wpi_strat)
res_wpi_adj<- rbind(T1_wpi_adj, T2_wpi_adj, T3_wpi_adj, HR1_wpi_adj, HR2_wpi_adj, HR3_wpi_adj)
 
#90% rain
res_HR90 <- rbind(HR1_90, HR2_90, HR3_90, HR1_strat_90, HR2_strat_90, HR3_strat_90)







p <- trichy_multiplot(res_prim, titlelab = ", unadjusted")
#p
 

p_adj <- trichy_multiplot(res_prim_adj, titlelab = ", adjusted")
#p_adj
 

pH2S <- trichy_multiplot(res_H2S, titlelab = ", unadjusted")
#pH2s

pH2S_adj <- trichy_multiplot(res_H2S_adj, titlelab = ", adjusted")

p_min <- trichy_multiplot(res_minmax_temp[grepl("min",res_minmax_temp$exposure),], titlelab = ", \nminimum daily temperature")
p_max <- trichy_multiplot(res_minmax_temp[grepl("max",res_minmax_temp$exposure),], titlelab = ", \nmaximum daily temperature")

p_minmax <- plot_grid(p_min, p_max, labels = c("A", "B"), ncol = 1)

p90 <- trichy_multiplot(res_HR90, titlelab = ", \n90th percentile of all days threshold")



#Stratified by intervention plots
d<-res_wpi
d$int <- str_split_fixed(d$strat, "_",2)[,2]
d$int[d$strat=="0"]<-"0"
d$int[d$strat=="1"]<-"1"
d$int[d$int =="0"] <- "Control village"
d$int[d$int =="1"] <- "Intervention village"


d$strat <- c(rep(c("2","3","4"),6),rep("Unstratified",6),rep(rep(c("T1","T2","T3"),each=2),3))

titlelab = ", \nstratified by intervention arm"

   #plot parameters
    theme_set(theme_bw())
    cols=c("(ref.)"="#919191",
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
      temp_blank <- d[1:6,]
      temp_blank$quartile <- 1
      temp_blank$strat <- "(ref.)"
      temp_blank$lag <- rep(c("7","14","21"), 2)
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
        ptemp <- PR_plotfun_strat(df=d1, xlab="Quartile contrast", title=paste0("Temperature - diarrhea association",titlelab))
      }

      if(nrow(d2)>0){
        prain <- PR_plotfun_strat(df=d2, xlab="Long-term rainfall strata", title=paste0("Heavy rainfall - diarrhea association",titlelab))
      }

      if(nrow(d3)>0){
        ptempH2s <- PR_plotfun_strat(df=d3, xlab="Quartile contrast", title=paste0("Temperature - drinking water H2S association",titlelab))
      }

      if(nrow(d4)>0){
        pH2S <- PR_plotfun_strat(df=d4, xlab="Long-term rainfall strata", title=paste0("Heavy rainfall - drinking water H2S association",titlelab))
      }

    plist <- list(ptemp, prain, ptempH2s,pH2S)
    names(plist) <- seq_along(plist)
    plist[sapply(plist, is.null)] <- NULL

    p_int <- plot_grid(plist[[1]], plist[[2]],  labels = c("A", "B"), ncol = 1)




  
pdf(file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Trichy_GAMM_figures.pdf")

  p
  p_adj
  pH2S
  pH2S_adj
  p_int
  p_minmax
  p90
  
dev.off()

  



#Wider intervention- stratified figure for reviewer

pdf(file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Trichy_intervention_stratified_figure.pdf", width=12)

  p_int
  
dev.off()


