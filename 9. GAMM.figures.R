

rm(list=ls())
library(tidyverse)
library(colourpicker)
library(cowplot)
library(stringr)

    #plot parameters
    theme_set(theme_bw())
    cols=c("(ref.)"="#919191",
           "Q1"="#919191",
           "Unstratified"="#000000",
           "Q2"="#ED8021",
           "Q3"="#E15324",
           "Q4"="#D62728",
           "Low"="#56B4E9",
           "Medium"="#4C8FE5",
           "High"="#426AE3")


#------------------------------
# Function to combine plots 
# into one figure
#------------------------------
PR_plotfun <- function(df, xlab="",title="", yticks){

  df$strat <- factor(df$strat, levels=unique(df$strat))
  df$lag <- factor(df$lag, levels=unique(df$lag))

  p <-ggplot(df, aes(x=strat)) +
      geom_point(aes(y=PR, fill=strat, color=strat), size = 1) +
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
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      facet_wrap(~lag,  scales = "fixed") +
      ggtitle(title)
  return(p)
}



trichy_multiplot <- function(d){


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
          if(d$outcome!="H2S"){temp_blank$lag <- c("7","14","21")}else{
            temp_blank$lag <- c("1","7","14")
          }
      temp_blank$PR <- 1
      temp_blank$ci.lb <- temp_blank$ci.ub <-temp_blank$Xvar <- NA
      d <- rbind(temp_blank,d)
    }

      #Clean up labels
    if(d$outcome!="H2S"){
    d$lag[d$lag=="7"] <- "One week lag"
    d$lag[d$lag=="14"] <- "Two week lag"
    d$lag[d$lag=="21"] <- "Three week lag"
    }else{
    d$lag[d$lag=="1"] <- "One week lag"
    d$lag[d$lag=="7"] <- "Two week lag"
    d$lag[d$lag=="14"] <- "Three week lag"
    }
    d$strat[d$strat=="2"] <- "Q2"
    d$strat[d$strat=="3"] <- "Q3"
    d$strat[d$strat=="4"] <- "Q4"
    if("Q2" %in% d$strat){d$strat[d$strat=="(ref.)"]<-"Q1"}

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
        ptemp <- PR_plotfun(df=d1, xlab="Quartile contrast", title=paste0("Prevalence ratios of diarrhea across quartiles of weekly mean temperature"), yticks=c(0.5,1,2,4,8))
      }

      if(nrow(d2)>0){
        prain <- PR_plotfun(df=d2, xlab="60-day rain trend", title=paste0("Heavy rainfall - diarrhea association, unstratified and stratified by 60-day rain trends"), yticks=c(0.125,0.25,0.5,1,2,4))
      }

      if(nrow(d3)>0){
        ptempH2s <- PR_plotfun(df=d3, xlab="Quartile contrast", title=expression('Prevalence ratios of drinking water H'[2]*'S across quartiles of weekly mean temperature'), yticks=c(0.75, 0.87, 1, 1.15, 1.33))
      }

      if(nrow(d4)>0){
        pH2S <- PR_plotfun(df=d4, xlab="60-day rain trend", title=expression('Heavy rainfall - drinking water H'[2]*'S association, unstratified and stratified by 60-day rain trends'), yticks=c(0.75, 0.87, 1, 1.15, 1.33))
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

 return(list(plot=p, facets=plist))
}



load("C:/Users/andre/Dropbox/Trichy analysis/Results/GAMM_results.Rdata")



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

res_H2S_adj <- rbind(h2s.T1_adj, h2s.T2_adj, h2s.T3_adj,
                     h2s.HR1_adj, h2s.HR2_adj, h2s.HR3_adj,
h2s.HR1_strat_adj, h2s.HR2_strat_adj, h2s.HR3_strat_adj)





#primary plots
p <- trichy_multiplot(res_prim_adj)


pH2S <- trichy_multiplot(res_H2S_adj)



setwd("C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/")
save(p, pH2S, file ="gamm_plot_facets.Rdata")





#--------------------------------
# Make manuscript plots by combining 
# with spline plots
#--------------------------------

library(cowplot)
library(stringr) 

#plot parameters
theme_set(theme_bw())
cols=c("(ref.)"="#919191",
       "Q1"="#919191",
       "Unstratified"="#000000",
       "Q2"="#ED8021",
       "Q3"="#E15324",
       "Q4"="#D62728",
       "Low"="#56B4E9",
       "Medium"="#4C8FE5",
       "High"="#426AE3")

load("spline_plot_facets.Rdata")
load("gamm_plot_facets.Rdata")


ptemp <- plot_grid(temp, p$facets[[1]], pH2S$facets[[1]],  labels = c("A", "B","C"), ncol = 1)

prain <- plot_grid(rain, p$facets[[2]], pH2S$facets[[2]],  labels = c("A", "B","C"), ncol = 1)

ggsave(ptemp , file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Figure3.pdf",width=10.4,height=8.32)    
ggsave(prain , file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Figure4.pdf",width=10.4,height=8.32)      


