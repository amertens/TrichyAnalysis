



#Trichy analysis figure creation
rm(list=ls())
library(xtable)
library(tidyverse)
library(colourpicker)
library(cowplot)

source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")


#plot parameters
# colfunc <- colorRampPalette(c("#ED8021","#D62728"))
# colfunc(3)
# colourPicker()
tempcol <- c("#000000","#ED8021", "#E15324", "#D62728")
raincol <- c("#000000",  "#56B4E9", "#4C8FE5", "#426AE3")
theme_set(theme_bw())



#Load results
setwd("C:/Users/andre/Dropbox/Trichy analysis/Results/")
load("temp.Results.unadjusted.Rdata")
load("temp.Results.adjusted.Rdata")
load("rain.Results.unadjusted.Rdata")
load("rain.Results.adjusted.Rdata")

#sensitivity results
load("SE.clustering.results.Rdata")
load("treatment_subgroup_results.Rdata")


#------------------------------
# Format results data.frames
#------------------------------

plotdf<-plotdf_format(temp.unadj, Hrain.unadj, Hrain.unadjLT1,Hrain.unadjLT2, Hrain.unadjLT3, H2S.Hrain.unadj, H2s.unadjLT1, H2s.unadjLT2, H2s.unadjLT3)


plotdf_seclust<-plotdf_format(SEclust_res[[1]],
                              SEclust_res[[2]],
                              SEclust_res[[3]],
                              SEclust_res[[4]],
                              SEclust_res[[5]],
                              SEclust_res[[6]],
                              SEclust_res[[7]],
                              SEclust_res[[8]],
                              SEclust_res[[9]])


plotdf_ctrl<-plotdf_format(ctrl_subgroup_res[[1]],
                              ctrl_subgroup_res[[2]],
                              ctrl_subgroup_res[[3]],
                              ctrl_subgroup_res[[4]],
                              ctrl_subgroup_res[[5]],
                              ctrl_subgroup_res[[6]],
                              ctrl_subgroup_res[[7]],
                              ctrl_subgroup_res[[8]],
                              ctrl_subgroup_res[[9]])


plotdf_trt<-plotdf_format(trt_subgroup_res[[1]],
                              trt_subgroup_res[[2]],
                              trt_subgroup_res[[3]],
                              trt_subgroup_res[[4]],
                              trt_subgroup_res[[5]],
                              trt_subgroup_res[[6]],
                              trt_subgroup_res[[7]],
                              trt_subgroup_res[[8]],
                              trt_subgroup_res[[9]])
                              

  
#------------------------------
# Create primary results plot
#------------------------------
  
  
  p <- trichy_panel_plot(plotdf$raindf, plotdf$H2Sdf, plotdf$tempdf)
  p
  
  ggsave(p, file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Figure3.pdf")
  
  
#------------------------------
# Create sensitivity results plot
#------------------------------
  
  p_seclust <- trichy_panel_plot(plotdf_seclust$raindf, plotdf_seclust$H2Sdf, plotdf_seclust$tempdf)

  ggsave(p_seclust, file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Clusterid_sensitivity_figure.pdf")

  p_ctrl <- trichy_panel_plot(plotdf_ctrl$raindf, plotdf_ctrl$H2Sdf, plotdf_ctrl$tempdf)
  p_trt <- trichy_panel_plot(plotdf_trt$raindf, plotdf_trt$H2Sdf, plotdf_trt$tempdf)
  
  
#------------------------------
# Compare plots from different 
# sensitivity analyses
#------------------------------
      
  
pdf(file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Sensitivity_figures.pdf")

  p
  p_seclust
  p_ctrl
  p_trt
  
dev.off()
  