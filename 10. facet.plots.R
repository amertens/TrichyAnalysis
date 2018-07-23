

rm(list=ls())
library(tidyverse)
library(colourpicker)
library(cowplot)
library(stringr)

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

load("spline_plot_facets.Rdata")
load("gamm_plot_facets.Rdata")

    
# 
#     
# ptemp <- plot_grid(temp, p$facets[[1]],  labels = c("A", "B"), ncol = 1)
# ptemp
# 
# 
#     
# prain <- plot_grid(rain, p$facets[[2]],  labels = c("A", "B"), ncol = 1)
# prain
# 
# 
# ggsave(ptemp , file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Figure2.pdf",width=10.4,height=8.32)    
#     
# ggsave(prain , file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Figure3.pdf",width=10.4,height=8.32)      
# 
# 
#   
    
ptemp <- plot_grid(temp, p$facets[[1]], pH2S$facets[[1]],  labels = c("A", "B","C"), ncol = 1)
#ptemp


prain <- plot_grid(rain, p$facets[[2]], pH2S$facets[[2]],  labels = c("A", "B","C"), ncol = 1)
#prain

ggsave(ptemp , file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Figure2.pdf",width=10.4,height=8.32)    
ggsave(prain , file="C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/Figure3.pdf",width=10.4,height=8.32)      
