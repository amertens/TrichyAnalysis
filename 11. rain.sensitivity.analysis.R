
#conditional acceptance update

rm(list=ls())

library(tidyverse)
library(washb)
library(caret)
library(mgcv)
library(stringr)


source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")

try(setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data"))

load("analysis_datasets.Rdata")

d <- d_HRsens


#Check for duplicate study weeks
d %>% group_by(individ, round, stdywk) %>% summarise(n=n()) %>% ungroup() %>% summarise(mean(n))


#Make unique studyweek so AR1 works
d <- d %>%  mutate(stdywk2 = stdywk) %>% arrange(individ, stdywk) %>% group_by(individ, stdywk) %>%
  mutate(stdywk2 = stdywk2 + (row_number()-1)) %>% ungroup() %>% mutate(stdywk=stdywk2) %>%
  subset(., select=-c(stdywk2)) %>% as.data.frame()




#time lag sensitivity analysis, 80% heavy rain
Hrain.unadj.sense80 <- NULL
for(i in c(8, 15, 22)){
  res <- trichy_gamm(d, A=paste0("HeavyRain.lag",i))
  Hrain.unadj.sense80 <- rbind(Hrain.unadj.sense80, res)
}
Hrain.unadj.sense80


Hrain.unadj.sense70 <- NULL
for(i in c(8, 15, 22)){
  res <- trichy_gamm(d, A=paste0("HeavyRain70.lag",i))
  Hrain.unadj.sense70 <- rbind(Hrain.unadj.sense70, res)
}
Hrain.unadj.sense70


Hrain.unadj.sense90 <- NULL
for(i in c(8, 15, 22)){
  res <- trichy_gamm(d, A=paste0("HeavyRain90.lag",i))
  Hrain.unadj.sense90 <- rbind(Hrain.unadj.sense90, res)
}
Hrain.unadj.sense90




#Create plot dataframe
plotdf <- rbind(
  data.frame(Hrain.unadj.sense70, Percentile = c(rep("70th percentile", 3))),
  data.frame(Hrain.unadj.sense80, Percentile = c(rep("80th percentile", 3))),
  data.frame(Hrain.unadj.sense90, Percentile = c(rep("90th percentile", 3))))



plotdf$lag <- as.numeric(str_split_fixed(plotdf$exposure, "Q|lag",2)[,2])
plotdf$lag_label <- rep(c("One week", "Two weeks", "Three weeks"),3)
plotdf$lag_label <- factor(plotdf$lag_label, levels=unique(plotdf$lag_label))

#Create supplimentary figure
theme_set(theme_bw())

setwd("C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/")
png("HRain_sensitivity_figure.png",width=8.5,height=4.249999, units="in", res=800)

ggplot(plotdf, aes(x = factor(lag_label))) +
  geom_pointrange(aes(y=PR, ymin = ci.lb, ymax = ci.ub), color="grey30", size = 0.75) +
  # geom_point(aes(y=PR), size = 1) +
  # geom_linerange(aes( ymin=ci.lb, ymax=ci.ub), alpha=0.5, size = 3) +
  geom_hline(yintercept = 1) +
  facet_grid( . ~ Percentile) +
  scale_y_log10(breaks=c(0.5, 0.75, 1, 1.5, 2), labels=c(0.5, 0.75, 1, 1.5, 2)) +
  labs(y = "Prevalence ratio",
       x = "Length of lag") +
  #scale_x_continuous(breaks=c(7,14,21,28)) +
  theme(strip.background = element_blank())

dev.off()






#Old figure
# 
# rm(list=ls())
# 
# library(tidyverse)
# library(washb)
# library(caret)
# library(mgcv)
# library(stringr)
# 
# 
# source("C:/Users/andre/Documents/TrichyAnalysis/0. TrichyFunctions.R")
# 
# try(setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data"))
# 
# load("analysis_datasets.Rdata")
# 
# d <- d_HRsens
# 
# 
# #Check for duplicate study weeks
# d %>% group_by(individ, round, stdywk) %>% summarise(n=n()) %>% ungroup() %>% summarise(mean(n))
# 
# 
# #Make unique studyweek so AR1 works
# d <- d %>%  mutate(stdywk2 = stdywk) %>% arrange(individ, stdywk) %>% group_by(individ, stdywk) %>%
#   mutate(stdywk2 = stdywk2 + (row_number()-1)) %>% ungroup() %>% mutate(stdywk=stdywk2) %>%
#   subset(., select=-c(stdywk2)) %>% as.data.frame()
# 
# 
# # HR1<-trichy_gamm(d, A="HeavyRain.lag7")
# # HR1
# 
# 
# 
# #time lag sensitivity analysis, 80% heavy rain
# Hrain.unadj.sense80 <- NULL
# for(i in 7:28){
#   res <- trichy_gamm(d, A=paste0("HeavyRain.lag",i))
#   Hrain.unadj.sense80 <- rbind(Hrain.unadj.sense80, res)
# }
# Hrain.unadj.sense80
# 
# 
# Hrain.unadj.sense70 <- NULL
# for(i in 7:28){
#   res <- trichy_gamm(d, A=paste0("HeavyRain70.lag",i))
#   Hrain.unadj.sense70 <- rbind(Hrain.unadj.sense70, res)
# }
# Hrain.unadj.sense70
# 
# 
# Hrain.unadj.sense90 <- NULL
# for(i in 7:28){
#   res <- trichy_gamm(d, A=paste0("HeavyRain90.lag",i))
#   Hrain.unadj.sense90 <- rbind(Hrain.unadj.sense90, res)
# }
# Hrain.unadj.sense90
# 
# Hrain.unadj.sense90a <- NULL
# for(i in 7:28){
#   res <- trichy_gamm(d, A=paste0("HeavyRain90a.lag",i))
#   Hrain.unadj.sense90a <- rbind(Hrain.unadj.sense90a, res)
# }
# Hrain.unadj.sense90a
# 
# 
# #Save results
# save(Hrain.unadj.sense80, Hrain.unadj.sense70, Hrain.unadj.sense90, Hrain.unadj.sense90a, file="C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data/HR_sens_res.Rdata")
# 
# 
# #Create plot dataframe
# plotdf <- rbind(
#   #data.frame(Hrain.unadj.sense90a, Percentile = c(rep("90th percentile of all days", 22))),
#   data.frame(Hrain.unadj.sense70, Percentile = c(rep("70th percentile", 22))),
#   data.frame(Hrain.unadj.sense80, Percentile = c(rep("80th percentile", 22))),
#   data.frame(Hrain.unadj.sense90, Percentile = c(rep("90th percentile", 22))))
# 
#   
# plotdf$lag<-plotdf$lagjitter<-0:21
# plotdf$lagjitter[1:22] <- plotdf$lagjitter[1:22] - 0.3
# plotdf$lagjitter[45:66] <- plotdf$lagjitter[45:66] + 0.3
# 
# plotdf$lag <- as.numeric(str_split_fixed(plotdf$exposure, "Q|lag",2)[,2])
# 
# 
# #
# plotdf %>% group_by(Percentile) %>% filter(lag <22) %>% 
#   summarize(mean( PR), var(PR))
# 
# #Add 7 days to the lag to match up plots with the 1,2,3 wk 
# #terminology in manuscript
# #plotdf$lag <- plotdf$lag+7
# 
# 
# 
# 
# #Create supplimentary figure
# theme_set(theme_bw())
# 
# setwd("C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/")
# #pdf("HRain_sensitivity_figure.pdf",width=10,height=5)
# png("HRain_sensitivity_figure.png",width=8.5,height=4.249999, units="in", res=800)
# 
# ggplot(plotdf, aes(x = lag)) +
#   geom_pointrange(aes(y=PR, ymin = ci.lb, ymax = ci.ub), color="grey30") +
#     geom_hline(yintercept = 1) +
#     facet_grid( . ~ Percentile) +
#     scale_y_log10(breaks=c(0.5, 0.75, 1, 1.5, 2), labels=c(0.5, 0.75, 1, 1.5, 2)) +
#     labs(y = "Prevalence ratio",
#          x = "Lag period (days)") +
#     scale_x_continuous(breaks=c(7,14,21,28)) +
#     theme(strip.background = element_blank())
# 
# dev.off()
# 
# 
