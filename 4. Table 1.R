

rm(list=ls())
try(setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data"))
sink("./table1.results.txt")
library(tidyverse)
library(haven)


mean_sd <- function(x){
  cat("\n", round(mean(x, na.rm=T),1),
      "  ", round(sd(x, na.rm=T),1),"\n") 
}
N_perc <- function(x){
  print(table(x))
  cat("\n")
  for(i in unique(x)){
    cat(i,":     ")
    y<-ifelse(x==i,1,0)
    cat(sum(y, na.rm=T),
        " (", round(mean(y*100, na.rm=T),1),")\n", sep = "") 
  }
}


load("analysis_datasets.Rdata")

#Load stata dataset for attribute labels
svy <- read_dta(file = "C:/Users/andre/Dropbox/Trichy analysis/Data/trichy_long.dta") 

# num obs
nrow(d)

# num children
length(unique(d$individ))

#Num visits
d %>% group_by(individ) %>% summarize(N=n()) %>% 
  ungroup() %>% summarize(n(), mean(N), sd(N))
mean_sd(as.numeric(d$round))

#Proportion of water samples taken by round
d %>% group_by(hhid,round) %>% slice(1) %>% ungroup() %>% 
  group_by(round) %>% summarize(N=n(), x=round(sum(!is.na(H2S))/900*100,2))

# child age
#Does it make sense to include child age in the table 1?
mean_sd(d$age)

# num households
length(unique(d$hhid))


# Pre-intervention covariates:	
#Subset data to first obs per child
dfull <- d
d <- d %>% arrange(round)%>% group_by(individ) %>% slice(1)

#   Received Sanitation Intervention	 (49.2)
N_perc(d$wpi)


# Child
# Female	(49.4)
N_perc(d$sex)

# Breastfeeding	(25.8)
N_perc(d$bfcur)
N_perc(dfull$bfcur)

# Maternal
d <- dfull %>% arrange(round)%>% group_by(hhid) %>% slice(1)

# Age	26.5 ± 
mean_sd(d$momage)

# Literate	(83.8)
cat("mom literate\n")
N_perc(d$pclit)

# No education	(13.0)
# Primary School	(13.2)
# Middle School	(27.2)
# High school	(32.8)
# Any secondary	(13.7)
N_perc(d$pceduc)


# Works	(47.8)
cat("mom works\n")
N_perc(d$momwork)

# Household	
# Number of people in household	4.9 ±
mean_sd(as.numeric(as.character(d$totp)))

# Number of children in household	X.X ±
d %>% group_by(hhid) %>% summarize(x=length(unique(individ))) %>% ungroup() %>% summarize(mean(x), sd(x))

# Has electricity	(90.1)
cat("electricity\n")
N_perc(d$elec)

# Has a dirt floor	(31.4)
cat("dirt floor\n")
N_perc(d$soil)

# Owns home
cat("owns home\n")
N_perc(d$ownhome)

# Owns land	(94.3)
cat("owns land\n")
N_perc(d$land)

# Animals
cat("buffalo\n")
N_perc(d$buff)
cat("cow\n")
N_perc(d$cow)
cat("ox\n")
N_perc(d$ox)
cat("calf\n")
N_perc(d$calf)
cat("goat\n")
N_perc(d$goat)
cat("chicken\n")
N_perc(d$chick)
cat("dog or cat\n")
N_perc(d$dogcat)


# Drinking water	
# Private Tap	(30.0)
# Public Tap	(62.9)
# Private Well	2.9)
# Public Well	(4.2)
attributes(svy$primwat)$`label`
N_perc(d$primwat)


# Soap observed by handwashing station	(49.3)


# Sanitation	
# Daily defecating in open	(83.4)
attributes(svy$sanOD)$`label`
N_perc(d$sanOD)

# Latrine	(41.0)
attributes(svy$sanlat)$`label`
N_perc(d$sanlat)

# Latrine has soap	(47.1)
# Smell of feces in home	(10.3)
# Feces observed in home	(27.4)

# Household assets
attributes(svy$thatch)$`label`
N_perc(d$thatch)
 
attributes(svy$tv)$`label`
N_perc(d$tv)

attributes(svy$cell)$`label`
N_perc(d$cell)

attributes(svy$bank)$`label`
N_perc(d$bank)

attributes(svy$moto)$`label`
N_perc(d$moto)

attributes(svy$bike)$`label`
N_perc(d$bike)

attributes(svy$mosq)$`label`
N_perc(d$mosq)

attributes(svy$stove)$`label`
N_perc(d$stove)

attributes(svy$fuel1)$`label`
N_perc(d$fuel1)


attributes(svy$anygrp)$`label`
N_perc(d$anygrp)

attributes(svy$shgw)$`label`
N_perc(d$shgw)

attributes(svy$agri)$`label`
N_perc(d$agri)


attributes(svy$elec)$`label`
N_perc(d$elec)

attributes(svy$land)$`label`
N_perc(d$land)

attributes(svy$sc)$`label`
N_perc(d$sc)

attributes(svy$kitchen)$`label`
N_perc(d$kitchen)

attributes(svy$kitchvent)$`label`
N_perc(d$kitchvent)

attributes(svy$hwflies)$`label`
N_perc(d$hwflies)

attributes(svy$hwwater)$`label`
N_perc(d$hwwater)

attributes(svy$hwsoap)$`label`
N_perc(d$hwsoap)

attributes(svy$hwash)$`label`
N_perc(d$hwash)

attributes(svy$hwtow)$`label`
N_perc(d$hwtow)

attributes(svy$hwsink)$`label`
N_perc(d$hwsink)


d <- d %>% group_by(vilid) %>% slice(1) %>% as.data.frame()

print("village OD")
mean_sd(d$villOD)

sink()
