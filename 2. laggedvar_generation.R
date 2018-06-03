
rm(list=ls())

library(dplyr)
library(foreign)

###Generating lagged weather variables
setwd("C:/Users/andre/Dropbox/Trichy analysis/Data")
df<-read.dta("Trichy_Weather_Dec07-Apr09_formatted.dta")

head(df)

#Order data
df<-arrange(df, year, month, day) %>% 
    subset(., select=c(year, month, day, maxtemp, mintemp, rain))


#Generate average temp
df$avetemp<-(df$maxtemp+df$mintemp)/2



#Generate lagged ave temp
lag.n <- function(df,name,var, n) {
    varname <- paste(name, n , sep=".")
    df[[varname]] <- with(df,  lag(var, n))
    return(df)
}


dayslag<-100

for(i in 1:dayslag){
  df<-lag.n(df=df,name="at",var=df$avetemp,n=i)
}


#Generate lagged min and max temp
for(i in 1:dayslag){
  df<-lag.n(df=df,name="maxt",var=df$maxtemp,n=i)
}
for(i in 1:dayslag){
  df<-lag.n(df=df,name="mint",var=df$mintemp,n=i)
}



#Generate lagged rain
for(i in 1:dayslag){
  df<-lag.n(df=df,name="rain",var=df$rain,n=i)
}






#save data
weather<-df

save(weather, file="C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data/LaggedWeather.Rdata")

