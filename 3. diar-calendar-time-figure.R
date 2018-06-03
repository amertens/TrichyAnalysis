#---------------------------------
# diar-calendar-time.R
#
# summarize diarrhea prevalence
# by calendar time in the trichy study
# and overlay rainfall and temperature trends
#---------------------------------

#---------------------------------------
# preamble
#---------------------------------------
rm(list=ls())
library(washb)
library(lubridate)
library(scales)
library(dplyr)

#---------------------------------------
# Load the analysis dataset
#---------------------------------------
setwd("C:/Users/andre/Dropbox/Trichy analysis/Data/Cleaned data")
load("bl_covariates.Rdata")
load("LaggedWeather.Rdata")

weather<-weather %>% subset(.,select=c(year, month, day, rain, avetemp))
head(weather)

df<-df %>% subset(.,select=c(intdate, diar7d, individ))
head(df)



#Set to date format for merge
weather$intdate<-as.Date(paste(weather$month,weather$day,weather$year, sep="/"),"%m/%d/%Y")
weather<-subset(weather, select=-c(year,month,day))

ad<-merge(df,weather, by="intdate", all.x = F, all.y = F)  
head(ad)


class(ad$intdate)



#---------------------------------------
# create a calendar time aggregator
# by month
#---------------------------------------
ad$caldate <- as.Date(ad$intdate,format="%d%b%Y")

head(ad$caldate)
ad$my <- floor_date(ad$caldate,"month")

#create by week for weather
weather$wy <- floor_date(weather$intdate,"week")
weather$my <- floor_date(weather$intdate,"month")
weather$yr <- floor_date(weather$intdate,"year")


#---------------------------------------
# monthly mean diarrhea w/ robust 95% SEs
#---------------------------------------
mys <- unique(ad$my)
diar <- sapply(mys,
                 function(x) {
                   washb_mean(ad$diar7d[ad$my==x],
                              id=ad$individ[ad$my==x],
                              print=F
                              )
                 }
                 )
diar <- data.frame(t(diar))
diar$month <- as.Date(mys)
diar$year <- floor_date(diar$month,"year")
colnames(diar) <- c("n","mean","sd","se","lb","ub","month","year")
diar <- diar[order(diar$month),c("month","year","n","mean","sd","se","lb","ub")]


  xtics <- unique(diar$month) #Save every month 
  xticsyr <- unique(diar$year) 

# summarize dist of obs per month
summary(diar$n)

#---------------------------------------
# temperature monthly means w/ robust 95% SEs
#---------------------------------------
mys <- unique(weather$my)
temp <- sapply(mys,
                 function(x) {
                   washb_mean(weather$avetemp[weather$my==x],
                              id=c(1:nrow(weather[weather$my==x,])),
                              print=F
                   )
                 }
)
temp <- data.frame(t(temp))
temp$month <- as.Date(mys)
temp$year<-floor_date(temp$month,"year")

colnames(temp) <- c("n","mean","sd","se","lb","ub","month","year")
temp <- temp[order(temp$month),c("month","year","n","mean","sd","se","lb","ub")]

 
# summarize dist of obs per month
summary(temp$n)

  xtics <- unique(temp$month) #Save every month 
  xticsyr <- unique(temp$year) 


#---------------------------------------
# rainfall weekly means w/ robust 95% SEs
#---------------------------------------
mys <- unique(weather$wy)
rain <- sapply(mys,
                 function(x) {
                   washb_mean(weather$rain[weather$wy==x],
                              id=c(1:nrow(weather[weather$wy==x,])),
                              print=F
                   )
                 }
)
rain <- data.frame(t(rain))
rain$month <- as.Date(mys)
rain$year<-floor_date(rain$month,"year")
colnames(rain) <- c("n","mean","sd","se","lb","ub","month","year")
rain <- rain[order(rain$month),c("month","year","n","mean","sd","se","lb","ub")]

 
# summarize dist of obs per week
summary(rain$n)


#---------------------------------------
# make plot w/ color + shaded 95% CIs
#---------------------------------------
pdf("C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/DiarWeatherByMonth-ci.pdf",width=8.5,height=3.541666)

# general plotting parameters

#Set up 1x2 plot
op <- par(mfrow=c(1,2),mar=c(5,6,2,4)+0.1)
#op <- par(mfrow=c(1,2))
xtics2 <- as.Date(c("2008-06-01","2009-02-01"))

cblue <- "#3366AA"
corange <- "#EEA722"
cgrey<-"black"
cols <- c(cgrey,"#56B4E9","#D55E00")

ytics <- seq(0, 0.05, by=0.01)

# empty plot
plot(temp$month,
     c(0,diar$mean),
     type="n",
     ylim=range(ytics),ylab="Diarrhea prevalence (%)",yaxt="n",
     xlim=range(xtics),xlab="Date",xaxt="n",
     bty="n",las=1, cex.lab=0.75)#,
     #mar=c(2,3,2,5)+0.1)
  axis(1,at=xtics, las=2,labels=format(xtics,"%m"),cex.axis=0.6)
  mtext(format(xtics2,"%Y"),at=xtics2,side=1,line=1.5,cex=0.75)
  axis(2,at=ytics,las=1,labels=sprintf("%1.0f",ytics*100))
    legend("topleft",c("Diarrhea","Rainfall"),ncol=1,lty=c(1,1),col=cols[1:2],bty="n",cex=0.7,lwd=1.5)
   
# diarrhea shaded 95% CIs
  polygon(x=c(diar$month,
              rev(diar$month)),
          y=c(diar$lb,
              rev(diar$ub)),
          col=alpha(cols[1],alpha=0.2),border=NA)
  
  
  # diarrhea lines + points
lines(diar$month,diar$mean,lty=1,col=cols[1])
points(diar$month,diar$mean,col=cols[1],pch=19,cex=0.5)

  
  #overlaid empty plot
  par(new=T)
  ytics <- seq(0, 60, by=10)
    #rainfall barplot
barplot(rain$mean, main=NA, 
  	xlab=NA, col=cols[2], ylab=NA,
  	border=NA, axes=F)
  axis(4,at=ytics,las=1,labels=sprintf("%1.0f",ytics))
  mtext(side = 4, line = 2, cex=0.75, 'Mean weekly rainfall (mm)')
  mtext("A)", side=2, line=4.25, adj = 1, las=1, font=2, at=60)


  ytics <- seq(0, 0.05, by=0.01)
# empty plot
plot(temp$month,
     c(0,diar$mean),
     type="n",
     ylim=range(ytics),ylab="Diarrhea prevalence (%)",yaxt="n",
     xlim=range(xtics),xlab="Date",xaxt="n",
     bty="n",las=1, cex.lab=0.75)#,
     #mar=c(2,4,2,6)+0.1)
  axis(1,at=xtics,las=2,labels=format(xtics,"%m"),cex.axis=0.6)
  mtext(format(xtics2,"%Y"),at=xtics2,side=1,line=1.5,cex=0.75)
  axis(2,at=ytics,las=1,labels=sprintf("%1.0f",ytics*100))
   legend("topright",c("Diarrhea","Temperature"),ncol=1,lty=c(1,1),col=cols[c(1,3)],bty="n",cex=0.7,lwd=1.5)
  
# diarrhea shaded 95% CIs
  polygon(x=c(diar$month,
              rev(diar$month)),
          y=c(diar$lb,
              rev(diar$ub)),
          col=alpha(cols[1],alpha=0.2),border=NA)
  
  
  # diarrhea lines + points
lines(diar$month,diar$mean,lty=1,col=cols[1])
points(diar$month,diar$mean,col=cols[1],pch=19,cex=0.5)

  par(new=T)
  ytics <- seq(25, 35, by=5)
  plot(temp$month,
     temp$mean,
     type="n",
     ylim=range(ytics),ylab=NA,yaxt="n",
     xlim=range(xtics),xlab=NA,xaxt="n",
     bty="n",las=1)
  
  # temperature shaded 95% CIs
  polygon(x=c(temp$month,
              rev(temp$month)),
          y=c(temp$lb,
              rev(temp$ub)),
          col=alpha(cols[3],alpha=0.2),border=NA)
  axis(4,at=ytics,las=1,labels=sprintf("%1.0f",ytics))
  mtext(side = 4, line = 2, cex=0.75,  'Mean monthly temperature (C)')
    mtext("B)", side=2, line=4.25, adj = 1, las=1, font=2, at=35)

  # temperature lines + points
lines(temp$month,temp$mean,lty=1,col=cols[3])
points(temp$month,temp$mean,col=cols[3],pch=19,cex=0.5)
  
    
dev.off()











