##################################
# Creating Inset Maps in ggplot2 #
##################################

rm(list=ls())
library(tidyverse)
library(raster)
library(grid)
library(gridExtra)
library(ggrepel)
library(maps)
library(FField)
library(rasterVis)
library(sp)
library(rgdal)
library(maptools)  ## For wrld_simpl
library(foreign)



#village locations
setwd("C:/Users/andre/Dropbox/Trichy analysis/Trichy data")
gps<-read.csv("VillageGPS.csv", header=T)

#Find mediod gps location for each village
gps <- gps %>% group_by(vilname) %>% summarise(lat=median(north, na.rm=T), long=median(east, na.rm=T)) %>% as.data.frame()
gps$group<-1
head(gps)

#Add in airport
airport <- data.frame(vilname="airport",lat=10.7605952,	long=78.7088753, group=1, season="January - June")

#Get country shapefiles
india<-getData("GADM", country="IND", level=0) # download India level 0 map for ucdavis site
tamnadu<-getData("GADM", country="IND", level=2) # download India level 2 map for ucdavis site
clipregion<-(tamnadu[tamnadu$NAME_2=="Tiruchirappalli" | tamnadu$NAME_2=="Karur" | tamnadu$NAME_2=="Namakkal" | tamnadu$NAME_2=="Perambalur",]) # subset province of Marinduque from PHL map
trichy<-(tamnadu[tamnadu$NAME_1=="Tamil Nadu",]) # subset province of Marinduque from PHL map
munnames<-coordinates(trichy) # get center coordinates of municipalities of Marinduque
munnames<-data.frame(munnames) # convert matrix format munnames object to data.frame
munnames$label<-trichy@data$NAME_2

# Extent rectangle for inset map
pol<-data.frame(xmin=75,xmax=85 ,ymin=13 ,ymax=20)



#Merge in rainfall averages
setwd("C:/Users/andre/Documents/Trichy analysis/Raw Data")
rain<-read.dta("Trichy_Weather_Dec07-Apr09_formatted.dta")
rain<-rain %>% filter(year==2008) %>% subset(., select=c(month, rain)) %>% group_by(month) %>% summarize(rain=sum(rain))
dryrain<-sum(rain[1:6,2])
wetrain<-sum(rain[7:12,2])


setwd("C:/Users/andre/Documents/Trichy analysis/Raw Data/WorldClim precip rasters")
xlim=c(70,85) 
ylim=c(6,14)
worldclim <- list()
worldclim$m1 <- raster("wc2.0_30s_prec_01.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m2 <- raster("wc2.0_30s_prec_02.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m3 <- raster("wc2.0_30s_prec_03.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m4 <- raster("wc2.0_30s_prec_04.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m5 <- raster("wc2.0_30s_prec_05.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m6 <- raster("wc2.0_30s_prec_06.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m7 <- raster("wc2.0_30s_prec_07.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m8 <- raster("wc2.0_30s_prec_08.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m9 <- raster("wc2.0_30s_prec_09.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m10 <- raster("wc2.0_30s_prec_10.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m11 <- raster("wc2.0_30s_prec_11.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])
worldclim$m12 <- raster("wc2.0_30s_prec_12.tif", xmn=xlim[1], xmx=xlim[2], ymn=ylim[1], ymx=ylim[2])

##India polygon
data(wrld_simpl)
SPDF <- subset(wrld_simpl, NAME=="India")

#Function to convert raster to dataframe
raster.to.df<- function(r, poly, xlim=c(70,85), ylim=c(6,14)){
  
## crop and mask
r <- crop(r, extent(poly))

r.spdf <- as(r, "SpatialPixelsDataFrame")
r.df <- as.data.frame(r.spdf)
colnames(r.df) <- c("value","long","lat")
r.df$group<-1

r.df <- r.df %>% 
  filter(long > xlim[1] & long < xlim[2]) %>% 
  filter(lat > ylim[1] & lat < ylim[2])
return(r.df)
}

#clip rasters
worldclim <- lapply(worldclim, function(x) raster.to.df(r=x, poly=clipregion))



dry <- worldclim$m1
dry$value <- dry$value + worldclim$m2$value + worldclim$m3$value + worldclim$m4$value + worldclim$m5$value + worldclim$m6$value

wet <- worldclim$m7
wet$value <- wet$value + worldclim$m8$value + worldclim$m9$value + worldclim$m10$value + worldclim$m11$value + worldclim$m12$value

dry$season <- "January - June"
wet$season <- "July - December"
plotdf<-rbind(dry, wet)

studydry <- studywet <- gps
studydry$season <- "January - June"
studywet$season <- "July - December"
studydry$studyrain <- dryrain
studywet$studyrain <- wetrain
studygps<-rbind(studydry, studywet)
studygps <- studygps %>% subset(., select=c(long, lat, group,season, studyrain)) %>% mutate(value=NA)

df <- bind_rows(cbind(plotdf,geom="raster"), cbind(studygps,geom="points"))





# Main Map
xrange=c(78,79)
yrange=c(10.5,11.5)
p1<-ggplot() + geom_raster(aes(x=long, y=lat,fill=value),data=df[df$geom=="raster",]) +
geom_polygon(data=trichy, aes(long,lat, group=group), colour="grey40", fill=NA)+
geom_point(aes(x = long, y = lat, fill=studyrain, group=group), colour="black",pch=21, data = df[df$geom=="points",], size=5)   + 
geom_text(data=munnames, aes(x=X1+0.02, y=X2+0.08,label=label), size=3, colour="grey20")+
coord_equal()+theme_bw()+xlab("")+ylab("")+
theme(axis.text =element_blank(), strip.background= element_blank(), strip.text.x = element_text(size = 16),axis.ticks = element_blank()) +
coord_cartesian(xlim=xrange, ylim=yrange)  


p1 <- p1 + scale_fill_gradient(limits = c(100,1100), low = "#EBEB86", high = "#0048E3", guide_colourbar(title = "Total seasonal\nrainfall (mm)")) 

#add in airport
p1 <- p1 + geom_point(aes(x = long, y = lat),pch=13,  data = airport, size=5)  + facet_wrap(~season)
p1 <- p1 + geom_text(aes(x = long, y = lat - 0.08, label = "Tiruchchirappalli\nAirport"), colour="grey20", size=3, data = airport, fontface='bold',inherit.aes=FALSE)



#Add India map inset
p2<-ggplot()+geom_polygon(data=india, aes(long,lat,group=group),colour="grey10",fill="#fff7bc")+  
coord_equal()+theme_bw()+labs(x=NULL,y=NULL)+
geom_rect(data = pol, aes(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2]), alpha=0, colour="red", size = 0.5, linetype=1)+
theme(axis.text.x =element_blank(),axis.text.y= element_blank(), axis.ticks=element_blank(),axis.title.x =element_blank(),
axis.title.y= element_blank())

#Save map
setwd("C:/Users/andre/Documents/Trichy analysis/Figures and Tables/")
png(file="rainfallmap.png",w=2000,h=1200, res=300)
#tiff(file = "C:/Users/andre/Documents/Trichy analysis/Figures and Tables/rainfallmap.tiff",w=2000,h=1200, units = "px", res=300) 

grid.newpage()
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.3, height = 0.3, x = 0.7127, y = 0.19425) #plot area for the inset map
print(p1,vp=v1) 
print(p2,vp=v2)
dev.off()






