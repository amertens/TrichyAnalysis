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
setwd("C:/Users/andre/Dropbox/Trichy analysis/Data")
gps<-read.csv("VillageGPS.csv", header=T)

# Village names Melpathu and Thudaiyur are the same village ID/intervention village with the same coordinates, so combine here
gps$vilname <- as.character(gps$vilname)
gps$vilname[gps$vilname=="Thudaiyur"] <- "Melpathu"

#Change spellings/punctuation to be consistent with survey spellings
gps$vilname[gps$vilname=="Thiruvallarai"] <- "Thiruvellarai"
gps$vilname[gps$vilname=="V A Samudram"] <- "V.A.Samudram"
gps$vilname[gps$vilname=="Ponnusangampatti"] <- "Ponnusangam Patti"
gps$vilname[gps$vilname=="O. Krishnapuram"] <- "Krishnapuram"
gps$vilname[gps$vilname=="Mettuppatti"] <- "Mettupatti"

gps$vilname[gps$vilname=="Ayinapatti"] <- "Ayinapatty"
gps$vilname[gps$vilname=="Devarappampatti"] <- "Thevarappampatti"
gps$vilname[gps$vilname=="E.Batherpettai"] <- "E.Badarpettai"
gps$vilname[gps$vilname=="Keela Karthigaipatti"] <- "Keelakarthigaipatti"
gps$vilname[gps$vilname=="Mela Karthigaipatti"] <- "Melakarthikaipatti"
gps$vilname[gps$vilname=="Mela Naduvalur"] <- "Naduvalur Melur"
gps$vilname[gps$vilname=="Melakothampatti"] <- "Mela Kothampatti"
gps$vilname[gps$vilname=="O. Krishnapuram"] <- "Krishnapuram"


#Find mediod gps location for each village
gps <- gps %>% 
  group_by(vilname) %>% summarise(lat=median(north, na.rm=T), long=median(east, na.rm=T)) %>% as.data.frame()
gps$group<-1
head(gps)


#Merge in treatment allocation
svy<-read.dta("trichy_long.dta") %>%
  subset(select=c("vilid","vilname","wpi")) %>%
  group_by(vilid) %>% slice(1)

df <- merge(gps, svy, by="vilname", all.x=T, all.y=T)
df$`Intervention\nstatus` <- factor(ifelse(df$wpi==1,"Intervention", "Control"))


#Add in airport
airport <- data.frame(vilname="airport",lat=10.7605952,	long=78.7088753, group=1)

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



##India polygon
data(wrld_simpl)
SPDF <- subset(wrld_simpl, NAME=="India")




#greyscale and colorblind friendly colorss
#cols=c("#56B4E9","#D55E00")
cols=c("#f7fcb9","#31a354")

# Main Map
xrange=c(78,79)
yrange=c(10.5,11.5)
p1<-ggplot() +
  geom_polygon(data=trichy, aes(long,lat, group=group), colour="grey40", fill=NA)+
  geom_point(aes(x = long, y = lat, fill=`Intervention\nstatus`, group=`Intervention\nstatus`), colour="grey20", pch=21, data = df, size=3, alpha=0.7)   + 
  scale_fill_manual(values = cols) +
  geom_text(data=munnames, aes(x=X1+0.02, y=X2+0.08,label=label), size=3, colour="grey20")+
  coord_equal()+theme_bw()+xlab("")+ylab("")+
  theme(axis.text =element_blank(), strip.background= element_blank(), strip.text.x = element_text(size = 16),axis.ticks = element_blank(),
        legend.position="right",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed(1, xlim=xrange, ylim=yrange) 



#add in airport
p1 <- p1 + geom_point(aes(x = long, y = lat),pch=13,  data = airport, size=5) 
p1 <- p1 + geom_text(aes(x = long, y = lat + 0.05, label = "Tiruchirappalli\nAirport"), colour="grey20", size=3, data = airport, fontface='bold',inherit.aes=FALSE)

p1 <- p1 + theme(legend.position = c(0.17, 0.8), 
           legend.key = element_rect(size = 0.5),
           legend.title = element_text(size = 10), legend.text = element_text(size = 8),
           legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"),
           legend.background = element_rect(color = "grey10", size = 0.3, linetype = "solid"))

#Add India map inset
p2<-ggplot()+geom_polygon(data=india, aes(long,lat,group=group),colour="grey10",fill="#fff7bc")+  
  coord_equal()+theme_bw()+labs(x=NULL,y=NULL)+
  geom_rect(data = pol, aes(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2]), alpha=0, colour="red", size = 0.75, linetype=1)+
  theme(axis.text.x =element_blank(),axis.text.y= element_blank(), axis.ticks=element_blank(),axis.title.x =element_blank(),
        axis.title.y= element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

#Save map
setwd("C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/")



pdf(file="Figure1.pdf",width=6,height=4.8)

grid.newpage()
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.4, height = 0.4, x = 0.295, y = 0.239, clip="off") #plot area for the inset map
print(p1,vp=v1) 
print(p2,vp=v2)
dev.off()




