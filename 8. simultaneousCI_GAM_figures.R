rm(list=ls())
library("ggplot2")
library("ggthemes")
theme_set(theme_bw(base_size = 8))


#Set jitter height
h=0.2


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


ulabplot <- function(title) {
  plot(1,1,type="n",
       xaxt="n",xlab="",xlim=c(0,1),
       yaxt="n",ylab="",bty="n",ylim=c(0,1)
  )
  text(1,0.5,title,adj=1,cex=1.5)
}


# grab a color blind friendly palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1,3,7)]

setwd("C:/Users/andre/Dropbox/Trichy analysis/Results/")
load("C:/Users/andre/Dropbox/Trichy analysis/Results/temp.unadjusted.GAMfits.Rdata")
load("C:/Users/andre/Dropbox/Trichy analysis/Results/rain.unadjusted.GAMfits.Rdata")
load("C:/Users/andre/Dropbox/Trichy analysis/Results/rain.unadjusted.GAMfits.stratified.Rdata")
load("C:/Users/andre/Dropbox/Trichy analysis/Results/temp.quartiles.Rdata")
load("C:/Users/andre/Dropbox/Trichy analysis/Results/HRthres.Rdata")




d.temp<-rbind(
data.frame(fit.temp7.unadj, lag="One week lag"),
data.frame(fit.temp14.unadj, lag="Two week lag"),
data.frame(fit.temp21.unadj, lag="Three week lag")
)

#Convert prev to percent
d.temp$fit<-d.temp$fit*100
d.temp$uprS<-d.temp$uprS*100
d.temp$lwrS<-d.temp$lwrS*100

#merge in quartile info
d.temp$quartiles<-rep(NA, nrow(d.temp))
d.temp$quartiles[d.temp$lag=="One week lag"]<-temp_quartiles[,1] 
d.temp$quartiles[d.temp$lag=="Two week lag"]<-temp_quartiles[,2] 
d.temp$quartiles[d.temp$lag=="Three week lag"]<-temp_quartiles[,3] 

temp<-ggplot(d.temp, aes(x = newd)) +
    geom_line(aes(y=fit), color="#D55E00") +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "#D55E00") +
    geom_jitter(aes(y=Y*6,x=X), height = h, width=0.2,  alpha = 0.1, size=0.5) +
    geom_vline(aes(xintercept=quartiles, color = "#D55E00"), alpha = 1, linetype = 2) +
    facet_grid( ~ lag) +
    labs(y = "Diarrhea prevalence (%)",
         x = "Weekly daytime average (C)",
         title = "(A) Temperature") +
    theme(strip.background = element_blank(),
          legend.position="none")




d.rain<-rbind(
data.frame(fit.rain7.unadj, lag="One week lag"),
data.frame(fit.rain14.unadj, lag="Two week lag"),
data.frame(fit.rain21.unadj, lag="Three week lag")
)



#Convert prev to percent
d.rain$fit<-d.rain$fit*100
d.rain$uprS<-d.rain$uprS*100
d.rain$lwrS<-d.rain$lwrS*100
rain<-ggplot(d.rain, aes(x = newd)) +
    geom_line(aes(y=fit), color="#56B4E9") +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "#56B4E9") +
    geom_jitter(aes(y=Y*6,x=X), height = h, width=0.2,  alpha = 0.1, size=0.5)+
    facet_grid( ~ lag) +
    labs(y = "Diarrhea prevalence (%)",
         x = "Weekly 24 hour accumulation average (log mm)",
         title = "(B) Rainfall") +
    theme(strip.background = element_blank())


setwd("C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/")
#pdf("GAMcurves.pdf",width=10,height=8)
pdf("GAMcurves5.2.pdf",width=5.2,height=4.16)

multiplot(temp, rain, cols=1)

dev.off()


#Stratified rainfall

d.rain.strat<-rbind(
data.frame(fit.rain7.unadjLT1, lag="One week lag", strat="Low long-term rainfall"),
data.frame(fit.rain14.unadjLT1, lag="Two week lag", strat="Low long-term rainfall"),
data.frame(fit.rain21.unadjLT1, lag="Three week lag", strat="Low long-term rainfall"),
data.frame(fit.rain7.unadjLT2, lag="One week lag", strat="Medium long-term rainfall"),
data.frame(fit.rain14.unadjLT2, lag="Two week lag", strat="Medium long-term rainfall"),
data.frame(fit.rain21.unadjLT2, lag="Three week lag", strat="Medium long-term rainfall"),
data.frame(fit.rain7.unadjLT3, lag="One week lag", strat="High long-term rainfall"),
data.frame(fit.rain14.unadjLT3, lag="Two week lag", strat="High long-term rainfall"),
data.frame(fit.rain21.unadjLT3, lag="Three week lag", strat="High long-term rainfall")
)

#Convert prev to percent
d.rain.strat$fit<-d.rain.strat$fit*100
d.rain.strat$uprS<-d.rain.strat$uprS*100
d.rain.strat$lwrS<-d.rain.strat$lwrS*100

rainstrat<-ggplot(d.rain.strat, aes(x = newd)) +
    coord_cartesian(ylim = c(0, 9.5)) +
    geom_line(aes(y=fit), color="#56B4E9") +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "#56B4E9") +
    geom_jitter(aes(y=Y*8,x=X), height = h*1.5, width=0.2,  alpha = 0.1, size=0.5)+
    facet_grid( strat~ lag) +
    labs(y = "Diarrhea prevalence (%)",
         x = "Weekly 24 hour Rainfall accumulation average (log mm)") +
    theme(strip.background = element_blank())



setwd("C:/Users/andre/Dropbox/Trichy analysis/Figures and Tables/")
#pdf("GAMcurves_stratified.pdf",width=10,height=8)
pdf("GAMcurves_stratified5.2.pdf",width=5.2,height=4.16)

rainstrat
dev.off()















