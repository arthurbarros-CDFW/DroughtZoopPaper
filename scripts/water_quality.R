library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(cowplot)
library(viridis)
library(DroughtData)

drought_wq <- read.csv("data/Drought_Year_Temp_Secchi_Salinity.csv")

WY<-read_excel("Data/Water years.xlsx", sheet="yearassignments")

drought_wq<-drought_wq%>%inner_join(WY)
drought_wq<-drought_wq%>%filter(Region!="North")
drought_wq$Region<-factor(drought_wq$Region,levels=c("Suisun Bay","Suisun Marsh","Confluence","SouthCentral"))
regions<-unique(drought_wq$Region)

for(i in 1:length(regions)){
  r<-regions[[i]]
  d<-drought_wq%>%filter(Region==r)
  pm1<-ggplot(d, aes(x = Drought, y = Salinity, fill=Drought)) +
    geom_boxplot(aes())+
    theme_bw()+
    drt_color_pal_drought(aes_type = "fill")+
    ggtitle(r)+
    theme(text = element_text(size=24))+
    xlab("Year Type")+ylab("Avg Salinity (ppt)")
  pm1
  save_plot(paste("figures/WQ/salinity",r,".png"),pm1,base_height = 8,base_width = 12)
}

m1<-aov(Salinity~Drought*Region,data=drought_wq)
summary(m1)
tm1<-TukeyHSD(m1)
tm1

for(i in 1:length(regions)){
  r<-regions[[i]]
  d<-drought_wq%>%filter(Region==r)
  pm2<-ggplot(d, aes(x = Drought, y = Temperature, fill=Drought)) +
    geom_boxplot(aes())+
    theme_bw()+
    drt_color_pal_drought(aes_type = "fill")+
    ggtitle(r)+
    theme(text = element_text(size=24))+
    xlab("Year Type")+ylab("Avg temperature (c)")
  pm2
  save_plot(paste("figures/WQ/temperature",r,".png"),pm2,base_height = 8,base_width = 12)
}


m2<-aov(Temperature~Drought*Region,data=drought_wq)
summary(m2)
tm2<-TukeyHSD(m2)
tm2

for(i in 1:length(regions)){
  r<-regions[[i]]
  d<-drought_wq%>%filter(Region==r)
  pm3<-ggplot(d, aes(x = Drought, y = Secchi, fill=Drought)) +
    geom_boxplot(aes())+
    theme_bw()+
    drt_color_pal_drought(aes_type = "fill")+
    ggtitle(r)+
    theme(text = element_text(size=24))+
    xlab("Year Type")+ylab("Avg Secchi (cm)")
  pm3
  save_plot(paste("figures/WQ/secchi",r,".png"),pm3,base_height = 8,base_width = 12)
}


m3<-aov(Secchi~Drought*Region,data=drought_wq)
summary(m3)
tm3<-TukeyHSD(m3)
tm3
