rm( list = ls()) #clear env

library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(cowplot)
library(DroughtData)
library(sf)
library(rgeos)
library(spacetools)
library(terra)
library(gdistance)


stratum<-readRDS("outputs/stratum.rds")
Zoop_BPUE<-readRDS("data/Zoop_BPUE.rds")
StationLookUp <- read_excel("data/StationLookUp.xlsx")
station_distances<-readRDS("data/station_distances_df.rds")

project_crs<-st_crs(spacetools::Delta,asText=TRUE) #create project CRS from WW_Delta from deltamapr package

stations<-StationLookUp%>%filter(Current==1)
stations<-dplyr::select(stations,StationNZ,latdec,longdec)
stations<-rename(stations,Station=StationNZ)
stations <- stations[complete.cases(stations), ]
stations$longdec<--(stations$longdec)

GG<-data.table::data.table(
  Station="GG",
  latitude=-122.477,
  longitude=37.819539)
GG<-st_as_sf(GG,coords=c("latitude","longitude"),crs=project_crs)
GG <- st_transform(GG, CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs "))
NZ086<-stations%>%filter(Station=="NZ092")
NZ086<-st_as_sf(NZ086,coords=c("longdec","latdec"),crs=project_crs)
NZ086 <- st_transform(NZ086, CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs "))
GG<-st_coordinates(GG)
NZ086<-st_coordinates(NZ086)

load(file="data/Delta_transitioned.rda")
sp<-gdistance::shortestPath(Delta_transitioned,GG,NZ086,output="SpatialLines")
sp<-sp::spTransform(sp, CRS=CRS("+proj=longlat +datum=NAD83 +no_defs"))
sp<-SpatialLinesDataFrame(sp, data = data.frame(ID = 1))

stations<-sf::st_as_sf(stations,coords = c("longdec", "latdec"),crs=project_crs)
stations<-stations%>%left_join(station_distances)

rosies_table<-read.csv("Data/Rosies_regions.csv",fileEncoding="UTF-8-BOM")
stratum<-dplyr::select(stratum,SubRegion,geometry)%>%inner_join(rosies_table)
stratum<-stratum%>%filter(Region!='North')
regions<-unique(stratum$Region)

p<-ggplot()+
  geom_sf(data=stratum,aes(fill=Region))+
  geom_sf(data=spacetools::Delta)+
  geom_line(mapping=aes(x=long,y=lat),data=sp,size=1.5,linetype='twodash')+
  geom_sf(data=stations$geometry,size=5,pch=21,fill='red')+
  #geom_sf_label(data=stations$geometry,nudge_y = .015,label=round(stations$distance_km,digits=1))+
  #geom_sf_label(data=stations$geometry,label=stations$Station)+
  coord_sf(xlim=c(-122.2, -121.38),ylim=c(37.9, 38.25))+
  theme_bw()+
  theme(legend.position = "none",axis.title = element_blank(),axis.text=element_text(size=15)) 
p
save_plot("figures/distance_map_v2.png",p,base_height = 6,base_width = 10)

