rm(list = ls()) #clear env
#load packages
library(tidyverse)
library(tibble)
library(dplyr)
library(ggplot2)
library(plotrix)g
library(purrr)
library(cowplot)
library(plyr)
library(readxl)
library(deltamapr)
library(spacetools)
library(sf)

#################################
#DROUGHT Center of distribution metrics
#################################

dist_center<-function(data,taxa){
  d=data%>%
    filter(Taxlifestage==taxa)
  d$km_freq<-d$BPUE*d$distance_km
  d1<-d%>%
    group_by(Survey,seasons,Year,Taxlifestage,Regime)%>%
    dplyr::summarize(km_center=sum(km_freq)/sum(BPUE),Mean_X2=mean(X2))
}

