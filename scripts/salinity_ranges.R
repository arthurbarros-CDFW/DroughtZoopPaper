#Taxa salinity ranges
#GAMS and GLMS of fecundity data
#trying to look at time changes in mysid fecundity?
rm( list = ls()) #clear env

library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(cowplot)
library(mgcv)
library(splines)
library(visreg)
library(pscl)
library(glmmTMB)
library(MASS)
library(DHARMa)
library(ggeffects)
library(brms)
library(DroughtData)
library(rstan)
library(Hmisc)


Zoop_BPUE<-readRDS("data/Zoop_BPUE.rds")
ZPEnvData<-read_excel("data/ZPEnvData.xlsx")
StationLookUp <- read_excel("data/StationLookUp.xlsx")

#filter Zoop_BPUE for EMP
EMP_BPUE<-Zoop_BPUE%>%filter(Source=="EMP")

#set dates, years, months
EMP_BPUE$SampleDate<-as.Date(EMP_BPUE$Date,"%m/%d/%Y")
EMP_BPUE$month<-as.numeric(format(as.Date(EMP_BPUE$SampleDate), "%m"))
EMP_BPUE$year<-as.numeric(format(as.Date(EMP_BPUE$SampleDate), "%Y"))

winter<-c("12","01","02")
EMP_BPUE<-EMP_BPUE%>%filter(!month %in%winter)
EMP_BPUE<-EMP_BPUE%>%
  group_by(Station,SampleDate,year,month,Taxlifestage)%>%
  dplyr::summarise(BPUE=sum(BPUE))

#set stations for env data
ZPEnvData_stations<-ZPEnvData%>%inner_join(StationLookUp)
ZPEnvData_stations$SampleDate<-as.Date(ZPEnvData_stations$SampleDate,"%m/%d/%Y")
ZPEnvData_stations$Station<-ZPEnvData_stations$StationNZ
ZPEnvData_stations<-ZPEnvData_stations%>%filter(Core==1)
ZPEnvData<-dplyr::select(ZPEnvData_stations,SampleDate,Station,secchi,beg_surf_temp,beg_surf_sc)

#filter for target taxa
target_DROUGHT_taxa<-c("Daphnia Adult",
                       "Pseudodiaptomus forbesi Adult",
                       "Limnoithona tetraspina Adult",
                       "Hyperacanthomysis longirostris Adult"
)

EMP_BPUE<-EMP_BPUE%>%
  filter(Taxlifestage%in%target_DROUGHT_taxa & year>1993)
EMP_BPUE<-dplyr::select(EMP_BPUE,SampleDate,Station,Taxlifestage,BPUE)
EMP_BPUE<-unique(EMP_BPUE)

#join env and bpue data
Env_BPUE<-EMP_BPUE%>%inner_join(ZPEnvData)
Env_BPUE<-Env_BPUE%>%
  filter(!is.na(BPUE))
Env_BPUE<-unique(Env_BPUE)

Env_BPUE_wide<-Env_BPUE%>%pivot_wider(names_from = Taxlifestage,values_from=BPUE,values_fill=0)
Env_BPUE<-Env_BPUE_wide%>%pivot_longer(cols=8:11,names_to="Taxlifestage",values_to="BPUE")

WY<-read_excel("Data/Water years.xlsx", sheet="yearassignments")
WY$year<-WY$Year
Env_BPUE<-Env_BPUE%>%inner_join(WY)

#Make bar plots of BPUE by salinity bins
taxa<-(unique(Env_BPUE$Taxlifestage))
target_months<-c(5,6,7,8,9,10,11)
model_data<-Env_BPUE
source("scripts/ec_to_sal.R")
model_data$salinity<-ec_2_sal(25,model_data$beg_surf_sc)
bin_size<-2
for(i in 1:length(taxa)){
  t<-taxa[[i]]
  d<-model_data%>%filter(Taxlifestage==t & month)
  d$bin_sal<-cut(d$salinity,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
  d<-d%>%filter(!is.na(bin_sal))
  d<-d%>%
    group_by(year,bin_sal)%>%
    dplyr::summarise(BPUE=mean(BPUE))
  p<-ggplot(d, aes(x = bin_sal, y = log(BPUE+1))) +
    geom_boxplot(aes())+
    theme_bw()+
    drt_color_pal_drought(aes_type = "fill")+
    ggtitle(t)+
    theme(text = element_text(size=24),axis.text = element_text(size=12))+
    xlab("Salinity (ppt) bin")+ylab("ln(BPUE+1)")
  p
  save_plot(paste("figures/sal_zones/bins",t,".png"),p,base_height = 8,base_width = 12)
}

###########
#try out Rosies suggested weighted means
plot_list<-list()
xm<-c(0,0,0,0)
sd<-c(0,0,0,0)
d2<-data.frame(taxa,xm,sd)
for(i in 1:length(taxa)){
  t<-taxa[[i]]
  d<-model_data%>%filter(Taxlifestage==t)
  d2[i,2]<-wtd.mean(d$salinity,d$BPUE)
  var<-wtd.var(d$salinity,d$BPUE)
  d2[i,3]<-sqrt(var)
}
write.csv(d2,"outputs/salinity_ranges.csv",row.names = F)
p<-ggplot(d2,aes(x=taxa,y=xm))+
  geom_point(stat="identity",color='black',size=2)+
  geom_errorbar(aes(ymin=xm-sd,ymax=xm+sd,width=.2))+
  coord_flip()+
  labs(y="Mean salinity",x="")+
  theme_light()
p
save_plot("figures/sal_zones/bins.png",p,base_height = 4,base_width = 6)

#Okay, examining the salinity ranges these are the prefered ranges based on weighted
#means and standard deviance:
#Daphnia: 0 - 
#H longirostris: 1 - 7ppt
#L tetraspina 1 - 15ppt
#P forbesi <5ppt
sal_min<-c(0.0,3.2,1.6,0.0)
sal_max<-c(0.5,6.7,9.1,1.8)
Taxlifestage<-taxa
taxa_sal_range<-data.frame(Taxlifestage,sal_min,sal_max)
#so next for each taxa we will do a boxplot and TukeyHSD for drought year type within target salinity range

for(i in 1:length(taxa)){
  t<-taxa[[i]]
  sal<-taxa_sal_range%>%filter(Taxlifestage==t)
  d<-model_data%>%filter(Taxlifestage==t)%>%
    filter(between(salinity, sal$sal_min,sal$sal_max))
  d<-d%>%
    filter(Drought!="N")%>%
    group_by(month,year,Drought)%>%
    dplyr::summarise(BPUE=mean(BPUE))
  p<-ggplot(d, aes(x = Drought, y = log(BPUE+1), fill=Drought)) +
    geom_boxplot(aes())+
    theme_bw()+
    drt_color_pal_drought(aes_type = "fill")+
    ggtitle(t)+
    theme(text = element_text(size=24),legend.position = "none")+
    xlab("Year Type")+ylab("ln(BPUE+1)")
  p
  save_plot(paste("figures/sal_zones/",t,".png"),p,base_height = 8,base_width = 8)
  m<-aov(BPUE~Drought,data=d)
  summary(m)
  TukeyHSD(m)
  capture.output(summary(m),file=paste("figures/sal_zones/",t,".txt"))
  capture.output(TukeyHSD(m),file=paste("figures/sal_zones/",t,"TUKEY.txt"))
}

############Next lets plot where the salinity zone was?
ez_distances<-readRDS("data/ez_distances.rds")
ez_distances<-ez_distances%>%
  rename(SampleDate=Date, distance_km=ez_km)
ez_distances$SampleDate<-as.Date(ez_distances$SampleDate,"%Y-%m-%d")
ez_distances<-dplyr::select(ez_distances,SampleDate,Station,distance_km)
station_distances<-readRDS("data/station_distances_df.rds")
station_distances<-dplyr::select(station_distances,Station,distance_km)

sal_ez_distances<-model_data%>%
  inner_join(ez_distances)
sal_distances<-model_data%>%
  inner_join(station_distances)

sal_distances<-sal_distances%>%rbind(sal_ez_distances)

#plot salinity zone distribution for each species and drought/wet year
for(i in 1:length(taxa)){
  t<-taxa[[i]]
  sal<-taxa_sal_range%>%filter(Taxlifestage==t)
  d<-sal_distances%>%filter(Taxlifestage==t)%>%
    filter(between(salinity, sal$sal_min,sal$sal_max) & Drought!="N")
  d<-unique(dplyr::select(d,SampleDate,Station,salinity,Drought,distance_km))
  p<-ggplot(d, aes(x = Drought, y = distance_km, fill=Drought)) +
    annotate("rect",ymin=30,ymax=42,xmin=-Inf,xmax=Inf,alpha=.3,fill="blue")+
    annotate("rect",ymin=42,ymax=56,xmin=-Inf,xmax=Inf,alpha=.3,fill="green")+
    annotate("rect",ymin=56,ymax=75,xmin=-Inf,xmax=Inf,alpha=.3,fill="yellow")+
    annotate("rect",ymin=75,ymax=90,xmin=-Inf,xmax=Inf,alpha=.3,fill="orange")+
    annotate("rect",ymin=90,ymax=105,xmin=-Inf,xmax=Inf,alpha=.3,fill="red")+
    annotate("rect",ymin=105,ymax=140,xmin=-Inf,xmax=Inf,alpha=.3,fill="purple")+
    annotate("text", y = 36, x = "D", label = "San Pablo \n Bay",size=3,vjust=3.5)+
    annotate("text", y = 49, x = "D", label = "Carquinez \n Strait",size=3,vjust=3.5)+
    annotate("text", y = 66, x = "D", label = "Suisun",size=3,vjust=8)+
    annotate("text", y = 82, x = "D", label = "West \n Delta",size=3,vjust=3.5)+
    annotate("text", y = 98, x = "D", label = "Central Delta",size=3,vjust=8)+
    annotate("text", y = 122, x = "D", label = "East Delta",size=3,vjust=8)+
    geom_boxplot(aes())+
    theme_bw()+
    drt_color_pal_drought(aes_type = "fill")+
    ggtitle(paste(t," target salinity zone"),subtitle =paste(sal$sal_min," ppt - ",sal$sal_max," ppt",sep=""))+
    theme(text = element_text(size=16),legend.position = "none")+
    xlab("Year Type")+ylab("River km")+
    coord_flip()
  p
  save_plot(paste("figures/sal_zones/",t,"distances.png"),p,base_height = 5,base_width = 8)
}
