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

Zoop_BPUE<-readRDS("data/Zoop_BPUE.rds")
ZPEnvData<-read_excel("data/ZPEnvData.xlsx")
StationLookUp <- read_excel("data/StationLookUp.xlsx")

#filter Zoop_BPUE for EMP
EMP_BPUE<-Zoop_BPUE%>%filter(Source=="EMP")

#set dates, years, months
#EMP_BPUE$SampleDate<-as.Date(EMP_BPUE$Date,"%m/%d/%Y")
#It looks like this read in the dates as dates already, so you don't need this step

EMP_BPUE$SampleDate<-EMP_BPUE$Date
EMP_BPUE$month<-as.numeric(format(as.Date(EMP_BPUE$SampleDate), "%m"))
EMP_BPUE$year<-as.numeric(format(as.Date(EMP_BPUE$SampleDate), "%Y"))
#you can use the 'lubridate' package, which has the 'month' and 'year' functions. I can show you sometime if you like.


#need to filter for target taxa within target gear/size class
Taxlifestage<-c("Bosmina longirostris Adult",
                                    "Pseudodiaptomus forbesi Adult",
                                    "Limnoithona tetraspina Adult",
                                    "Hyperacanthomysis longirostris Adult")
SizeClass<-c("Meso","Meso","Micro","Macro")
target_sizes<-data.frame(Taxlifestage,SizeClass)
EMP_BPUE<-EMP_BPUE%>%
  inner_join(target_sizes)%>%
  filter(year>1993)
EMP_BPUE<-EMP_BPUE%>%
  ungroup()
EMP_BPUE<-dplyr::select(EMP_BPUE,-SizeClass)%>%ungroup()

#set stations for env data
ZPEnvData_stations<-ZPEnvData%>%inner_join(StationLookUp)

#this got read in correctly, so not needed
#ZPEnvData_stations$SampleDate<-as.Date(ZPEnvData_stations$SampleDate,"%m/%d/%Y")

ZPEnvData_stations$Station<-ZPEnvData_stations$StationNZ
ZPEnvData_stations<-ZPEnvData_stations%>%filter(Core==1)
ZPEnvData<-dplyr::select(ZPEnvData_stations,SampleDate,Station,secchi,beg_surf_temp,beg_surf_sc)

#join env and bpue data
Env_BPUE<-EMP_BPUE%>%inner_join(ZPEnvData)
Env_BPUE<-Env_BPUE%>%
  filter(!is.na(BPUE))

#pivot and fill in 0s
Env_BPUE_wide<-Env_BPUE%>%pivot_wider(names_from = Taxlifestage,values_from=BPUE,values_fill=0)
Env_BPUE<-Env_BPUE_wide%>%pivot_longer(cols=14:17,names_to="Taxlifestage",values_to="BPUE")

#bring in WY drought data
WY<-read_excel("Data/Water years.xlsx", sheet="yearassignments")
WY$year<-WY$Year
Env_BPUE<-Env_BPUE%>%inner_join(WY)

target_DROUGHT_taxa<-c("Bosmina longirostris Adult",
                       "Pseudodiaptomus forbesi Adult",
                       "Limnoithona tetraspina Adult",
                       "Hyperacanthomysis longirostris Adult")
Env_BPUE<-Env_BPUE%>%filter(Taxlifestage %in% target_DROUGHT_taxa)

##################
#Hists of BPUE
#################
taxa<-(unique(Env_BPUE$Taxlifestage))
for(i in 1:length(taxa)){
  t<-taxa[i]
  d<-Env_BPUE%>%filter(Taxlifestage==t)
  p1<-ggplot(d,aes(x=BPUE))+
    geom_histogram(color="black", fill="white")+
    ggtitle(paste("BPUE",t))
  p1
  p2<-ggplot(d,aes(x=log(BPUE+.01)))+
    geom_histogram(color="black", fill="white")+
    ggtitle(paste("logBPUE",t," 1993-2020"))
  p2
  ggsave(plot=p1,paste("figures/histograms/",t,".png"),height=4,width=6)
  ggsave(plot=p2,paste("figures/histograms/",t,"_log.png"),height=4,width=6)
}
#looking at histograms all the taxa are zero-inflated as expected

##################
#Gams
##################
target_months<-c(5,6,7,8,9,10,11) #looking at May-Nov, when zooplankton populations are highest
model_data<-Env_BPUE

#this is a nifty script! You can also try using the ec2pss function from the 'wql' package
source("scripts/ec_to_sal.R")
model_data$salinity<-ec_2_sal(25,model_data$beg_surf_sc)

#scale temperature by monthly mean
#temp_monthlymean<-model_data%>%   #removing for now
  #group_by(year,month)%>%
  #dplyr::summarise(month_temp=mean(beg_surf_temp,na.rm=T))
#model_data<-model_data%>%
  inner_join(temp_monthlymean)
#model_data$scaled_temp<-model_data$beg_surf_temp-model_data$month_temp

for(i in 1:length(taxa)){
  t<-taxa[i]
  d<-model_data%>%filter(Taxlifestage==t & !is.na(salinity) & year>1993 & month %in% target_months)
  d<-dplyr::select(d,BPUE,salinity,month,Station,Drought)
  d$Station<-as.factor(d$Station)
  d$Drought<-as.factor(d$Drought)
  
  m1<-gam(BPUE~s(salinity),family="nb",data=d)
  summary(m1)
  m2<-gam(BPUE~s(salinity)+s(month,k=5),family='nb',data=d)
  summary(m2)
  #Now with a random effect of station
  m3<-gamm(BPUE~s(salinity)+s(month,k=5),random = list(Station = ~1), niterPQL=40,family='nb',data=d)
  summary(m3[[1]])
  summary(m3[[2]])
  
  #try two seperate models for p/a binomial and just presence nb
  d_pa<-d
  d_pa$Presence<-ifelse(d_pa$BPUE>0,1,0)
  d_p<-d_pa%>%filter(Presence==1)
  m4.1<-glm(Presence~salinity+month,family="binomial",data=d_pa)
  summary(m4.1)
  m4.2<-gam(BPUE~s(salinity)+s(month,k=5),random = list(Station = ~1), niterPQL=40,family='nb',data=d_p)
  summary(m4.2)
 
  capture.output(summary(m1),file=paste("outputs/model_outputs/",t,"_m1.txt"))
  capture.output(summary(m2),file=paste("outputs/model_outputs/",t,"_m2.txt"))
  capture.output(summary(m3[[2]]),file=paste("outputs/model_outputs/",t,"_m3.txt"))
  capture.output(summary(m4.1),file=paste("outputs/model_outputs/",t,"_m4_1.txt"))
  capture.output(summary(m4.2),file=paste("outputs/model_outputs/",t,"_m4_2.txt"))
  
  #Anova or AIC to test compare model fits?
  m_anova<-anova(m1,m2,m3, m4.1)
  m_anova
  capture.output(m_anova,file=paste("outputs/model_outputs/",t,"_m_anova.txt"))
  m_aic<-AIC(m1,m2,m4.1,m4.2)
  m_aic
  capture.output(m_aic,file=paste("outputs/model_outputs/",t,"_m_aic.txt"))
  
  png(paste("figures/gamcheck/",t,"_m4_1.png"))
  par(mfrow=c(2,2))
  gam.check(m4.1)
  dev.off()
  
  png(paste("figures/gamcheck/",t,"_m4_2.png"))
  par(mfrow=c(2,2))
  gam.check(m4.2)
  dev.off()
  
}

#another method for random effects:
  m5 = gamm(BPUE~s(salinity)+s(month,k=5) + s(Station, bs = "re"), niterPQL=40,family='nb',data=d)
  summary(m5[[1]])
  summary(m5[[2]])
  plot(m5[[2]])
  