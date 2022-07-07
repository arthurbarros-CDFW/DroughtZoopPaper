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
target_DROUGHT_taxa<-c("Bosmina longirostris Adult",
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
Env_BPUE<-Env_BPUE_wide%>%pivot_longer(cols=8:21,names_to="Taxlifestage",values_to="BPUE")

WY<-read_excel("Data/Water years.xlsx", sheet="yearassignments")
WY$year<-WY$Year
Env_BPUE<-Env_BPUE%>%inner_join(WY)

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
target_months<-c(5,6,7,8,9,10,11)
model_data<-Env_BPUE
source("scripts/ec_to_sal.R")
model_data$salinity<-ec_2_sal(25,model_data$beg_surf_sc)
#scale temperature by monthly mean
temp_monthlymean<-model_data%>%
  group_by(year,month)%>%
  dplyr::summarise(month_temp=mean(beg_surf_temp,na.rm=T))
model_data<-model_data%>%
  inner_join(temp_monthlymean)
model_data$scaled_temp<-model_data$beg_surf_temp-model_data$month_temp
#########################################
#single smoother GAMs for salinity to find range of target taxa

for(i in 1:length(taxa)){
  t<-taxa[i]
  d<-model_data%>%filter(Taxlifestage==t)
  d<-d%>%
    filter(!is.na(salinity) & month%in%target_months)
  #run and plot model
  m_sal<-mgcv::gam(BPUE~s(salinity),family="nb",data=d)
  capture.output(summary(m_sal),file=paste("outputs/salinity_gams/",t,"ppt_model.txt"))
  summary(m_sal)
  png(paste("figures/single_gams/",t,"_ppt_smooth_summer.png"))
  par(mfrow=c(1,1))
  plot(m_sal,scale=0,main = paste(t," May-Nov"))
  abline(0,0)
  dev.off()
  
  #evaluate model fit
  png(paste("figures/single_gams/",t,"_ppt_fit_summer.png"))
  layout(mat = matrix(c(1,2,3,3), 2, 2, byrow=TRUE))
  par(mar=c(4.5,4.5,1.5,1.5)+0.1, oma=rep(0,4))
  plot(d$BPUE ~ predict(m_sal), las=1, pch=20, col=gray(0.2,0.2),main = t)
  abline(0,1)
  hist(residuals(m_sal), las=1, main=""); box(which="plot")
  plot(residuals(m_sal) ~ d$salinity, las=1, pch=20, col=gray(0.2,0.2))
  abline(h=0)
  dev.off()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~ Visualize model uncertainty:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Determine the confidence interval to present in your figure and determine the
  # t-distribution value of which to multiply the standard error
  interval_value = 0.95
  t_distribution_probability = 1.0 - (1-interval_value)/2
  CI_factor = qt(t_distribution_probability, Inf)
  
  #~ Generate a data frame with the sequence from the minimum to the maximum
  # observed value of the predictor variable: day of water year
  m_sal_range = range(na.omit(d$salinity))
  new_dat_1 = data.frame(salinity = seq(from = m_sal_range[1],
                                        to = m_sal_range[2], by=1))
  
  #~ Use the *predict()* function to estimate model fit and standard error
  preds = predict(m_sal, newdata = new_dat_1, se.fit = TRUE,type = "response")
  
  # Multiply the standard error fit by the t-distribution value to estimate 
  # the upper and lower confidence interval limits
  fit = preds$fit
  upper = fit-CI_factor*preds$se.fit
  lower = fit+CI_factor*preds$se.fit
  
  
  #~~   Visualize the model
  
  # Plot the raw data and make the axes and labels aesthetically pleasing
  png(paste("figures/single_gams/",t,"_ppt_predict_summer.png"))
  par(mfrow=c(1,1),mar=c(2.5, 4.5, 2,4.5)+0.1, oma=c(1.5,0,0.5,0))
  plot(d$BPUE ~ d$salinity, pch=20, col=gray(0.1,0.2),
       las=1, ylab="", xlab="",main = t)
  mtext(text = "BPUE", side = 4, line = 3)
  mtext(text = "Salinity (ppt)", side = 1, line = 2.25)
  
  # Use the polygon function to add the CI
  polygon(x = c(new_dat_1$salinity, max(new_dat_1$salinity), 
                rev(new_dat_1$salinity), new_dat_1$salinity[1]),
          y = c(lower, upper[length(upper)],
                rev(upper), lower[1]), border=NA,
          col=rgb(20,200,20,alpha=150,maxColorValue=255))
  
  # Use the lines() function to add the mean model prediction
  lines(x = new_dat_1$salinity, y = fit, lwd=2,
        col=rgb(20,150,20,alpha=255,maxColorValue=255))
  dev.off()
}

############testing
test<-EMP_BPUE%>%filter(Taxlifestage=="Hyperacanthomysis longirostris Adult")