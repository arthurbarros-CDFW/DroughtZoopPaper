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
library(lme4)
library(stats)

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
Taxlifestage<-c("Daphnia Adult",
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

target_DROUGHT_taxa<-c("Daphnia Adult",
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
  #inner_join(temp_monthlymean)
#model_data$scaled_temp<-model_data$beg_surf_temp-model_data$month_temp

# Function to generate posterior predictions from a gam model
# From https://stats.stackexchange.com/questions/190348/can-i-use-bootstrapping-to-estimate-the-uncertainty-in-a-maximum-value-of-a-gam
predict_posterior<-function(model, newdata, n=1e4, seed=999){
  Xp <- predict(model, newdata=newdata, type="lpmatrix", exclude=c("s(month)", "s(Station)"), newdata.guaranteed=TRUE) ## map coefs to fitted curves
  beta <- coef(model)
  Vb   <- vcov(model) ## posterior mean and cov of coefs
  set.seed(seed)
  mrand <- mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
  pred<-matrix(nrow=nrow(newdata), ncol=n)
  ilink <- family(model)$linkinv
  for (i in seq_len(n)) { 
    pred[,i]   <- ilink(Xp %*% mrand[i, ])
  }
  return(pred)
}

models<-list()
plot_list<-list()
for(i in 1:length(taxa)){
  t<-taxa[i]
  d<-model_data%>%filter(Taxlifestage==t & !is.na(salinity) & year>1993 & month %in% target_months)
  d<-dplyr::select(d,BPUE,salinity,month,Station,Drought)
  d$Station<-as.factor(d$Station)
  d$Drought<-as.factor(d$Drought)
#  m1<-gam(BPUE~s(salinity),family="nb",data=d)
#  summary(m1)
#  m2<-gam(BPUE~s(salinity)+s(month,k=5),family='nb',data=d)
#  summary(m2)
#  #Now with a random effect of station
#  m3<-gamm(BPUE~s(salinity)+s(month,k=5),random = list(Station = ~1), niterPQL=40,family='nb',data=d)
  
  #try two seperate models for p/a binomial and just presence nb
  d_pa<-d
  d_pa$Presence<-ifelse(d_pa$BPUE>0,1,0)
  d_p<-d_pa%>%filter(Presence==1)
  m4.1<-gam(Presence~s(salinity),family="binomial",data=d_pa)
  summary(m4.1)
  models[[taxa[i]]][["presence"]]<-m4.1
  m4.2<-if(t=="Daphnia Adult"){
    gam(BPUE~s(salinity)+s(month,k=5)+ s(Station, bs='re'),niterPQL=40,family='nb',data=d_p)
  }else{
          gam(BPUE~s(salinity)+s(month,k=5), niterPQL=40,family='nb',data=d_p)
    }
  
  summary(m4.2)
  models[[taxa[i]]][["abundance"]]<-m4.2
 
  #capture.output(summary(m1),file=paste("outputs/model_outputs/",t,"_m1.txt"))
  #capture.output(summary(m2),file=paste("outputs/model_outputs/",t,"_m2.txt"))
  #capture.output(summary(m3[[2]]),file=paste("outputs/model_outputs/",t,"_m3.txt"))
  capture.output(summary(m4.1),file=paste("outputs/model_outputs/",t,"_m4_1.txt"))
  capture.output(summary(m4.2),file=paste("outputs/model_outputs/",t,"_m4_2.txt"))
  
  #Anova or AIC to test compare model fits?
  #m_anova<-anova(m1,m2,m3, m4.1)
  #m_anova
  #capture.output(m_anova,file=paste("outputs/model_outputs/",t,"_m_anova.txt"))
  #m_aic<-AIC(m1,m2,m4.1,m4.2)
  #m_aic
  #capture.output(m_aic,file=paste("outputs/model_outputs/",t,"_m_aic.txt"))
  
  png(paste("figures/gamcheck/",t,"_m4_1.png"))
  par(mfrow=c(2,2))
  gam.check(m4.1)
  dev.off()
  
  png(paste("figures/gamcheck/",t,"_m4_2.png"))
  par(mfrow=c(2,2))
  gam.check(m4.2)
  dev.off()
  
  png(paste("figures/gamcheck/",t,"_salspline.png"))
  par(mfrow=c(1,1))
  plot.gam(m4.2,select=1)
  abline(h=0)
  dev.off()
  
  # Generate and combine model predictions
  # Set up salinity vector to predict over
  newdata<-data.frame(salinity=seq(quantile(d$salinity, 0.05), quantile(d$salinity, 0.95), length.out=100))
  # predict for presence model
  m4.1.post<-predict_posterior(m4.1, newdata=newdata)
  # predict for abundance model
  m4.2.post<-predict_posterior(m4.2, newdata=newdata)
  # Multiply them together
  post<-m4.1.post*m4.2.post
  #Collapse to 95% confidence intervals
  preds<-apply(post, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))
  preds_tidy<-newdata%>%
    bind_cols(
      tibble(l95=preds["2.5%",], median=preds["50%",], u95=preds["97.5%",])
    )
  t<-gsub(" Adult","",t)
  p<-ggplot(preds_tidy, aes(x=salinity, ymin=l95, y=median, ymax=u95))+
    geom_ribbon(alpha=0.4)+
    geom_line()+
    ggtitle(t)+
    labs(x="salinity(ppt)",y="BPUE")+
    theme_bw()+
    theme(text = element_text(size = 20))
  p
  plot_list[[i]]=p
  save_plot(paste("figures/gamcheck/",t,"_predict.png",sep=""),p,base_height = 5,base_width = 8)
}
all_plots<-plot_grid(plotlist=plot_list,align = "v",ncol=1)
all_plots
save_plot("figures/gamcheck/all_gams.png", all_plots,ncol=1,base_height = 14,base_width = 9)
