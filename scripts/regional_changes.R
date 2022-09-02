#modeling and figures for regional drought taxa differences
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

library(DroughtData)
WY<-read_excel("data/Water years.xlsx", sheet="yearassignments")

data<-read.csv("data/Drought_region_taxa.csv")%>%
  left_join(WY, by=c("water_year"="Year"))%>%
  mutate(Yr_type=factor(Yr_type, levels=c("Wet", "Above Normal", "Below Normal", "Dry", "Critical")),
         Index_s=(Index-mean(Index))/sd(Index))%>%
  pivot_longer(cols = c(-Region,-water_year, -Index, -Index_s, -Yr_type, -Drought), names_to="Taxa", values_to="BPUE")
nonlog_taxa<-c("Acartiella.sinensis.Adult","All.taxa")

data<-data%>%filter(Region!="North" & Drought!="N")
regions<-unique(data$Region)

aov3_eqn1 = function(region){
  m=aov(log(BPUE)~Drought, data=filter(d1, Region==region))
  pvalue <- list(pvalue = format(summary(m)[[1]][["Pr(>F)"]][1], digits = 3))
}

target_taxa<-c("Daphnia.Adult","Pseudodiaptomus.forbesi.Adult","Hyperacanthomysis.longirostris.Adult",
               "Limnoithona.tetraspina.Adult")
#change factor levels and names for regions
data$Region<-ifelse(data$Region=="SouthCentral","South Central",data$Region)
data$Region<-factor(data$Region,levels=c("Suisun Bay","Suisun Marsh","Confluence","South Central","North"))

taxa<-unique(data$Taxa)
plot_list<-list()
pvalues_list<-list()
for( t in 1:length(target_taxa)){
  tx<-target_taxa[[t]]
  d1<-data%>%filter(Taxa==tx,water_year>1994)
  d1<-d1%>%filter(!is.na(BPUE))
  d1$Region<-as.factor(d1$Region)
  
  bpue_drought_eq1<-map(set_names(unique(d1$Region)), aov3_eqn1)

  bpue_drought_eq1<-data.frame(bpue_drought_eq1)
  colnames(bpue_drought_eq1)<-c("Confluence","South Central","Suisun Bay","Suisun Marsh")
  
  bpue_drought_eq1<- pivot_longer(as.data.frame(bpue_drought_eq1),cols=everything(),names_to = "Region",values_to="p1")
  
  bpue_drought_pvalues<-bpue_drought_eq1
  bpue_drought_pvalues<-cbind(Taxa=tx,bpue_drought_pvalues)
  
  p<-ggplot(d1,aes(Drought,log(BPUE),fill=Drought))+
    geom_boxplot()+
    ggtitle(tx)+
    ylab("ln(BPUE)")+
    theme_classic()+
    facet_grid(cols=vars(factor(Region,levels=c("Suisun Bay","Suisun Marsh","Confluence","South Central"))))+
    drt_color_pal_drought(aes_type = "fill")+
    theme(plot.title = element_text(size=26),axis.title = element_text(size=24),axis.text = element_text(size=24),legend.position = "none",strip.text.x = element_text(size = 20))
  p
  save_plot(paste("figures/regional_differences/",tx,".png",sep=""),p,base_height = 8,base_width = 16)
  plot_list[[t]]=p
  pvalues_list[[t]] <- bpue_drought_pvalues # add it to your list
}
all_plots<-plot_grid(plotlist=plot_list,labels = "AUTO", align = "v",ncol=1)
all_plots
save_plot("figures/regional_differences/all_plots.png",all_plots,base_height = 18,base_width = 12)
all_pvalues<-do.call(rbind,pvalues_list)
write.csv(all_pvalues,"outputs/all_pvalues.csv",row.names = F)

###############
#Lets use the above data to make a heatmap showing percent change

region_means<-data%>%group_by(Taxa,Region,Drought)%>%
  dplyr::summarise(mean_BPUE=mean(BPUE,na.rm=T))
drought_change<-region_means%>%pivot_wider(names_from = Drought,values_from=mean_BPUE)
drought_change$BPUE_change<-drought_change$D-drought_change$W
drought_change$pct_change<-(drought_change$BPUE_change/drought_change$W)*100
drought_change<-drought_change%>%
  filter(Taxa%in%target_taxa)

drought_change<-drought_change%>%
  inner_join(all_pvalues)
#drought_change$pct_change<-ifelse(drought_change$p1<0.05,drought_change$pct_change,NA)
drought_change$Region<-factor(drought_change$Region,levels=c("Suisun Bay","Suisun Marsh","Confluence","South Central"))

drought_change$Taxa<-gsub("[.]"," ",drought_change$Taxa)
drought_change$Taxa<-gsub(" Adult","",drought_change$Taxa)

p<-ggplot(drought_change,aes(x=Taxa,y=Region,fill=pct_change))+
  geom_tile()+
  scale_fill_gradientn(name="% change",colours=c("blue","white","red"))+
  theme(axis.text.x = element_text(angle = 45, vjust =.9, hjust=.9),axis.text = element_text(size=14,face="bold"),axis.title = element_blank())+
  geom_label(aes(label=paste("% change: ",round(pct_change,1))),size=4,fill="white")+
  geom_label(aes(label=paste("p-value: ",p1)),size=4,nudge_y = -.2,fill="white")
p
save_plot("figures/regional_differences/heat_map.png",p,base_height = 9,base_width = 9)
