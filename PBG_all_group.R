#Author: Joshua Ajowele####
#This script is for plant biomass and species composition response to fire and grazing heterogeneity
#Date: Feb 6, 2026 Last modified: 

#Load library####
library(tidyverse)
library(ggthemes)
library(readr)
library(performance)
library(car)
library(lme4)
library(lmerTest)
library(nlme)
library(see)
library(patchwork)
library(phia)
### Standard Error function
SE_function<-function(x,na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

#Import plant biomass data####
biomass_diskpasture<-read_csv("Data_PBG_species/PBG031.csv")%>%
  #remove years before 2013 because 2013 was the year the first cycle of PBG
  #treatment was completed across all unit
  filter(Recyear%in%2013:2023)
Diskht_data<-read.csv("Data_PBG_species/PBG032.csv")%>%
  filter(Recyear%in%2013:2023)
#BIOMASS####
##modify data####
biomass_reg <- biomass_diskpasture%>%
  #creating row number in order to visualise outliers in the graph
  mutate(row_num= 1:length(Lvgrass))%>%
  #converting to the appropriate scale g/m2
  mutate(Lvgrass_clean = Lvgrass*10,
         Pdead_clean= Pdead*10,
         Forbs_clean = Forbs*10,
         Woody_clean= Woody*10)%>%
  mutate(total_biomass= Lvgrass_clean+Pdead_clean+Forbs_clean+Woody_clean)%>%
  mutate(total_nonwoody=Lvgrass_clean+Pdead_clean+Forbs_clean)%>%
  mutate(forbs_grass=Lvgrass_clean+Forbs_clean)

##model to predict biomass from diskht####
#Run the regression for each year and have a table of regression equation and t, P value. No need for figures.
m2013<-lm(total_biomass~Diskht, data=biomass_reg[biomass_reg$Recyear==2013,])
summ_2013<-summary(m2013)
check_model(m2013)
coeff_2013<-as.data.frame(summ_2013$coefficients)%>%
  mutate(year=2013, adjR2=0.847, transf="no")
m2014<-lm(total_biomass~Diskht, data=biomass_reg[biomass_reg$Recyear==2014,])
summary(m2014)
check_model(m2014)
summ_2014<-summary(m2014)
coeff_2014<-as.data.frame(summ_2014$coefficients)%>%
  mutate(year=2014, adjR2=0.858, transf="no")
m2015<-lm(sqrt(total_biomass)~Diskht, data=biomass_reg[biomass_reg$Recyear==2015,])
summary(m2015)
check_model(m2015)
summ_2015<-summary(m2015)
coeff_2015<-as.data.frame(summ_2015$coefficients)%>%
  mutate(year=2015, adjR2=0.725, transf="sqrt")
m2016<-lm(sqrt(total_biomass)~Diskht, data=biomass_reg[biomass_reg$Recyear==2016,])
summary(m2016)
check_model(m2016)
summ_2016<-summary(m2016)
coeff_2016<-as.data.frame(summ_2016$coefficients)%>%
  mutate(year=2016, adjR2=0.808, transf="sqrt")
m2017<-lm(total_biomass~Diskht, data=biomass_reg[biomass_reg$Recyear==2017,])
summary(m2017)
check_model(m2017)
summ_2017<-summary(m2017)
coeff_2017<-as.data.frame(summ_2017$coefficients)%>%
  mutate(year=2017, adjR2=0.705, transf="no")
m2018<-lm(sqrt(total_biomass)~Diskht, data=biomass_reg[biomass_reg$Recyear==2018,])
summary(m2018)
check_model(m2018)
summ_2018<-summary(m2018)
coeff_2018<-as.data.frame(summ_2018$coefficients)%>%
  mutate(year=2018, adjR2=0.668, transf="sqrt")
m2019<-lm(total_biomass~Diskht, data=biomass_reg[biomass_reg$Recyear==2019,])
summary(m2019)
check_model(m2019)
summ_2019<-summary(m2019)
coeff_2019<-as.data.frame(summ_2019$coefficients)%>%
  mutate(year=2019, adjR2=0.806, transf="no")
m2020<-lm(log(total_biomass)~Diskht, data=biomass_reg[biomass_reg$Recyear==2020,])
summary(m2020)
check_model(m2020)
summ_2020<-summary(m2020)
coeff_2020<-as.data.frame(summ_2020$coefficients)%>%
  mutate(year=2020, adjR2=0.655, transf="log")
m2021<-lm(sqrt(total_biomass)~Diskht, data=biomass_reg[biomass_reg$Recyear==2021,])
summary(m2021)
check_model(m2021)
summ_2021<-summary(m2021)
coeff_2021<-as.data.frame(summ_2021$coefficients)%>%
  mutate(year=2021, adjR2=0.595, transf="sqrt")
m2022<-lm(total_biomass~Diskht, data=biomass_reg[biomass_reg$Recyear==2022,])
summary(m2022)
check_model(m2022)
summ_2022<-summary(m2022)
coeff_2022<-as.data.frame(summ_2022$coefficients)%>%
  mutate(year=2022, adjR2=0.721, transf="no")
m2023<-lm(sqrt(total_biomass)~Diskht, data=biomass_reg[biomass_reg$Recyear==2023,])
summary(m2023)
check_model(m2023)
summ_2023<-summary(m2023)
coeff_2023<-as.data.frame(summ_2023$coefficients)%>%
  mutate(year=2023, adjR2=0.587, transf="sqrt")

###combine predicted coeff####
coeff_combo<-coeff_2013%>%
  bind_rows(coeff_2014,coeff_2015, coeff_2016, coeff_2017, coeff_2018, 
            coeff_2019, coeff_2020, coeff_2021, coeff_2022, coeff_2023)
#export as csv file
write.csv(coeff_combo, "biomass_diskht_reg.csv")
##estimate biomass from regression equations####
biomass_data<-Diskht_data%>%
  mutate(biomass=case_when(Recyear==2013~23.98823740*Diskht+36.97926816,
                           Recyear==2014~43.04014231*Diskht-141.48476093,
                           Recyear==2015~(0.74268901*Diskht+11.99752162)^2,
                           Recyear==2016~(0.82315555*Diskht+9.44507545)^2,
                           Recyear==2017~34.58429807*Diskht+68.30190206,
                           Recyear==2018~(0.71731156*Diskht+9.73263247)^2,
                           Recyear==2019~23.65683371*Diskht+101.37077552,
                           Recyear==2020~exp(0.07039011*Diskht+4.75326444),
                           Recyear==2021~(0.80860442*Diskht+10.42779960)^2,
                           Recyear==2022~40.00037444*Diskht+79.29615579,
                           Recyear==2023~(0.92699872*Diskht+12.61167802)^2
                           ),
         #create column to match with time since fire key
         year_watershed=paste(Recyear,Watershed,sep="_"))


#Create a watershed key for treatment and Unit column to merge with the raw data
watershed_key <- tibble(Watershed=levels(factor(biomass_data$Watershed)),
                        FireGrzTrt=c("ABG", "PBG", "PBG", "PBG", "ABG", "PBG", "PBG", "PBG"),
                        Unit=c("south", "south", "south", "south", "north", "north",
                               "north", "north"))
#creating a key for year since fire
YrSinceFire_key <- tibble(year_watershed= c("2013_C03A",
                                            "2013_C03B", "2013_C03C", "2013_C1SB",
                                            "2013_C3SA", "2013_C3SB", "2013_C01A",
                                            "2013_C3SC", "2014_C01A", "2014_C03A",
                                            "2014_C03B", "2014_C03C", "2014_C1SB",
                                            "2014_C3SA", "2014_C3SB", "2014_C3SC",
                                            "2015_C1SB", "2015_C03A", "2015_C03B",
                                            "2015_C3SC", "2015_C3SA", "2015_C3SB",
                                            "2015_C03C", "2015_C01A", "2016_C03A",
                                            "2016_C03B", "2016_C3SC", "2016_C01A",
                                            "2016_C03C", "2016_C1SB", "2016_C3SA",
                                            "2016_C3SB", "2017_C3SA", "2017_C3SB",
                                            "2017_C03C", "2017_C03A", "2017_C03B",
                                            "2017_C1SB", "2017_C3SC", "2017_C01A",
                                            "2018_C03A", "2018_C03B", "2018_C03C",
                                            "2018_C01A", "2018_C3SA", "2018_C3SB",
                                            "2018_C1SB", "2018_C3SC", "2019_C3SA",
                                            "2019_C3SB", "2019_C1SB", "2019_C03A",
                                            "2019_C03B", "2019_C03C", "2019_C3SC",
                                            "2019_C01A", "2020_C1SB", "2020_C3SA",
                                            "2020_C3SB", "2020_C3SC", "2020_C03A",
                                            "2020_C03B", "2020_C03C", "2020_C01A",
                                            "2021_C01A", "2021_C03C", "2021_C03A",
                                            "2021_C03B", "2021_C3SA", "2021_C3SB",
                                            "2021_C3SC", "2021_C1SB", 
                                            "2022_C01A", "2022_C1SB", "2022_C03A", "2022_C03B",
                                            "2022_C03C","2022_C3SA","2022_C3SB","2022_C3SC",
                                            "2023_C01A","2023_C1SB","2023_C03A","2023_C03B","2023_C03C",
                                            "2023_C3SA","2023_C3SB","2023_C3SC"),
                          time_fire= c("PBG0", "PBG2", "PBG1", "ABG0","PBG2", "PBG0",
                                         "ABG0", "PBG1","ABG0","PBG1", "PBG0","PBG2","ABG0",
                                         "PBG0","PBG1","PBG2","ABG0","PBG2","PBG1","PBG0","PBG1",
                                         "PBG2","PBG0","ABG0","PBG0","PBG2","PBG1","ABG0","PBG1",
                                         "ABG0","PBG2","PBG0","PBG0","PBG1","PBG2","PBG1","PBG0",
                                         "ABG0","PBG2","ABG0","PBG2","PBG1","PBG0","ABG0","PBG1",
                                         "PBG2","ABG0","PBG0","PBG2","PBG0","ABG0","PBG0","PBG2",
                                         "PBG1","PBG1","ABG0","ABG0","PBG0","PBG1","PBG2","PBG1",
                                         "PBG0","PBG2","ABG0","ABG0","PBG0","PBG2","PBG1","PBG1",
                                         "PBG2","PBG0","ABG0", 
                                       "ABG0","ABG0","PBG0","PBG2","PBG1","PBG2","PBG0","PBG1",
                                       "ABG0","ABG0","PBG1","PBG0","PBG2","PBG0","PBG1","PBG2"))

##combine keys with data and average to the plot level####
biomass_ready<-biomass_data%>%
  left_join(watershed_key, by="Watershed")%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  group_by(Recyear, Unit, Watershed, FireGrzTrt, Transect, Plotnum)%>%
  summarise(biomass=mean(biomass, na.rm=T))%>%
  group_by(Recyear, Unit, Watershed, FireGrzTrt,Transect)%>%
  summarise(biomass=mean(biomass, na.rm=T))
  
##Temporal analysis at the local/transect scale####
#prepare data
temp_biomdata<-biomass_ready%>%
  group_by(Unit, Watershed, FireGrzTrt, Transect)%>%
  summarise(temp_biom=mean(biomass, na.rm=T),
            temp_biom_sd=sd(biomass),
            temp_biom_cv=temp_biom_sd/temp_biom)
#convert charcter to factors
temp_biomdata$Unit=as.factor(temp_biomdata$Unit)
temp_biomdata$Watershed=as.factor(temp_biomdata$Watershed)  
temp_biomdata$Transect=as.factor(temp_biomdata$Transect)
temp_biomdata$FireGrzTrt=as.factor(temp_biomdata$FireGrzTrt)
#start analysis
temp_mean_m<-lmer(temp_biom~FireGrzTrt+(1|Unit/Watershed), data=temp_biomdata)
check_model(temp_mean_m)#outlier
#a good way to visually check outliers from raw values
ggplot(temp_biomdata[temp_biomdata$temp_biom<1000,],aes(FireGrzTrt, temp_biom))+
  geom_boxplot()
#remove outlier and refit model with log transformation to improve model 
temp_mean_m<-lmer(log(temp_biom)~FireGrzTrt+(1|Unit/Watershed), data=temp_biomdata[temp_biomdata$temp_biom<1000,])
anova(temp_mean_m)
check_model(temp_mean_m)
temp_sd_m<-lmer(log(temp_biom_sd)~FireGrzTrt+(1|Unit), data=temp_biomdata[temp_biomdata$temp_biom_sd<1000,])
anova(temp_sd_m)
check_model(temp_sd_m)
temp_cv_m<-lmer(log(temp_biom_cv)~FireGrzTrt+(1|Unit), data=temp_biomdata[temp_biomdata$temp_biom_cv<1,])
check_model(temp_cv_m)
anova(temp_cv_m)

###visuals####
TBM<-interactionMeans(temp_mean_m)%>%
  mutate(resp="biomass_mean")
TBSD<-interactionMeans(temp_sd_m)%>%
  mutate(resp="biomass_sd")
TBCV<-interactionMeans(temp_cv_m)%>%
  mutate(resp="biomass_cv")
temp_biom_vizdata<-TBM%>%
  bind_rows(TBSD,TBCV)%>%
  mutate(avg=exp(`adjusted mean`),
         up=exp(`adjusted mean`+`SE of link`),
         low=exp(`adjusted mean`-`SE of link`))

TBM_viz<-ggplot(temp_biom_vizdata[temp_biom_vizdata$resp=="biomass_mean",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Plant Biomass (gm"^-2*")"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
TBSD_viz<-ggplot(temp_biom_vizdata[temp_biom_vizdata$resp=="biomass_sd",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Plant Biomass SD (gm"^-2*")"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
TBCV_viz<-ggplot(temp_biom_vizdata[temp_biom_vizdata$resp=="biomass_cv",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Plant Biomass CV")+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#Figure to showw each transect across years
ggplot(biomass_ready[biomass_ready$biomass<1200,],aes(Recyear, biomass,col=FireGrzTrt, fill=Watershed, shape=Transect))+
  geom_point(size=5)+
  geom_path()+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Plant Biomass (gm"^-2*")"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )+
  facet_wrap(~Unit)
###combine figures####
TBM_viz+TBSD_viz+TBCV_viz+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")




##Spatial variability at pasture scale####
#wrangle data
spat_biomdata<-biomass_ready%>%
  group_by(Recyear,Unit,FireGrzTrt)%>%
  summarise(spat_biom=mean(biomass, na.rm=T),
            spat_biom_sd=sd(biomass),
            spat_biom_cv=spat_biom_sd/spat_biom)
#convert charcter to factors
spat_biomdata$Unit=as.factor(spat_biomdata$Unit)
spat_biomdata$FireGrzTrt=as.factor(spat_biomdata$FireGrzTrt)
spat_biomdata$Recyear=as.factor(spat_biomdata$Recyear)
#analysis
spat_biommean_m<-lmer(spat_biom~FireGrzTrt*Recyear+(1|Unit), data=spat_biomdata)#remove outliers and refit
check_outliers(spat_biommean_m)
spat_biommean_m<-lmer(log(spat_biom)~FireGrzTrt*Recyear+(1|Unit), data=spat_biomdata[spat_biomdata$spat_biom<1000,])
check_model(spat_biommean_m)
anova(spat_biommean_m)
#check interaction
testInteractions(spat_biommean_m, pairwise="FireGrzTrt",fixed="Recyear")


spat_biomsd_m<-lmer(log(spat_biom_sd)~FireGrzTrt*Recyear+(1|Unit), data=spat_biomdata[spat_biomdata$spat_biom_sd<1000,])
check_model(spat_biomsd_m)
anova(spat_biomsd_m)

spat_biomcv_m<-lmer(log(spat_biom_cv)~FireGrzTrt*Recyear+(1|Unit), data=spat_biomdata[spat_biomdata$spat_biom_cv<1,])
check_model(spat_biomcv_m)
anova(spat_biomcv_m)


##visuals####
biom_mean<-interactionMeans(spat_biommean_m)%>%
  mutate(biomass=exp(`adjusted mean`),
         biomass_up=exp(`adjusted mean`+`SE of link`),
         biomass_low=exp(`adjusted mean`-`SE of link`),
         variab="AVG")

biom_mean_fig<-ggplot(biom_mean,aes(Recyear, biomass,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(Recyear)))+
  geom_errorbar(aes(ymin=biomass_low,
                    ymax=biomass_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Plant Biomass (gm"^-2*")"))+
  xlab(label=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#sd
biom_sd<-interactionMeans(spat_biomsd_m)%>%
  mutate(biomass=exp(`adjusted mean`),
         biomass_up=exp(`adjusted mean`+`SE of link`),
         biomass_low=exp(`adjusted mean`-`SE of link`),
         variab="SD")

biom_sd_fig<-ggplot(biom_sd,aes(Recyear, biomass,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(Recyear)))+
  geom_errorbar(aes(ymin=biomass_low,
                    ymax=biomass_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Plant Biomass SD (gm"^-2*")"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#cv
biom_cv<-interactionMeans(spat_biomcv_m)%>%
  mutate(biomass=exp(`adjusted mean`),
         biomass_up=exp(`adjusted mean`+`SE of link`),
         biomass_low=exp(`adjusted mean`-`SE of link`),
         variab="CV")

biom_cv_fig<-ggplot(biom_cv,aes(Recyear, biomass,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(Recyear)))+
  geom_errorbar(aes(ymin=biomass_low,
                    ymax=biomass_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Plant Biomass CV")+
  xlab(label=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
###combine figures####
biom_mean_fig+biom_sd_fig+biom_cv_fig+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")
###average figures####
#wrangle data
combo_spat_biom<-biom_mean%>%
  bind_rows(biom_sd,biom_cv)%>%
  group_by(FireGrzTrt, variab)%>%
  summarise(biom_avg=mean(biomass, na.rm=T),
            biom_se=SE_function(biomass))
#visual
avg_spat_viz<-ggplot(combo_spat_biom[combo_spat_biom$variab=="AVG",],aes(FireGrzTrt, biom_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=biom_avg-biom_se,
                    ymax=biom_avg+biom_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Plant Biomass (gm"^-2*")"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
sd_spat_viz<-ggplot(combo_spat_biom[combo_spat_biom$variab=="SD",],aes(FireGrzTrt, biom_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=biom_avg-biom_se,
                    ymax=biom_avg+biom_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Plant Biomass (gm"^-2*")"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
cv_spat_viz<-ggplot(combo_spat_biom[combo_spat_biom$variab=="CV",],aes(FireGrzTrt, biom_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=biom_avg-biom_se,
                    ymax=biom_avg+biom_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Plant Biomass CV")+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
avg_spat_viz+sd_spat_viz+cv_spat_viz+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")
