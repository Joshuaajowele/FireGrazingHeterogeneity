#Author: Joshua Ajowele####
#This script is for plant biomass and species composition response to fire and grazing heterogeneity
#Date: Feb 6, 2026 Last modified: March 8, 2026

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
library(vegan)
library(codyn)
library(stringr)
library(readxl)
#install devtools 
#install.packages("devtools")
#library(devtools)
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
### Standard Error function
SE_function<-function(x,na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

#Import data####
#Import plant biomass data
biomass_diskpasture<-read_csv("Data_PBG_species/PBG031.csv")%>%
  #remove years before 2013 because 2013 was the year the first cycle of PBG
  #treatment was completed across all unit
  filter(Recyear%in%2013:2023)
Diskht_data<-read.csv("Data_PBG_species/PBG032.csv")%>%
  filter(Recyear%in%2013:2023)
#import plant composition data
pl_sp_comp<-read_csv("Data_PBG_species/PBG011.csv")%>%
  filter(RecYear%in%2013:2023)
#import plant lifeform class
lifeforms_data<- read_excel("Data_PBG_species/sp_list_update.xlsx")%>%
  mutate(life_form=paste(growthform,lifeform, sep="_"))%>%
  rename(SpeCode=code)%>%
  mutate(life_form=case_when(life_form=="a_s"~"a_g",
                             life_form=="b_f"~"p_f",
                             life_form=="p_s"~"p_g",
                             life_form=="p_o"~"p_f",
                             life_form=="P_w"~"p_w",
                             life_form=="a_m"~"a_f",
                             life_form=="p_m"~"p_f",
                             TRUE~life_form))
#import grasshopper data
grasshopperspcomp_df_raw <- read.csv("Data_PBG_species/PBG073.csv")%>%
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


#keys for treatment and Unit column to merge with the raw data####
watershed_key <- tibble(Watershed=c("C01A", "C03C","C03A","C03B","C1SB","C3SA","C3SB","C3SC"),
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
  group_by(Recyear, Unit, time_fire,Watershed, FireGrzTrt, Transect, Plotnum)%>%
  summarise(biomass=mean(biomass, na.rm=T))%>%
  group_by(Recyear, Unit, time_fire,Watershed, FireGrzTrt,Transect)%>%
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

###combine figures####
TBM_viz+TBSD_viz+TBCV_viz+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")




##Spatial variability at pasture scale####

###time_fire analysis to show differences in PBG patches####
biom_t_fire<-biomass_ready
biom_t_fire$Recyear=as.factor(biom_t_fire$Recyear)
biom_t_fire$Unit=as.factor(biom_t_fire$Unit)
biom_t_fire$time_fire=as.factor(biom_t_fire$time_fire)
biom_t_fire$Watershed=as.factor(biom_t_fire$Watershed)
#model
t_biom_m<-lmer(log(biomass)~time_fire*Recyear+(1|Unit/Watershed), data=biom_t_fire[biom_t_fire$biomass<3000,])
check_model(t_biom_m)
anova(t_biom_m)
testInteractions(t_biom_m,pairwise = "time_fire")
####visual####
t_biom_mean<-interactionMeans(t_biom_m)%>%
  mutate(biomass=exp(`adjusted mean`),
         biomass_up=exp(`adjusted mean`+`SE of link`),
         biomass_low=exp(`adjusted mean`-`SE of link`))

t_fire_biom_fig<-ggplot(t_biom_mean,aes(Recyear, biomass,col=time_fire))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(Recyear)))+
  geom_errorbar(aes(ymin=biomass_low,
                    ymax=biomass_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label=expression("Plant Biomass (gm"^-2*")"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
t_biom_mean_avg<-t_biom_mean%>%
  group_by(time_fire)%>%
  summarise(biomass_avg=mean(biomass, na.rm=T),
            biomass_se=SE_function(biomass))

t_fire_biom_avg_fig<-ggplot(t_biom_mean_avg,aes(time_fire, biomass_avg,col=time_fire))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=biomass_avg-biomass_se,
                    ymax=biomass_avg+biomass_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label=expression("Plant Biomass (gm"^-2*")"))+
  xlab(label="Year since last fire")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

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
  ylab(label=expression("Plant Biomass SD (gm"^-2*")"))+
  xlab(labe="Fire & grazing treatment")+
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
  xlab(label="Fire & grazing treatment")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
t_fire_biom_avg_fig+sd_spat_viz+cv_spat_viz+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")

#PLANT COMPOSITION####
##create required column####
#Creating a key for converting cover class to abundance 
cover_key <-data.frame(CoverClass=0:7, abundance=c(0,0.5,3,15,37.5,62.5,85,97.5))
#need to make watershed & transect all uppercase
pl_sp_comp$Watershed<-toupper(pl_sp_comp$Watershed)
pl_sp_comp$Transect<-toupper(pl_sp_comp$Transect)
#make all genera lowercase for consistency
pl_sp_comp$Ab_genus<-tolower(pl_sp_comp$Ab_genus)
pl_sp_comp$Ab_species<-tolower(pl_sp_comp$Ab_species)

##Temporal plant comp####
##plant species data####
pl_species_comp<-pl_sp_comp%>%
  left_join(cover_key, by="CoverClass")%>%
  #create unique name for species by combining binomial nomenclature
  mutate(sp=paste(Ab_genus,Ab_species, sep="_"))%>%
  group_by(Watershed, RecYear,Transect,Plot,SpeCode,sp)%>%
  #selecting the maximum cover for each species
  summarise(abundance=max(abundance, na.rm=T))%>%
  mutate(year_watershed=paste(RecYear, Watershed, sep="_"))%>%
  #merging key with data
  left_join(YrSinceFire_key, by="year_watershed")%>%
  left_join(watershed_key, by="Watershed")%>%
  group_by(Unit, FireGrzTrt, time_fire, Watershed, RecYear,Transect,SpeCode,sp)%>%
  #average
  summarise(abundance=mean(abundance, na.rm=T))%>%
  group_by(Unit, FireGrzTrt, time_fire, Watershed, RecYear,Transect)%>%
  mutate(total_abund=sum(abundance))%>%
  group_by(Unit, FireGrzTrt, time_fire, Watershed, RecYear,Transect,SpeCode,sp)%>%
  mutate(rel_abund=abundance/total_abund)%>%
  mutate(rep_id=paste(Unit, Watershed, FireGrzTrt, Transect, sep="_"),
         wsd_rep=paste(Watershed,Transect, sep="_"))


##community change calculation yearly in each transect####
time_change<-multivariate_change(pl_species_comp, species.var="sp",
                    abundance.var = "rel_abund",
                    replicate.var="rep_id",
                    time.var = "RecYear",
                    treatment.var="wsd_rep")

#combine to have datatset with tretament variablesfor analysis
t_change<-pl_species_comp%>%
  ungroup()%>%
  select(Unit, Watershed, Transect, FireGrzTrt, time_fire, wsd_rep)%>%
  distinct()%>%
  left_join(time_change, by="wsd_rep")%>%
  mutate(year_diff=paste(RecYear, RecYear2, sep="_"))
t_change$Unit<-as.factor(t_change$Unit)
t_change$FireGrzTrt<-as.factor(t_change$FireGrzTrt)
t_change$Watershed<-as.factor(t_change$Watershed)
t_change$year_diff<-as.factor(t_change$year_diff)

##model for analysis####
t_ch_m<-lmer(log(composition_change)~FireGrzTrt*year_diff+(1|Unit/Watershed), data=t_change)
anova(t_ch_m)
check_model(t_ch_m)
testInteractions(t_ch_m, pairwise = "FireGrzTrt", fixed="year_diff")
##visuals####
#wrangle data from model
comp_change_mean<-interactionMeans(t_ch_m)%>%
  mutate(comp_chang=exp(`adjusted mean`),
         chang_up=exp(`adjusted mean`+`SE of link`),
         chang_low=exp(`adjusted mean`-`SE of link`))

chang_mean_fig<-ggplot(comp_change_mean,aes(year_diff, comp_chang,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(year_diff)))+
  geom_errorbar(aes(ymin=chang_low,
                    ymax=chang_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Composition change"))+
  xlab(label="Year comparison")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#calculate average 
chang_avg<-comp_change_mean%>%
  group_by(FireGrzTrt)%>%
  summarise(comp_chang_avg=mean(comp_chang, na.rm=T),
            c_chang_se=SE_function(comp_chang))
#visual
avg_chang_viz<-ggplot(chang_avg,aes(FireGrzTrt, comp_chang_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=comp_chang_avg-c_chang_se,
                    ymax=comp_chang_avg+c_chang_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Composition change")+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )


##Spatial plant comp####
###transect scale diversity####
pl_rich_trsct <- community_structure(pl_species_comp, time.var = "RecYear", 
                               abundance.var = "rel_abund",
                               replicate.var = "rep_id", metric = "Evar")
#combine with treatmnet info
rich_trsct<-pl_species_comp%>%
  ungroup()%>%
  select(RecYear, Unit, Watershed, Transect, FireGrzTrt, rep_id)%>%
  distinct()%>%
  left_join(pl_rich_trsct, by=c("RecYear","rep_id"))

####analysis####
#convert to factors for analysis
rich_trsct$RecYear<-as.factor(rich_trsct$RecYear)
rich_trsct$Unit<-as.factor(rich_trsct$Unit)
rich_trsct$Watershed<-as.factor(rich_trsct$Watershed)
rich_trsct$FireGrzTrt<-as.factor(rich_trsct$FireGrzTrt)

rich_trsct_m<-lmer(richness~FireGrzTrt*RecYear+(1|Unit/Watershed), data=rich_trsct)
check_model(rich_trsct_m)
anova(rich_trsct_m)
testInteractions(rich_trsct_m, pairwise = "FireGrzTrt", fixed="RecYear")

even_trsct_m<-lmer(Evar~FireGrzTrt*RecYear+(1|Unit/Watershed), data=rich_trsct)
check_model(even_trsct_m)
anova(even_trsct_m)


#####visuals#####
#wrangle data from model
rich_t<-interactionMeans(rich_trsct_m)%>%
  mutate(rich =`adjusted mean`,
         r_up=`adjusted mean`+`SE of link`,
         r_low=`adjusted mean`-`SE of link`)

rich_trsct_fig<-ggplot(rich_t,aes(RecYear, rich,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(RecYear)))+
  geom_errorbar(aes(ymin=r_low,
                    ymax=r_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Richness"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
####time_fire transect####
rich_t_fire_tst<-pl_species_comp%>%
  ungroup()%>%
  select(RecYear, Unit, Watershed, Transect, time_fire, rep_id)%>%
  distinct()%>%
  left_join(pl_rich_trsct, by=c("RecYear","rep_id"))

####analysis####
#convert to factors for analysis
rich_t_fire_tst$RecYear<-as.factor(rich_t_fire_tst$RecYear)
rich_t_fire_tst$Unit<-as.factor(rich_t_fire_tst$Unit)
rich_t_fire_tst$Watershed<-as.factor(rich_t_fire_tst$Watershed)
rich_t_fire_tst$time_fire<-as.factor(rich_t_fire_tst$time_fire)

rich_t_fire_tst_m<-lmer(richness~time_fire*RecYear+(1|Unit/Watershed), data=rich_t_fire_tst)
check_model(rich_t_fire_tst_m)
anova(rich_t_fire_tst_m)

even_t_fire_tst_m<-lmer(Evar~time_fire*RecYear+(1|Unit/Watershed), data=rich_t_fire_tst)
check_model(even_t_fire_tst_m)
anova(even_t_fire_tst_m)
testInteractions(even_t_fire_tst_m, pairwise="time_fire")
ggplot(rich_t_fire_tst, aes(RecYear, richness, col=time_fire))+
  geom_boxplot()

###pasture scale####
####selecting PBG samples covering similar area as ABG####
#wrangle data
pl_comp_pasture_scale<-pl_sp_comp%>%
  left_join(cover_key, by="CoverClass")%>%
  #create unique name for species by combining binomial nomenclature
  mutate(sp=paste(Ab_genus,Ab_species, sep="_"))%>%
  group_by(Watershed, RecYear,Transect,Plot,SpeCode,sp)%>%
  #selecting the maximum cover for each species
  summarise(abundance=max(abundance, na.rm=T))%>%
  mutate(year_watershed=paste(RecYear, Watershed, sep="_"))%>%
  #merging key with data
  left_join(YrSinceFire_key, by="year_watershed")%>%
  left_join(watershed_key, by="Watershed")%>%
  filter(Transect!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Transect!="D")%>%
  filter(Watershed!="C03B" | Transect!="A")%>%
  filter(Watershed!="C03B" | Transect!="B")%>%
  filter(Watershed!="C03B" | Transect!="C")%>%
  filter(Watershed!="C03A" | Transect!="A")%>%
  filter(Watershed!="C03A" | Transect!="B")%>%
  filter(Watershed!="C03A" | Transect!="C")%>%
  filter(Watershed!="C3SC" | Transect!="A")%>%
  filter(Watershed!="C3SC" | Transect!="B")%>%
  filter(Watershed!="C3SC" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="B")%>%
  filter(Watershed!="C3SA" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="D")%>%
  filter(Watershed!="C3SB" | Transect!="A")%>%
  filter(Watershed!="C3SB" | Transect!="B")%>%
  group_by(Unit, FireGrzTrt, RecYear,SpeCode,sp)%>%
  #average
  summarise(abundance=mean(abundance, na.rm=T))%>%
  group_by(Unit, FireGrzTrt, RecYear)%>%
  mutate(total_abund=sum(abundance))%>%
  group_by(Unit, FireGrzTrt,  RecYear,SpeCode,sp)%>%
  mutate(rel_abund=abundance/total_abund)%>%
  mutate(rep_id=paste(Unit, FireGrzTrt, sep="_"))
###pasture/treatment scale diversity
pl_rich_past <- community_structure(pl_comp_pasture_scale, time.var = "RecYear", 
                                    abundance.var = "rel_abund",
                                    replicate.var = "rep_id", metric = "Evar")
#combine with treatmnet info
rich_past<-pl_comp_pasture_scale%>%
  ungroup()%>%
  select(RecYear, Unit, FireGrzTrt, rep_id)%>%
  distinct()%>%
  left_join(pl_rich_past, by=c("RecYear","rep_id"))
rich_past$RecYear=as.factor(rich_past$RecYear)
rich_past$FireGrzTrt=as.factor(rich_past$FireGrzTrt)
rich_past$Unit=as.factor(rich_past$Unit)

####analysis####
rich_past_m<-lmer(richness~FireGrzTrt*RecYear+(1|Unit), data=rich_past)
check_model(rich_past_m)
anova(rich_past_m)
even_past_m<-lmer(Evar~FireGrzTrt*RecYear+(1|Unit), data=rich_past)
check_model(even_past_m)
anova(even_past_m)
####visual####
rich_past_viz<-interactionMeans(rich_past_m)%>%
  mutate(rich =`adjusted mean`,
         r_up=`adjusted mean`+`SE of link`,
         r_low=`adjusted mean`-`SE of link`)%>%
  group_by(FireGrzTrt)%>%
  summarise(richn=mean(rich, na.rm=T),
            r_upp=mean(r_up, na.rm=T),
            r_lower=mean(r_low, na.rm=T))

rich_past_fig<-ggplot(rich_past_viz,aes(FireGrzTrt, richn,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=r_lower,
                    ymax=r_upp),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Richness"))+
  xlab(label="Fire & grazing treatment")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#evenness
even_past_viz<-interactionMeans(even_past_m)%>%
  mutate(Evar =`adjusted mean`,
         e_up=`adjusted mean`+`SE of link`,
         e_low=`adjusted mean`-`SE of link`)%>%
  group_by(FireGrzTrt)%>%
  summarise(Even=mean(Evar, na.rm=T),
            e_upp=mean(e_up, na.rm=T),
            e_lower=mean(e_low, na.rm=T))

evar_past_fig<-ggplot(even_past_viz,aes(FireGrzTrt, Even,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=e_lower,
                    ymax=e_upp),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Evenness"))+
  xlab(label="Fire & grazing treatment")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
####betadiversity####
pl_comp_wide<-pl_species_comp%>%
  filter(Transect!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Transect!="D")%>%
  filter(Watershed!="C03B" | Transect!="A")%>%
  filter(Watershed!="C03B" | Transect!="B")%>%
  filter(Watershed!="C03B" | Transect!="C")%>%
  filter(Watershed!="C03A" | Transect!="A")%>%
  filter(Watershed!="C03A" | Transect!="B")%>%
  filter(Watershed!="C03A" | Transect!="C")%>%
  filter(Watershed!="C3SC" | Transect!="A")%>%
  filter(Watershed!="C3SC" | Transect!="B")%>%
  filter(Watershed!="C3SC" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="B")%>%
  filter(Watershed!="C3SA" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="D")%>%
  filter(Watershed!="C3SB" | Transect!="A")%>%
  filter(Watershed!="C3SB" | Transect!="B")%>%
  ungroup()%>%
  select(-SpeCode,-rel_abund,-total_abund)%>%
  distinct()%>%
  mutate(unit_trt=paste(Unit, FireGrzTrt, sep="_"))%>%
  pivot_wider(names_from = sp, values_from = abundance, values_fill = 0)
  
# Separate pcomp and environmental columns
pl_sp_data <- pl_comp_wide %>%
  ungroup()%>%
  dplyr::select(-1:-9)
pl_env_data <- pl_comp_wide%>%dplyr::select(1:9)
#calculating betadiversity by unit
#####loop;creating a loop to do this####
year_vec_pl <- unique(pl_env_data$RecYear)
pl_perm_unit <- {}
pl_beta_unit <- {}


for(YEAR in 1:length(year_vec_pl)){
  vdist_temp_pl_unit <- vegdist(filter(pl_sp_data, pl_env_data$RecYear ==  year_vec_pl[YEAR]))
  permanova_temp_pl_unit <- adonis2(vdist_temp_pl_unit ~ subset(pl_env_data, RecYear == year_vec_pl[YEAR])$unit_trt)
  permanova_out_temp_pl_unit <- data.frame(RecYear = year_vec_pl[YEAR], 
                                           DF = permanova_temp_pl_unit$Df,
                                           F_value = permanova_temp_pl_unit$`F`,
                                           P_value = permanova_temp_pl_unit$`Pr(>F)`)
  pl_perm_unit <- rbind(pl_perm_unit,permanova_out_temp_pl_unit)
  
  bdisp_temp_pl_unit <- betadisper(vdist_temp_pl_unit, filter(pl_env_data, RecYear==year_vec_pl[YEAR])$unit_trt, type = "centroid")
  bdisp_out_temp_pl_unit <- data.frame(filter(pl_env_data, RecYear==year_vec_pl[YEAR]), distance = bdisp_temp_pl_unit$distances)
  pl_beta_unit <- rbind(pl_beta_unit, bdisp_out_temp_pl_unit)
  
  rm(vdist_temp_pl_unit, permanova_temp_pl_unit, permanova_out_temp_pl_unit, bdisp_temp_pl_unit, bdisp_out_temp_pl_unit)
}
#write.csv(pl_beta_unit, "Data_PBG_species/plant_beta_sep_unit.csv")
#write.csv(pl_perm_unit, "Data_PBG_species/plant_perm_sep_unit.csv")
#pl_beta_unit<-read_csv("Data_PBG_species/plant_beta_sep_unit.csv")
#####analysis####
pl_beta_unit$FireGrzTrt=as.factor(pl_beta_unit$FireGrzTrt)
pl_beta_unit$Unit=as.factor(pl_beta_unit$Unit)
pl_beta_unit$RecYear=as.factor(pl_beta_unit$RecYear)
pl_beta_m<-lmer(distance~FireGrzTrt*RecYear+(1|Unit), data=pl_beta_unit)
check_model(pl_beta_m)
anova(pl_beta_m)
testInteractions(pl_beta_m, pairwise = "FireGrzTrt", fixed="RecYear")
#####visual####
pl_beta_viz<-interactionMeans(pl_beta_m)%>%
  mutate(beta =`adjusted mean`,
         b_up=`adjusted mean`+`SE of link`,
         b_low=`adjusted mean`-`SE of link`)
pl_beta_fig<-ggplot(pl_beta_viz,aes(RecYear, beta,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(RecYear)))+
  geom_errorbar(aes(ymin=b_low,
                    ymax=b_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Betadiversity"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#summarize across years
pl_beta_viz_avg<-pl_beta_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(beta_avg=mean(beta, na.rm=T),
            b_upp=mean(b_up, na.rm=T),
            b_lower=mean(b_low, na.rm=T))
pl_beta_avg_fig<-ggplot(pl_beta_viz_avg,aes(FireGrzTrt, beta_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=b_lower,
                    ymax=b_upp),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Betadiversity"))+
  xlab(label=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

####time since fire composition####
#separate south and north unit
pl_t_fire_comp_s<-pl_species_comp%>%
  filter(Unit=="south")%>%
  ungroup()%>%
  select(-abundance, -total_abund, -rep_id,-wsd_rep,-SpeCode)%>%
  distinct()%>%
  pivot_wider(names_from = sp, values_from = rel_abund, values_fill = 0)
pl_t_fire_comp_n<-pl_species_comp%>%
  filter(Unit=="north")%>%
  ungroup()%>%
  select(-abundance, -total_abund, -rep_id,-wsd_rep,-SpeCode)%>%
  distinct()%>%
  pivot_wider(names_from = sp, values_from = rel_abund, values_fill = 0)
#split environmental and species data
#south
burn_time_sp_data_s <- pl_t_fire_comp_s %>%
  ungroup()%>%
  dplyr::select(-1:-6)
burn_time_env_data_s <- pl_t_fire_comp_s%>%dplyr::select(1:6)
#north
burn_time_sp_data_n <- pl_t_fire_comp_n %>%
  ungroup()%>%
  dplyr::select(-1:-6)
burn_time_env_data_n <- pl_t_fire_comp_n%>%dplyr::select(1:6)

#get nmds1 and 2
burn_time_mds_s <- metaMDS(burn_time_sp_data_s, distance = "bray",k=2) 
burn_time_mds_n <- metaMDS(burn_time_sp_data_n, distance = "bray",k=2) 
#combine NMDS1 and 2 with factor columns and create centroids
burn_time_mds_scores_s <- data.frame(burn_time_env_data_s, scores(burn_time_mds_s, display="sites"))%>%
  group_by(RecYear, Unit,time_fire,Watershed)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))
burn_time_mds_scores_n <- data.frame(burn_time_env_data_n, scores(burn_time_mds_n, display="sites"))%>%
  group_by(RecYear, Unit,time_fire,Watershed)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))

#####visual-plotting centroid through time####
pl_t_s<-ggplot(burn_time_mds_scores_s, aes(x=NMDS1_mean, y=NMDS2_mean, fill=time_fire, shape=Watershed))+
  geom_point(size=8, stroke=2)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  scale_shape_manual(values=c(21:24))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
pl_t_n<-ggplot(burn_time_mds_scores_n, aes(x=NMDS1_mean, y=NMDS2_mean, fill=time_fire, shape=Watershed))+
  geom_point(size=8, stroke=2)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  scale_shape_manual(values=c(21:24))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
 
#####stat #####
#compare by timesince fire include watershed and year
dist_south<-vegdist(burn_time_sp_data_s)
dist_north<-vegdist(burn_time_sp_data_n)
permanova_south<-adonis2(dist_south~burn_time_env_data_s$time_fire+burn_time_env_data_s$Watershed+as.factor(burn_time_env_data_s$RecYear), by="terms")
permanova_south
#multiple comparisons
pairwise_south<-pairwise.adonis2(dist_south~time_fire+Watershed+as.factor(RecYear),by="terms", data=burn_time_env_data_s)

#comparing time since fire within each watershed to validate
wat_south<-vegdist(burn_time_sp_data_s[burn_time_env_data_s$Watershed=="C03C",])
adonis2(wat_south~time_fire+as.factor(RecYear), by="terms", data=burn_time_env_data_s[burn_time_env_data_s$Watershed=="C03C",])
wat_south<-vegdist(burn_time_sp_data_s[burn_time_env_data_s$Watershed=="C03A",])
adonis2(wat_south~time_fire+as.factor(RecYear), by="terms", data=burn_time_env_data_s[burn_time_env_data_s$Watershed=="C03A",])
wat_south<-vegdist(burn_time_sp_data_s[burn_time_env_data_s$Watershed=="C03B",])
adonis2(wat_south~time_fire+as.factor(RecYear), by="terms", data=burn_time_env_data_s[burn_time_env_data_s$Watershed=="C03B",])

permanova_north<-adonis2(dist_north~burn_time_env_data_n$time_fire+burn_time_env_data_n$Watershed+as.factor(burn_time_env_data_n$RecYear), by="terms")
permanova_north
#multiple comparisons
pairwise.adonis2(dist_north~time_fire+Watershed+as.factor(RecYear),by="terms", data=burn_time_env_data_n)


####lifeform####
#wrangle data
lifeforms_data$gen=tolower(lifeforms_data$gen)
lifeforms_data$spec=tolower(lifeforms_data$spec)
lifeforms_data<-lifeforms_data%>%
  mutate(sp=paste(gen, spec, sep="_"))
pl_comp_past_life<-pl_comp_wide%>%
  pivot_longer(10:179, names_to = "sp", values_to = "abundance")%>%
  group_by(Unit,RecYear,FireGrzTrt,sp)%>%
  summarise(abundance=mean(abundance, na.rm=T))%>%
  left_join(lifeforms_data, by="sp")%>%
  mutate(life_form=case_when(sp=="carex_duriu"~"p_g",
                             TRUE~life_form))
  
pl_comp_life<-pl_comp_past_life%>%
  group_by(Unit,RecYear,FireGrzTrt,life_form)%>%
  summarise(abund=sum(abundance))%>%
  group_by(Unit,RecYear,FireGrzTrt)%>%
  mutate(t_abund=sum(abund),
         rel_abund=abund/t_abund)
pl_comp_life$Unit<-as.factor(pl_comp_life$Unit)
pl_comp_life$RecYear<-as.factor(pl_comp_life$RecYear)
pl_comp_life$FireGrzTrt<-as.factor(pl_comp_life$FireGrzTrt)
#####analysis####
#perennial grass
peren_grass<-lmer(rel_abund~FireGrzTrt*RecYear+(1|Unit), data=pl_comp_life[pl_comp_life$life_form=="p_g",])
check_model(peren_grass)
anova(peren_grass)
testInteractions(peren_grass,pairwise="FireGrzTrt", fixed="RecYear")
pgrass_data<-interactionMeans(peren_grass)%>%
  mutate(group="perennial graminoid")
#annual grass
ann_grass<-lmer(sqrt(rel_abund)~FireGrzTrt*RecYear+(1|Unit), data=pl_comp_life[pl_comp_life$life_form=="a_g",])
check_model(ann_grass)
anova(ann_grass)
agrass_data<-interactionMeans(ann_grass)%>%
  mutate(group="annual graminoid")

#perennial forbish
peren_forb<-lmer(rel_abund~FireGrzTrt*RecYear+(1|Unit), data=pl_comp_life[pl_comp_life$life_form=="p_f",])
check_model(peren_forb)
anova(peren_forb)
pforb_data<-interactionMeans(peren_forb)%>%
  mutate(group="perennial forb")
#annual forb
ann_forb<-lmer(rel_abund~FireGrzTrt*RecYear+(1|Unit), data=pl_comp_life[pl_comp_life$life_form=="a_f",])
check_model(ann_forb)
anova(ann_forb)
aforb_data<-interactionMeans(ann_forb)%>%
  mutate(group="annual forb")
#woody
woody<-lmer(rel_abund~FireGrzTrt*RecYear+(1|Unit), data=pl_comp_life[pl_comp_life$life_form=="p_w",])
check_model(woody)
anova(woody)
ggplot(pl_comp_life[pl_comp_life$life_form=="p_w",], aes(RecYear, rel_abund, col=FireGrzTrt))+
  geom_boxplot()
woody_data<-interactionMeans(woody)%>%
  mutate(group="woody")
#####visual####
group_data<-pgrass_data%>%
  bind_rows(agrass_data,aforb_data,pforb_data, woody_data)%>%
  mutate(upper=(`adjusted mean`+`SE of link`),
         lower=(`adjusted mean`-`SE of link`),
         grp_mean=case_when(group=="annual graminoid"~(`adjusted mean`)^2,
                                                      TRUE~`adjusted mean`),
         grp_upper=case_when(group=="annual graminoid"~(upper)^2,
                                                          TRUE~upper),
         grp_lower=case_when(group=="annual graminoid"~(lower)^2,
                                               TRUE~lower))
grp_data_ready<-group_data%>%
  group_by(group, FireGrzTrt)%>%
  summarise(grp_m=mean(grp_mean, na.rm=T),
            grp_upper=mean(grp_upper, na.rm=T),
            grp_lower=mean(grp_lower, na.rm=T))
Pl_group<-ggplot(grp_data_ready,aes(group, grp_m,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=grp_lower,
                    ymax=grp_upper),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Relative cover")+
  xlab(label="Plant group")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

#are woody relative cover increasing with year
woody_l<-lmer(rel_abund~FireGrzTrt+RecYear+(1|Unit), data=pl_comp_life[pl_comp_life$life_form=="p_w",])
woody_l<-lm(log(rel_abund)~FireGrzTrt*as.numeric(RecYear),data=pl_comp_life[pl_comp_life$life_form=="p_w",])
anova(woody_l)
check_model(woody_l)
summary(woody_l)
ggplot(data=pl_comp_life[pl_comp_life$life_form=="p_w",], aes(as.numeric(RecYear), rel_abund, col=FireGrzTrt))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Unit)
#result much pronounced in south unit which is the unit we have bird data for
data=pl_comp_life[pl_comp_life$life_form=="p_w",]
woody_l<-lm(rel_abund~FireGrzTrt*as.numeric(RecYear),data=data[data$Unit=="south",])
summary(woody_l)#0.889-Adjsuted r square
check_model(woody_l)
ggplot(data[data$Unit=="south",], aes(RecYear, rel_abund, col=FireGrzTrt))+
  geom_point(size=3)+
  geom_smooth(method="lm", aes(as.numeric(RecYear)))+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Relative woody cover")+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

#GRASSHOPPER####
#wrangle data
grasshopperspcomp_df <- grasshopperspcomp_df_raw %>% 
  mutate(Watershed = replace (Watershed, Watershed == "0c1a", "C01A"),
         Watershed = replace (Watershed, Watershed == "0c3a", "SC03A"),
         Watershed = replace (Watershed, Watershed == "0c3b", "C03B"),
         Watershed = replace (Watershed, Watershed == "0c3c", "C03C"),
         Watershed = replace (Watershed, Watershed == "c01a", "C01A"),
         Watershed = replace (Watershed, Watershed == "c03a", "C03A"),
         Watershed = replace (Watershed, Watershed == "c03b", "C03B"),
         Watershed = replace (Watershed, Watershed == "c03c", "C03C"),
         Watershed = replace (Watershed, Watershed == "c1sb", "C1SB"),
         Watershed = replace (Watershed, Watershed == "c3sa", "C3SA"),
         Watershed = replace (Watershed, Watershed == "c3sb", "C3SB"),
         Watershed = replace (Watershed, Watershed == "c3sc", "C3SC"),
         Watershed = replace (Watershed, Watershed == "c1a", "C01A"),
         Watershed = replace (Watershed, Watershed == "c3a", "C03A"),
         Watershed = replace (Watershed, Watershed == "c3b", "C03B"),
         Watershed = replace (Watershed, Watershed == "c3c", "C03C"))

#converting all repsite into uppercase
grasshopperspcomp_df$Repsite=toupper(grasshopperspcomp_df$Repsite)
#checking species 
species_list<-unique(grasshopperspcomp_df$Species)%>%
  tibble()
species_list<-unique(grasshopperspcomp_df$Species[grasshopperspcomp_df$Recyear>=2022])%>%
  tibble()
#merging key with dataset and removing katydid, etc
grassh_count_df <- grasshopperspcomp_df%>%
  left_join(watershed_key, by = "Watershed")%>%
  mutate(RecYear = Recyear)%>%
  filter(!Species%in%c("Oecanthinae spp.","Tettigoniidae spp.","Gryllidae spp.",
                   "Conocephalus spp.","Neoconocephalus robustus","Scudderia texensis",
                   "Arethaea constricta","Orchelimum spp.","Amblycorypha oblongifolia","Pediodectes haldemani",
                   "Amblycorypha rotundifolia","Neoconocephalus spp.","Neoconocephalus ensiger","Pediodectes nigromarginatus",
                   "Scudderia furcata","Scudderia spp."))%>%
  #using the max cover from the two sweeps done on each transect
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite,Species)%>%
  summarise(Total=max(Total))%>%
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite)%>%
  #calculate total count per repsite
  summarise(Tcount=sum(Total, na.rm=T))%>%
  ungroup()
##Temporal analysis at the local/transect scale####
#prepare data
g_temp_countdata<-grassh_count_df%>%
  group_by(Unit, Watershed, FireGrzTrt, Repsite)%>%
  summarise(temp_count=mean(Tcount, na.rm=T),
            temp_count_sd=sd(Tcount),
            temp_count_cv=temp_count_sd/temp_count)
#convert charcter to factors
g_temp_countdata$Unit=as.factor(g_temp_countdata$Unit)
g_temp_countdata$Watershed=as.factor(g_temp_countdata$Watershed)  
g_temp_countdata$Repsite=as.factor(g_temp_countdata$Repsite)
g_temp_countdata$FireGrzTrt=as.factor(g_temp_countdata$FireGrzTrt)
#start analysis
#temporal mean
g_temp_mean_m<-lmer(temp_count~FireGrzTrt+(1|Unit/Watershed), data=g_temp_countdata)
check_model(g_temp_mean_m)
anova(g_temp_mean_m)
#temporal variability
g_temp_sd_m<-lmer(log(temp_count_sd)~FireGrzTrt+(1|Unit/Watershed), data=g_temp_countdata)
check_model(g_temp_sd_m)
check_outliers(g_temp_sd_m)
anova(g_temp_sd_m)
#remove outlier from data
ggg<-g_temp_countdata%>%
   filter(Watershed!="C3SC"|Repsite!="D")
#run without outlier
g_temp_sd_m<-lmer(log(temp_count_sd)~FireGrzTrt+(1|Unit/Watershed), data=ggg)
check_model(g_temp_sd_m)
anova(g_temp_sd_m)#result similar with/without outlier

g_temp_cv_m<-lmer(log(temp_count_cv)~FireGrzTrt+(1|Unit/Watershed), data=g_temp_countdata)
check_model(g_temp_cv_m)
anova(g_temp_cv_m)

##visuals####
g_TBM<-interactionMeans(g_temp_mean_m)%>%
  mutate(resp="count_mean")
g_TBSD<-interactionMeans(g_temp_sd_m)%>%
  mutate(resp="count_sd")
g_TBCV<-interactionMeans(g_temp_cv_m)%>%
  mutate(resp="count_cv")
g_temp_count_vizdata<-g_TBM%>%
  bind_rows(g_TBSD,g_TBCV)%>%
  mutate(avg=case_when(resp=="count_sd"~exp(`adjusted mean`),
                      resp=="count_cv"~exp(`adjusted mean`),
                      TRUE~`adjusted mean`),
         up=case_when(resp=="count_sd"~exp(`adjusted mean`+`SE of link`),
                     resp=="count_cv"~exp(`adjusted mean`+`SE of link`),
                     TRUE~(`adjusted mean`+`SE of link`)),
         low=case_when(resp=="count_sd"~exp(`adjusted mean`-`SE of link`),
                      resp=="count_cv"~exp(`adjusted mean`-`SE of link`),
                      TRUE~(`adjusted mean`-`SE of link`)))

g_TBM_viz<-ggplot(g_temp_count_vizdata[g_temp_count_vizdata$resp=="count_mean",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper abundance"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
g_TBCV_viz<-ggplot(g_temp_count_vizdata[g_temp_count_vizdata$resp=="count_sd",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper SD"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

g_TBCV_viz<-ggplot(g_temp_count_vizdata[g_temp_count_vizdata$resp=="count_cv",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper CV"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
###combine figures####
g_TBM_viz+g_TBSD_viz+g_TBCV_viz+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")

##Spatial variability at pasture scale####
grass_count_tfire<-grassh_count_df%>%
  mutate(year_watershed=paste(RecYear, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")

###time_fire analysis to show differences in PBG patches####
grass_count_tfire$RecYear=as.factor(grass_count_tfire$RecYear)
grass_count_tfire$Unit=as.factor(grass_count_tfire$Unit)
grass_count_tfire$time_fire=as.factor(grass_count_tfire$time_fire)
grass_count_tfire$Watershed=as.factor(grass_count_tfire$Watershed)
#model
g_t_count_m<-lmer(log(Tcount)~time_fire*RecYear+(1|Unit/Watershed), data=grass_count_tfire)
check_model(g_t_count_m)
anova(g_t_count_m)
testInteractions(g_t_count_m,pairwise = "time_fire", fixed="RecYear")
testInteractions(g_t_count_m,pairwise = "time_fire")

g_t_count_mean<-interactionMeans(g_t_count_m)%>%
  mutate(count_m=exp(`adjusted mean`),
         count_up=exp(`adjusted mean`+`SE of link`),
         count_low=exp(`adjusted mean`-`SE of link`))

g_t_fire_count_fig<-ggplot(g_t_count_mean,aes(RecYear, count_m,col=time_fire))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(RecYear)))+
  geom_errorbar(aes(ymin=count_low,
                    ymax=count_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label=expression("Grasshopper abundance"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
g_t_c_mean_avg<-g_t_count_mean%>%
  group_by(time_fire)%>%
  summarise(count_avg=mean(count_m, na.rm=T),
            count_se=SE_function(count_m))

g_t_fire_c_avg_fig<-ggplot(g_t_c_mean_avg,aes(time_fire, count_avg,col=time_fire))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=count_avg-count_se,
                    ymax=count_avg+count_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label=expression("Grasshopper abundance"))+
  xlab(label="Year since last fire")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#wrangle data
g_spat_cdata<-grass_count_tfire%>%
  group_by(RecYear,Unit,FireGrzTrt)%>%
  summarise(spat_c=mean(Tcount, na.rm=T),
            spat_c_sd=sd(Tcount),
            spat_c_cv=spat_c_sd/spat_c)
#convert character to factors
g_spat_cdata$Unit=as.factor(g_spat_cdata$Unit)
g_spat_cdata$FireGrzTrt=as.factor(g_spat_cdata$FireGrzTrt)
g_spat_cdata$RecYear=as.factor(g_spat_cdata$RecYear)
#analysis
g_spat_c_m<-lmer(spat_c~FireGrzTrt*RecYear+(1|Unit), data=g_spat_cdata)
check_model(g_spat_c_m)
anova(g_spat_c_m)
testInteractions(g_spat_c_m, pairwise = "FireGrzTrt", fixed="RecYear")
#SD
g_spat_c_sd<-lmer(log(spat_c_sd)~FireGrzTrt*RecYear+(1|Unit), data=g_spat_cdata)
check_model(g_spat_c_sd)
anova(g_spat_c_sd)
#cv
g_spat_c_cv<-lmer(spat_c_cv~FireGrzTrt*RecYear+(1|Unit), data=g_spat_cdata)
check_model(g_spat_c_cv)
anova(g_spat_c_sd)
