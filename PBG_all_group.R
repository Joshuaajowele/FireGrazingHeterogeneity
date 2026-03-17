#Author: Joshua Ajowele####
#This script is for plant biomass and species composition response to fire and grazing heterogeneity
#Date: Feb 6, 2026 Last modified: March 11, 2026

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
#grasshopper feeding guild
feeding_df<-read_excel("Data_PBG_species/Grasshopper_guild.xlsx")
#Bird data
PBG_bird_data_raw<-read_excel("Data_PBG_species/Qy_ExportPBGSurveyData_Mar2023.xlsx")
PBG_bird_C1A<-read_excel("Data_PBG_species/Qy_ExportPBGSurveyData_C1A_Apr2023.xlsx")%>%
  filter(Year%in%2011:2021)
PBG_bird_updated<-read_excel("Data_PBG_species/Qy_ExportPBGSurveyData_Feb2026.xlsx")
unique(PBG_bird_C1A$Year)
#2019 data is missing for C1A-talk to Alice
PBG_bird_class<-read_excel("Data_PBG_species/KNZSpecies.xlsx")
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
###combine####
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
         Watershed = replace (Watershed, Watershed == "0c3a", "C03A"),
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
unique(grasshopperspcomp_df$Watershed)
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
g_TBSD_viz<-ggplot(g_temp_count_vizdata[g_temp_count_vizdata$resp=="count_sd",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
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
anova(g_spat_c_cv)
testInteractions(g_spat_c_cv,pairwise = "FireGrzTrt")
##visuals####
gcount_mean<-interactionMeans(g_spat_c_m)%>%
  mutate(count_m=(`adjusted mean`),
         c_up=(`adjusted mean`+`SE of link`),
         c_low=(`adjusted mean`-`SE of link`),
         variab="AVG")

count_mean_fig<-ggplot(gcount_mean,aes(RecYear, count_m,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(RecYear)))+
  geom_errorbar(aes(ymin=c_low,
                    ymax=c_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper abundance"))+
  xlab(label=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
gcount_sd<-interactionMeans(g_spat_c_sd)%>%
  mutate(count_m=exp(`adjusted mean`),
         c_up=exp(`adjusted mean`+`SE of link`),
         c_low=exp(`adjusted mean`-`SE of link`),
         variab="SD")

count_sd_fig<-ggplot(gcount_sd,aes(RecYear, count_m,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(RecYear)))+
  geom_errorbar(aes(ymin=c_low,
                    ymax=c_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper SD"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
gcount_cv<-interactionMeans(g_spat_c_cv)%>%
  mutate(count_m=(`adjusted mean`),
         c_up=(`adjusted mean`+`SE of link`),
         c_low=(`adjusted mean`-`SE of link`),
         variab="CV")

count_cv_fig<-ggplot(gcount_cv,aes(RecYear, count_m,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(RecYear)))+
  geom_errorbar(aes(ymin=c_low,
                    ymax=c_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper CV"))+
  xlab(label=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
###combine figures####
count_mean_fig+count_sd_fig+count_cv_fig+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")
###average figures####
#wrangle data
combo_spat_gcount<-gcount_mean%>%
  bind_rows(gcount_sd,gcount_cv)%>%
  group_by(FireGrzTrt, variab)%>%
  summarise(gcount_avg=mean(count_m, na.rm=T),
            c_se=SE_function(count_m))
#visual
g_avg_spat_viz<-ggplot(combo_spat_gcount[combo_spat_gcount$variab=="AVG",],aes(FireGrzTrt, gcount_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=gcount_avg-c_se,
                    ymax=gcount_avg+c_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper abundance"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

g_sd_spat_viz<-ggplot(combo_spat_gcount[combo_spat_gcount$variab=="SD",],aes(FireGrzTrt, gcount_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=gcount_avg-c_se,
                    ymax=gcount_avg+c_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper SD"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

g_cv_spat_viz<-ggplot(combo_spat_gcount[combo_spat_gcount$variab=="CV",],aes(FireGrzTrt, gcount_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=gcount_avg-c_se,
                    ymax=gcount_avg+c_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper CV"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
###combine####
g_t_fire_c_avg_fig+g_sd_spat_viz+g_cv_spat_viz+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")

#G composition#####
#merging key with dataset and removing katydid, etc
grassh_comm_df <- grasshopperspcomp_df%>%
  left_join(watershed_key, by = "Watershed")%>%
  mutate(RecYear = Recyear)%>%
  filter(!Species%in%c("Oecanthinae spp.","Tettigoniidae spp.","Gryllidae spp.",
                       "Conocephalus spp.","Neoconocephalus robustus","Scudderia texensis",
                       "Arethaea constricta","Orchelimum spp.","Amblycorypha oblongifolia","Pediodectes haldemani",
                       "Amblycorypha rotundifolia","Neoconocephalus spp.","Neoconocephalus ensiger","Pediodectes nigromarginatus",
                       "Scudderia furcata","Scudderia spp."))%>%
  #using the max count from the two diff survey done on each transect
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite,Species)%>%
  summarise(gcount=max(Total))
grassh_tcomm_df<-grassh_comm_df%>%
  group_by(Unit, FireGrzTrt,Watershed, RecYear,Repsite)%>%
  mutate(total_count=sum(gcount))%>%
  group_by(Unit, FireGrzTrt, Watershed, RecYear,Repsite,Species)%>%
  mutate(rel_abund=gcount/total_count)%>%
  mutate(rep_id=paste(Unit, Watershed, FireGrzTrt, Repsite, sep="_"),
         wsd_rep=paste(Watershed,Repsite, sep="_"))

##community change calculation yearly in each transect####
g_time_change<-multivariate_change(grassh_tcomm_df, species.var="Species",
                                 abundance.var = "rel_abund",
                                 replicate.var="rep_id",
                                 time.var = "RecYear",
                                 treatment.var="wsd_rep")

#combine to have datatset with tretament variables for analysis
g_t_change<-grassh_tcomm_df%>%
  ungroup()%>%
  select(Unit, Watershed, Repsite, FireGrzTrt, wsd_rep)%>%
  distinct()%>%
  left_join(g_time_change, by="wsd_rep")%>%
  mutate(year_diff=paste(RecYear, RecYear2, sep="_"))
g_t_change$Unit<-as.factor(g_t_change$Unit)
g_t_change$FireGrzTrt<-as.factor(g_t_change$FireGrzTrt)
g_t_change$Watershed<-as.factor(g_t_change$Watershed)
g_t_change$year_diff<-as.factor(g_t_change$year_diff)

##model for analysis####
g_t_ch_m<-lmer(composition_change~FireGrzTrt*year_diff+(1|Unit/Watershed), data=g_t_change)
check_model(g_t_ch_m)
anova(g_t_ch_m)
testInteractions(g_t_ch_m, pairwise = "FireGrzTrt")
##visuals####
#wrangle data from model
g_comp_change_mean<-interactionMeans(g_t_ch_m)%>%
  mutate(comp_chang=(`adjusted mean`),
         chang_up=(`adjusted mean`+`SE of link`),
         chang_low=(`adjusted mean`-`SE of link`))

g_chang_mean_fig<-ggplot(g_comp_change_mean,aes(year_diff, comp_chang,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(year_diff)))+
  geom_errorbar(aes(ymin=chang_low,
                    ymax=chang_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper composition change"))+
  xlab(label="Year comparison")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#calculate average 
g_chang_avg<-g_comp_change_mean%>%
  group_by(FireGrzTrt)%>%
  summarise(comp_chang_avg=mean(comp_chang, na.rm=T),
            c_chang_se=SE_function(comp_chang))
#visual
g_avg_chang_viz<-ggplot(g_chang_avg,aes(FireGrzTrt, comp_chang_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=comp_chang_avg-c_chang_se,
                    ymax=comp_chang_avg+c_chang_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Grasshopper composition change")+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
##Spatial grassh comp####
###transect scale diversity####
g_rich_trsct <- community_structure(grassh_tcomm_df, time.var = "RecYear", 
                                     abundance.var = "rel_abund",
                                     replicate.var = "rep_id", metric = "Evar")
#combine with treatmnet info
g_rich_trsct<-grassh_tcomm_df%>%
  ungroup()%>%
  select(RecYear, Unit, Watershed, Repsite, FireGrzTrt, rep_id)%>%
  distinct()%>%
  left_join(g_rich_trsct, by=c("RecYear","rep_id"))

####analysis####
#convert to factors for analysis
g_rich_trsct$RecYear<-as.factor(g_rich_trsct$RecYear)
g_rich_trsct$Unit<-as.factor(g_rich_trsct$Unit)
g_rich_trsct$Watershed<-as.factor(g_rich_trsct$Watershed)
g_rich_trsct$FireGrzTrt<-as.factor(g_rich_trsct$FireGrzTrt)

g_rich_trsct_m<-lmer(richness~FireGrzTrt*RecYear+(1|Unit/Watershed), data=g_rich_trsct)
check_model(g_rich_trsct_m)
anova(g_rich_trsct_m)
testInteractions(g_rich_trsct_m, pairwise = "FireGrzTrt")

g_even_trsct_m<-lmer(Evar~FireGrzTrt*RecYear+(1|Unit/Watershed), data=g_rich_trsct)
check_model(g_even_trsct_m)
anova(g_even_trsct_m)
testInteractions(g_even_trsct_m, pairwise = "FireGrzTrt")

#####visuals#####
#wrangle data from model
g_rich_t<-interactionMeans(g_rich_trsct_m)%>%
  mutate(rich =`adjusted mean`,
         r_up=`adjusted mean`+`SE of link`,
         r_low=`adjusted mean`-`SE of link`)

g_rich_trsct_fig<-ggplot(g_rich_t,aes(RecYear, rich,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(RecYear)))+
  geom_errorbar(aes(ymin=r_low,
                    ymax=r_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper richness"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
####time_fire transect####
grassh_tcomm_df$RecYear<-as.factor(grassh_tcomm_df$RecYear)
g_rich_t_fire_tst<-grassh_tcomm_df%>%
  ungroup()%>%
  mutate(year_watershed=paste(RecYear, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  select(RecYear, Unit, Watershed, Repsite, time_fire, rep_id)%>%
  distinct()%>%
  left_join(g_rich_trsct, by=c("RecYear","rep_id","Watershed","Unit","Repsite"))

####analysis####
#convert to factors for analysis
g_rich_t_fire_tst$Unit<-as.factor(g_rich_t_fire_tst$Unit)
g_rich_t_fire_tst$Watershed<-as.factor(g_rich_t_fire_tst$Watershed)
g_rich_t_fire_tst$time_fire<-as.factor(g_rich_t_fire_tst$time_fire)

g_rich_t_fire_tst_m<-lmer(richness~time_fire*RecYear+(1|Unit/Watershed), data=g_rich_t_fire_tst)
check_model(g_rich_t_fire_tst_m)
anova(g_rich_t_fire_tst_m)
testInteractions(g_rich_t_fire_tst_m, pairwise="time_fire", fixed="RecYear")
testInteractions(g_rich_t_fire_tst_m, pairwise="time_fire")

g_even_t_fire_tst_m<-lmer(Evar~time_fire*RecYear+(1|Unit/Watershed), data=g_rich_t_fire_tst)
check_model(g_even_t_fire_tst_m)
anova(g_even_t_fire_tst_m)
testInteractions(g_even_t_fire_tst_m, pairwise="time_fire", fixed="RecYear")
testInteractions(g_even_t_fire_tst_m, pairwise="time_fire")
ggplot(g_rich_t_fire_tst, aes(RecYear, richness, col=time_fire))+
  geom_boxplot()
ggplot(g_rich_t_fire_tst, aes(RecYear, Evar, col=time_fire))+
  geom_boxplot()

###pasture scale####
####selecting PBG samples covering similar area as ABG####
#wrangle data
g_comp_pasture<-grassh_comm_df%>%
  filter(Repsite!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Repsite!="D")%>%
  filter(Watershed!="C03B" | Repsite!="A")%>%
  filter(Watershed!="C03B" | Repsite!="B")%>%
  filter(Watershed!="C03B" | Repsite!="C")%>%
  filter(Watershed!="C03A" | Repsite!="A")%>%
  filter(Watershed!="C03A" | Repsite!="B")%>%
  filter(Watershed!="C03A" | Repsite!="C")%>%
  filter(Watershed!="C3SC" | Repsite!="A")%>%
  filter(Watershed!="C3SC" | Repsite!="B")%>%
  filter(Watershed!="C3SC" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="B")%>%
  filter(Watershed!="C3SA" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="D")%>%
  filter(Watershed!="C3SB" | Repsite!="A")%>%
  filter(Watershed!="C3SB" | Repsite!="B")%>%
  group_by(Unit, FireGrzTrt, RecYear,Species)%>%
  #average
  summarise(abundance=mean(gcount, na.rm=T))%>%
  group_by(Unit, FireGrzTrt, RecYear)%>%
  mutate(total_abund=sum(abundance))%>%
  group_by(Unit, FireGrzTrt,  RecYear,Species)%>%
  mutate(rel_abund=abundance/total_abund)%>%
  mutate(rep_id=paste(Unit, FireGrzTrt, sep="_"))
###pasture/treatment scale diversity
g_rich_past <- community_structure(g_comp_pasture, time.var = "RecYear", 
                                    abundance.var = "rel_abund",
                                    replicate.var = "rep_id", metric = "Evar")
#combine with treatmnet info
g_rich_past<-g_comp_pasture%>%
  ungroup()%>%
  select(RecYear, Unit, FireGrzTrt, rep_id)%>%
  distinct()%>%
  left_join(g_rich_past, by=c("RecYear","rep_id"))
g_rich_past$RecYear=as.factor(g_rich_past$RecYear)
g_rich_past$FireGrzTrt=as.factor(g_rich_past$FireGrzTrt)
g_rich_past$Unit=as.factor(g_rich_past$Unit)

####analysis####
g_rich_past_m<-lmer(richness~FireGrzTrt*RecYear+(1|Unit), data=g_rich_past)
check_model(g_rich_past_m)
anova(g_rich_past_m)
testInteractions(g_rich_past_m, pairwise="FireGrzTrt")
g_even_past_m<-lmer(Evar~FireGrzTrt*RecYear+(1|Unit), data=g_rich_past)
check_model(g_even_past_m)
anova(g_even_past_m)
testInteractions(g_even_past_m, pairwise="FireGrzTrt")
####visual####
g_rich_past_viz<-interactionMeans(g_rich_past_m)%>%
  mutate(rich =`adjusted mean`,
         r_up=`adjusted mean`+`SE of link`,
         r_low=`adjusted mean`-`SE of link`)%>%
  group_by(FireGrzTrt)%>%
  summarise(richn=mean(rich, na.rm=T),
            r_upp=mean(r_up, na.rm=T),
            r_lower=mean(r_low, na.rm=T))

g_rich_past_fig<-ggplot(g_rich_past_viz,aes(FireGrzTrt, richn,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=r_lower,
                    ymax=r_upp),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper richness"))+
  xlab(label="Fire & grazing treatment")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#evenness
g_even_past_viz<-interactionMeans(g_even_past_m)%>%
  mutate(Evar =`adjusted mean`,
         e_up=`adjusted mean`+`SE of link`,
         e_low=`adjusted mean`-`SE of link`)%>%
  group_by(FireGrzTrt)%>%
  summarise(Even=mean(Evar, na.rm=T),
            e_upp=mean(e_up, na.rm=T),
            e_lower=mean(e_low, na.rm=T))

g_evar_past_fig<-ggplot(g_even_past_viz,aes(FireGrzTrt, Even,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=e_lower,
                    ymax=e_upp),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper evenness"))+
  xlab(label="Fire & grazing treatment")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
####betadiversity####
g_comp_wide<-grassh_tcomm_df%>%
  mutate(year_watershed=paste(RecYear, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  filter(Repsite!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Repsite!="D")%>%
  filter(Watershed!="C03B" | Repsite!="A")%>%
  filter(Watershed!="C03B" | Repsite!="B")%>%
  filter(Watershed!="C03B" | Repsite!="C")%>%
  filter(Watershed!="C03A" | Repsite!="A")%>%
  filter(Watershed!="C03A" | Repsite!="B")%>%
  filter(Watershed!="C03A" | Repsite!="C")%>%
  filter(Watershed!="C3SC" | Repsite!="A")%>%
  filter(Watershed!="C3SC" | Repsite!="B")%>%
  filter(Watershed!="C3SC" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="B")%>%
  filter(Watershed!="C3SA" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="D")%>%
  filter(Watershed!="C3SB" | Repsite!="A")%>%
  filter(Watershed!="C3SB" | Repsite!="B")%>%
  ungroup()%>%
  select(-total_count:-year_watershed)%>%
  distinct()%>%
  mutate(unit_trt=paste(Unit, FireGrzTrt, sep="_"))%>%
  pivot_wider(names_from = Species, values_from = gcount, values_fill = 0)

# Separate pcomp and environmental columns
g_sp_data <- g_comp_wide %>%
  ungroup()%>%
  dplyr::select(-1:-7)
g_env_data <- g_comp_wide%>%dplyr::select(1:7)
#calculating betadiversity by unit
#####loop;creating a loop to do this####
year_vec_g <- unique(g_env_data$RecYear)
g_perm_unit <- {}
g_beta_unit <- {}


for(YEAR in 1:length(year_vec_g)){
  vdist_temp_g_unit <- vegdist(filter(g_sp_data, g_env_data$RecYear ==  year_vec_g[YEAR]))
  bdisp_temp_g_unit <- betadisper(vdist_temp_g_unit, filter(g_env_data, RecYear==year_vec_g[YEAR])$unit_trt, type = "centroid")
  bdisp_out_temp_g_unit <- data.frame(filter(g_env_data, RecYear==year_vec_g[YEAR]), distance = bdisp_temp_g_unit$distances)
  g_beta_unit <- rbind(g_beta_unit, bdisp_out_temp_g_unit)
  
  rm(vdist_temp_g_unit, bdisp_temp_g_unit, bdisp_out_temp_g_unit)
}
#write.csv(g_beta_unit, "Data_PBG_species/g_beta_sep_unit.csv")
#####analysis####
g_beta_unit$FireGrzTrt=as.factor(g_beta_unit$FireGrzTrt)
g_beta_unit$Unit=as.factor(g_beta_unit$Unit)
g_beta_unit$RecYear=as.factor(g_beta_unit$RecYear)
g_beta_m<-lmer(log(distance)~FireGrzTrt*RecYear+(1|Unit), data=g_beta_unit)
check_model(g_beta_m)
anova(g_beta_m)
testInteractions(g_beta_m, pairwise = "FireGrzTrt")
#####visual####
g_beta_viz<-interactionMeans(g_beta_m)%>%
  mutate(beta =exp(`adjusted mean`),
         b_up=exp(`adjusted mean`+`SE of link`),
         b_low=exp(`adjusted mean`-`SE of link`))
g_beta_fig<-ggplot(g_beta_viz,aes(RecYear, beta,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(RecYear)))+
  geom_errorbar(aes(ymin=b_low,
                    ymax=b_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper betadiversity"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#summarize across years
g_beta_viz_avg<-g_beta_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(beta_avg=mean(beta, na.rm=T),
            b_upp=mean(b_up, na.rm=T),
            b_lower=mean(b_low, na.rm=T))
g_beta_avg_fig<-ggplot(g_beta_viz_avg,aes(FireGrzTrt, beta_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=b_lower,
                    ymax=b_upp),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Grasshopper betadiversity"))+
  xlab(label=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

####G NMDS####
g_t_fire_comp_s<-grassh_tcomm_df%>%
  mutate(year_watershed=paste(RecYear, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  filter(Unit=="south")%>%
  ungroup()%>%
  select(-gcount, -total_count, -rep_id,-wsd_rep,-year_watershed)%>%
  distinct()%>%
  pivot_wider(names_from = Species, values_from = rel_abund, values_fill = 0)
g_t_fire_comp_n<-grassh_tcomm_df%>%
  mutate(year_watershed=paste(RecYear, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  filter(Unit=="north")%>%
  ungroup()%>%
  select(-gcount, -total_count, -rep_id,-wsd_rep,-year_watershed)%>%
  distinct()%>%
  pivot_wider(names_from = Species, values_from = rel_abund, values_fill = 0)
#split environmental and species data
#south
g_burn_time_sp_data_s <- g_t_fire_comp_s %>%
  ungroup()%>%
  dplyr::select(-1:-6)
g_burn_time_env_data_s <- g_t_fire_comp_s%>%dplyr::select(1:6)
#north
g_burn_time_sp_data_n <- g_t_fire_comp_n %>%
  ungroup()%>%
  dplyr::select(-1:-6)
g_burn_time_env_data_n <- g_t_fire_comp_n%>%dplyr::select(1:6)

#get nmds1 and 2
g_burn_time_mds_s <- metaMDS(g_burn_time_sp_data_s, distance = "bray",k=2) 
g_burn_time_mds_n <- metaMDS(g_burn_time_sp_data_n, distance = "bray",k=2) 
#combine NMDS1 and 2 with factor columns and create centroids
g_burn_time_mds_scores_s <- data.frame(g_burn_time_env_data_s, scores(g_burn_time_mds_s, display="sites"))%>%
  group_by(RecYear, Unit,time_fire,Watershed)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))
g_burn_time_mds_scores_n <- data.frame(g_burn_time_env_data_n, scores(g_burn_time_mds_n, display="sites"))%>%
  group_by(RecYear, Unit,time_fire,Watershed)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))

#####visual-plotting centroid through time####
g_t_s<-ggplot(g_burn_time_mds_scores_s, aes(x=NMDS1_mean, y=NMDS2_mean, fill=time_fire, shape=Watershed))+
  geom_point(size=8, stroke=2)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  scale_shape_manual(values=c(21:24))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
g_t_n<-ggplot(g_burn_time_mds_scores_n, aes(x=NMDS1_mean, y=NMDS2_mean, fill=time_fire, shape=Watershed))+
  geom_point(size=8, stroke=2)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  scale_shape_manual(values=c(21:24))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

#####stat #####
#compare by timesince fire include watershed and year
g_dist_south<-vegdist(g_burn_time_sp_data_s)
g_dist_north<-vegdist(g_burn_time_sp_data_n)
g_permanova_south<-adonis2(g_dist_south~g_burn_time_env_data_s$time_fire+g_burn_time_env_data_s$Watershed+as.factor(g_burn_time_env_data_s$RecYear), by="terms")
g_permanova_south
#multiple comparisons
g_pairwise_south<-pairwise.adonis2(g_dist_south~time_fire+Watershed+as.factor(RecYear),by="terms", data=g_burn_time_env_data_s)

#comparing time since fire within each watershed to validate
g_wat_south<-vegdist(g_burn_time_sp_data_s[g_burn_time_env_data_s$Watershed=="C03C",])
adonis2(g_wat_south~time_fire+as.factor(RecYear), by="terms", data=g_burn_time_env_data_s[g_burn_time_env_data_s$Watershed=="C03C",])
g_wat_south<-vegdist(g_burn_time_sp_data_s[g_burn_time_env_data_s$Watershed=="C03A",])
adonis2(g_wat_south~time_fire+as.factor(RecYear), by="terms", data=g_burn_time_env_data_s[g_burn_time_env_data_s$Watershed=="C03A",])
g_wat_south<-vegdist(g_burn_time_sp_data_s[g_burn_time_env_data_s$Watershed=="C03B",])
adonis2(g_wat_south~time_fire+as.factor(RecYear), by="terms", data=g_burn_time_env_data_s[g_burn_time_env_data_s$Watershed=="C03B",])

g_permanova_north<-adonis2(g_dist_north~g_burn_time_env_data_n$time_fire+g_burn_time_env_data_n$Watershed+as.factor(g_burn_time_env_data_n$RecYear), by="terms")
g_permanova_north
#multiple comparisons
pairwise.adonis2(g_dist_north~time_fire+Watershed+as.factor(RecYear),by="terms", data=g_burn_time_env_data_n)
g_wat_north<-vegdist(g_burn_time_sp_data_n[g_burn_time_env_data_n$Watershed=="C3SA",])
adonis2(g_wat_north~time_fire+as.factor(RecYear), by="terms", data=g_burn_time_env_data_n[g_burn_time_env_data_n$Watershed=="C3SA",])
g_wat_north<-vegdist(g_burn_time_sp_data_n[g_burn_time_env_data_n$Watershed=="C3SB",])
adonis2(g_wat_north~time_fire+as.factor(RecYear), by="terms", data=g_burn_time_env_data_n[g_burn_time_env_data_n$Watershed=="C3SB",])
g_wat_north<-vegdist(g_burn_time_sp_data_n[g_burn_time_env_data_n$Watershed=="C3SC",])
adonis2(g_wat_north~time_fire+as.factor(RecYear), by="terms", data=g_burn_time_env_data_n[g_burn_time_env_data_n$Watershed=="C3SC",])

####lifeform####
#wrangle data
g_comp_past_life<-g_comp_wide%>%
  pivot_longer(8:60, names_to = "Species", values_to = "abundance")%>%
  group_by(Unit,RecYear,FireGrzTrt,Species)%>%
  summarise(abundance=mean(abundance, na.rm=T))%>%
  left_join(feeding_df, by="Species")%>%
  #removing unknown
  filter(feeding!="unknown")
g_comp_life<-g_comp_past_life%>%
  group_by(Unit,RecYear,FireGrzTrt,feeding)%>%
  summarise(abund=sum(abundance))%>%
  group_by(Unit,RecYear,FireGrzTrt)%>%
  mutate(t_abund=sum(abund),
         rel_abund=abund/t_abund)
g_comp_life$Unit<-as.factor(g_comp_life$Unit)
g_comp_life$RecYear<-as.factor(g_comp_life$RecYear)
g_comp_life$FireGrzTrt<-as.factor(g_comp_life$FireGrzTrt)
#####analysis####
#perennial grass
gh_grass<-lmer(rel_abund~FireGrzTrt*RecYear+(1|Unit), data=g_comp_life[g_comp_life$feeding=="grass",])
check_model(gh_grass)
anova(gh_grass)
testInteractions(gh_grass,pairwise="FireGrzTrt")
gh_grass_data<-interactionMeans(gh_grass)%>%
  mutate(group="Grass feeders")
#annual grass
gh_forb<-lmer(rel_abund~FireGrzTrt*RecYear+(1|Unit), data=g_comp_life[g_comp_life$feeding=="forb",])
check_model(gh_forb)
anova(gh_forb)
testInteractions(gh_forb,pairwise="FireGrzTrt",fixed="RecYear")
gh_forb_data<-interactionMeans(gh_forb)%>%
  mutate(group="Forb feeders")

#perennial forbish
gh_mix<-lmer(rel_abund~FireGrzTrt*RecYear+(1|Unit), data=g_comp_life[g_comp_life$feeding=="mix",])
check_model(gh_mix)
anova(gh_mix)
gh_mix_data<-interactionMeans(gh_mix)%>%
  mutate(group="Mix feeders")

#####visual####
g_group_data<-gh_mix_data%>%
  bind_rows(gh_grass_data,gh_forb_data)%>%
  mutate(grp_mean=`adjusted mean`,
           upper=(`adjusted mean`+`SE of link`),
         lower=(`adjusted mean`-`SE of link`))
#explore
# ggplot(g_comp_life[g_comp_life$feeding=="mix",],aes(RecYear, rel_abund, fill=FireGrzTrt))+
#   geom_boxplot()
g_grp_data_ready<-g_group_data%>%
  group_by(group, FireGrzTrt)%>%
  summarise(grp_m=mean(grp_mean, na.rm=T),
            grp_upper=mean(upper, na.rm=T),
            grp_lower=mean(lower, na.rm=T))
g_group<-ggplot(g_grp_data_ready,aes(group, grp_m,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=grp_lower,
                    ymax=grp_upper),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Relative cover")+
  xlab(label="Grasshopper group")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

#BIRD####
##data wrangling####


#merge data
PBG_bird_merge<-PBG_bird_data_raw%>%
  rename(SpeciesCode=AOUCODE)%>%
  full_join(PBG_bird_C1A, by=c("Year","TransectName","Duration","Date","Observer","StartTime",
                               "Watershed", "DistanceFromLine_Old","SpeciesCode","Count","SEX",
                               "ObsLocation","Notes","Angle","AngularDistance"))%>%
  full_join(PBG_bird_updated, by=c("Year","TransectName","Duration","Date","Observer","StartTime",
                               "Watershed", "DistanceFromLine_Old","SpeciesCode","Count","SEX",
                               "ObsLocation","Notes","Angle","AngularDistance"))%>%
  filter(Year%in%2013:2023)%>%
  filter(Watershed!="C3SC")%>%
  filter(Watershed!="C3SB")%>%
  filter(Watershed!="C3SA")
unique(PBG_bird_merge$Watershed)
#Create a watershed key column to merge with the raw data
watershed_bird_key <- data.frame(Watershed = levels(factor(PBG_bird_merge$Watershed)),# create key to merge with raw data
                                 FireGrzTrt=c("ABG", "PBG", "PBG", "PBG"))
#need to create transect key for easier statistical analysis
transect_key<- data.frame(TransectName = levels(factor(PBG_bird_merge$TransectName)),# create key to merge with raw data
                          Transect=c("A","B","A","B","A","B","A","B"))
#merge key
PBG_bird_viz<-PBG_bird_merge%>%
  left_join(watershed_bird_key, by="Watershed")%>%
  left_join(transect_key, by="TransectName")
#check data for entry errors, 2018 C3A seems to be suspect: 
#I think it has an extra transect that should belong to C3B-90% sure-Corrected in raw data
#email Alice about it

#There are multiple observation by multiple observers 2011-2016
#or just multiple observation by one observer 2011, C3C, C1A
#solutions:group by time, species and take max count, regardless of observer
check_data<-PBG_bird_viz%>%
  filter(Year==2018)%>%
  filter(Watershed=="C3A")
#check unique species code name
codes_name<-data.frame(unique(PBG_bird_viz$SpeciesCode))
#convert to uppercase
PBG_bird_viz$SpeciesCode<-toupper(PBG_bird_viz$SpeciesCode)

#data ready for visuals
PBG_bird_viz_ready_all<-PBG_bird_viz%>%
  rename(total=Count)%>%
  #remove NA otherwise it will appear in max count
  filter(total!="NA")%>%
  #add up each species regardless of sex
  group_by(Year, FireGrzTrt, Watershed, TransectName, Transect, SpeciesCode, Date, Observer)%>%
  summarise(total=sum(total))%>%
  group_by(Year, FireGrzTrt, Watershed, TransectName, Transect, SpeciesCode)%>%
  summarise(total_max=max(total))%>%
  filter(Year!=2019)%>%#removing 2019 because it is missing for C1A
  #rename watershed for ease
  mutate(Watershed=case_when(Watershed=="C1A"~"C01A",
                             Watershed=="C3A"~"C03A",
                             Watershed=="C3B"~"C03B",
                             Watershed=="C3C"~"C03C"))
##Temporal analysis at the local/transect scale####
#prepare data
b_temp_countdata<-PBG_bird_viz_ready_all%>%
  group_by(Year,Watershed, FireGrzTrt, Transect)%>%
  summarise(Tcount=sum(total_max, na.rm=T))%>%
  group_by(Watershed, FireGrzTrt, Transect)%>%
  summarise(temp_count=mean(Tcount, na.rm=T),
            temp_count_sd=sd(Tcount),
            temp_count_cv=temp_count_sd/temp_count)
#convert charcter to factors
b_temp_countdata$Watershed=as.factor(b_temp_countdata$Watershed)  
b_temp_countdata$FireGrzTrt=as.factor(b_temp_countdata$FireGrzTrt)
#start analysis
#temporal mean
b_temp_mean_m<-lmer(temp_count~FireGrzTrt+(1|Watershed), data=b_temp_countdata)
check_model(b_temp_mean_m)
anova(b_temp_mean_m)
#temporal variability
b_temp_sd_m<-lmer(temp_count_sd~FireGrzTrt+(1|Watershed), data=b_temp_countdata)
check_model(b_temp_sd_m)
check_outliers(b_temp_sd_m)
anova(b_temp_sd_m)


b_temp_cv_m<-lmer(temp_count_cv~FireGrzTrt+(1|Watershed), data=b_temp_countdata)
check_model(b_temp_cv_m)
check_normality(b_temp_cv_m)
anova(b_temp_cv_m)

##visuals####
b_TBM<-interactionMeans(b_temp_mean_m)%>%
  mutate(resp="count_mean")
b_TBSD<-interactionMeans(b_temp_sd_m)%>%
  mutate(resp="count_sd")
b_TBCV<-interactionMeans(b_temp_cv_m)%>%
  mutate(resp="count_cv")
b_temp_count_vizdata<-b_TBM%>%
  bind_rows(b_TBSD,b_TBCV)%>%
  mutate(avg=`adjusted mean`,
         up=`adjusted mean`+`SE of link`,
         low=`adjusted mean`-`SE of link`)

b_TBM_viz<-ggplot(b_temp_count_vizdata[b_temp_count_vizdata$resp=="count_mean",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Bird abundance"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
b_TBSD_viz<-ggplot(b_temp_count_vizdata[b_temp_count_vizdata$resp=="count_sd",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Bird SD"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

b_TBCV_viz<-ggplot(b_temp_count_vizdata[b_temp_count_vizdata$resp=="count_cv",],aes(FireGrzTrt, avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=low,
                    ymax=up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Bird CV"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
###combine figures####
bird_temp_viz<-b_TBM_viz/b_TBSD_viz+b_TBCV_viz+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")

##Spatial variability at pasture scale####
b_count_tfire<-PBG_bird_viz_ready_all%>%
  group_by(Year,Watershed, FireGrzTrt, Transect)%>%
  summarise(Tcount=sum(total_max, na.rm=T))%>%
  mutate(year_watershed=paste(Year, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")

###time_fire analysis to show differences in PBG patches####
b_count_tfire$Year=as.factor(b_count_tfire$Year)
b_count_tfire$time_fire=as.factor(b_count_tfire$time_fire)
b_count_tfire$Watershed=as.factor(b_count_tfire$Watershed)
#model
b_t_count_m<-lmer(log(Tcount)~time_fire+Year+(1|Watershed), data=b_count_tfire)#since spatial sd and cv below cannot produce interaction
check_model(b_t_count_m)
anova(b_t_count_m)
testInteractions(b_t_count_m,pairwise = "time_fire", fixed="Year")
testInteractions(b_t_count_m,pairwise = "time_fire")

b_t_count_mean<-interactionMeans(b_t_count_m)%>%
  mutate(count_m=exp(`adjusted mean`),
         count_up=exp(`adjusted mean`+`SE of link`),
         count_low=exp(`adjusted mean`-`SE of link`))

b_t_fire_count_fig<-ggplot(b_t_count_mean,aes(Year, count_m,col=time_fire))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(Year)))+
  geom_errorbar(aes(ymin=count_low,
                    ymax=count_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label=expression("Bird abundance"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
b_t_c_mean_avg<-b_t_count_mean%>%
  group_by(time_fire)%>%
  summarise(count_avg=mean(count_m, na.rm=T),
            count_se=SE_function(count_m))

b_t_fire_c_avg_fig<-ggplot(b_t_c_mean_avg,aes(time_fire, count_avg,col=time_fire))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=count_avg-count_se,
                    ymax=count_avg+count_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label=expression("Bird abundance"))+
  xlab(label="Year since last fire")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

#wrangle data
b_spat_cdata<-b_count_tfire%>%
  group_by(Year, FireGrzTrt)%>%
  summarise(spat_c=mean(Tcount, na.rm=T),
            spat_c_sd=sd(Tcount),
            spat_c_cv=spat_c_sd/spat_c)
#convert character to factors
b_spat_cdata$FireGrzTrt=as.factor(b_spat_cdata$FireGrzTrt)
b_spat_cdata$Year=as.factor(b_spat_cdata$Year)
#analysis
b_spat_c_m<-lm(log(spat_c)~FireGrzTrt+Year, data=b_spat_cdata)#nan returned when interaction included
check_model(b_spat_c_m)
anova(b_spat_c_m)
summary(b_spat_c_m)
testInteractions(b_spat_c_m, pairwise = "FireGrzTrt")
#SD
b_spat_c_sd<-lm(spat_c_sd~FireGrzTrt+Year, data=b_spat_cdata)
check_model(b_spat_c_sd)
anova(b_spat_c_sd)
#cv
b_spat_c_cv<-lm(spat_c_cv~FireGrzTrt+Year, data=b_spat_cdata)
check_model(b_spat_c_cv)
anova(b_spat_c_cv)
testInteractions(b_spat_c_cv,pairwise = "FireGrzTrt")
##visuals####
bcount_mean<-interactionMeans(b_spat_c_m)%>%
  mutate(count_m=exp(`adjusted mean`),
         c_up=exp(`adjusted mean`+`std. error`),
         c_low=exp(`adjusted mean`-`std. error`),
         variab="AVG")

bcount_sd<-interactionMeans(b_spat_c_sd)%>%
  mutate(count_m=(`adjusted mean`),
         c_up=(`adjusted mean`+`std. error`),
         c_low=(`adjusted mean`+`std. error`),
         variab="SD")

bcount_cv<-interactionMeans(b_spat_c_cv)%>%
  mutate(count_m=(`adjusted mean`),
         c_up=(`adjusted mean`+`std. error`),
         c_low=(`adjusted mean`+`std. error`),
         variab="CV")
###average figures####
#wrangle data
combo_spat_bcount<-bcount_mean%>%
  bind_rows(bcount_sd,bcount_cv)%>%
  group_by(FireGrzTrt, variab)%>%
  summarise(bcount_avg=mean(count_m, na.rm=T),
            c_se=SE_function(count_m))
#visual
b_avg_spat_viz<-ggplot(combo_spat_bcount[combo_spat_bcount$variab=="AVG",],aes(FireGrzTrt, bcount_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=bcount_avg-c_se,
                    ymax=bcount_avg+c_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Bird abundance"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

b_sd_spat_viz<-ggplot(combo_spat_bcount[combo_spat_bcount$variab=="SD",],aes(FireGrzTrt, bcount_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=bcount_avg-c_se,
                    ymax=bcount_avg+c_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Bird SD"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

b_cv_spat_viz<-ggplot(combo_spat_bcount[combo_spat_bcount$variab=="CV",],aes(FireGrzTrt, bcount_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=bcount_avg-c_se,
                    ymax=bcount_avg+c_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Bird CV"))+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
###combine####
bird_spatial_viz<-b_t_fire_c_avg_fig/b_sd_spat_viz+b_cv_spat_viz+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")

#Bird composition#####
b_tcomm_df<-PBG_bird_viz_ready_all%>%
  group_by(FireGrzTrt,Watershed, Year,Transect)%>%
  mutate(total_count=sum(total_max))%>%
  group_by(FireGrzTrt, Watershed, Year,Transect,SpeciesCode)%>%
  mutate(rel_abund=total_max/total_count)%>%
  mutate(rep_id=paste(Watershed, FireGrzTrt, Transect, sep="_"),
         wsd_rep=paste(Watershed,Transect, sep="_"))
##community change calculation yearly in each transect####
b_time_change<-multivariate_change(b_tcomm_df, species.var="SpeciesCode",
                                   abundance.var = "rel_abund",
                                   replicate.var="rep_id",
                                   time.var = "Year",
                                   treatment.var="wsd_rep")

#combine to have datatset with tretament variables for analysis
b_t_change<-b_tcomm_df%>%
  ungroup()%>%
  select(Watershed, Transect, FireGrzTrt, wsd_rep)%>%
  distinct()%>%
  left_join(b_time_change, by="wsd_rep")%>%
  mutate(year_diff=paste(Year, Year2, sep="_"))

b_t_change$FireGrzTrt<-as.factor(b_t_change$FireGrzTrt)
b_t_change$Watershed<-as.factor(b_t_change$Watershed)
b_t_change$year_diff<-as.factor(b_t_change$year_diff)

##model for analysis####
b_t_ch_m<-lmer(composition_change~FireGrzTrt*year_diff+(1|Watershed), data=b_t_change)
check_model(b_t_ch_m)
anova(b_t_ch_m)

##visuals####
#wrangle data from model
b_comp_change_mean<-interactionMeans(b_t_ch_m)%>%
  mutate(comp_chang=(`adjusted mean`),
         chang_up=(`adjusted mean`+`SE of link`),
         chang_low=(`adjusted mean`-`SE of link`))

b_chang_mean_fig<-ggplot(b_comp_change_mean,aes(year_diff, comp_chang,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(year_diff)))+
  geom_errorbar(aes(ymin=chang_low,
                    ymax=chang_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Bird composition change"))+
  xlab(label="Year comparison")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#calculate average 
b_chang_avg<-b_comp_change_mean%>%
  group_by(FireGrzTrt)%>%
  summarise(comp_chang_avg=mean(comp_chang, na.rm=T),
            c_chang_se=SE_function(comp_chang))
#visual
b_avg_chang_viz<-ggplot(b_chang_avg,aes(FireGrzTrt, comp_chang_avg,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=comp_chang_avg-c_chang_se,
                    ymax=comp_chang_avg+c_chang_se),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Bird composition change")+
  xlab(labe=NULL)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

##Spatial bird comp####
###transect scale diversity####
b_rich_trsct <- community_structure(b_tcomm_df, time.var = "Year", 
                                    abundance.var = "rel_abund",
                                    replicate.var = "rep_id", metric = "Evar")
#combine with treatmnet info
b_rich_trsct<-b_tcomm_df%>%
  ungroup()%>%
  select(Year, Watershed, Transect, FireGrzTrt, rep_id)%>%
  distinct()%>%
  left_join(b_rich_trsct, by=c("Year","rep_id"))

####analysis####
#convert to factors for analysis
b_rich_trsct$Year<-as.factor(b_rich_trsct$Year)
b_rich_trsct$Watershed<-as.factor(b_rich_trsct$Watershed)
b_rich_trsct$FireGrzTrt<-as.factor(b_rich_trsct$FireGrzTrt)

b_rich_trsct_m<-lmer(richness~FireGrzTrt*Year+(1|Watershed), data=b_rich_trsct)
check_model(b_rich_trsct_m)
anova(b_rich_trsct_m)


b_even_trsct_m<-lmer(Evar~FireGrzTrt*Year+(1|Watershed), data=b_rich_trsct)
check_model(b_even_trsct_m)
anova(b_even_trsct_m)


#####visuals#####
#wrangle data from model
b_rich_t<-interactionMeans(b_rich_trsct_m)%>%
  mutate(rich =`adjusted mean`,
         r_up=`adjusted mean`+`SE of link`,
         r_low=`adjusted mean`-`SE of link`)

b_rich_trsct_fig<-ggplot(b_rich_t,aes(Year, rich,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path(aes(as.numeric(Year)))+
  geom_errorbar(aes(ymin=r_low,
                    ymax=r_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label=expression("Bird richness"))+
  xlab(label="Year")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
####time_fire transect####
b_tcomm_df$Year<-as.factor(b_tcomm_df$Year)
b_rich_t_fire_tst<-b_tcomm_df%>%
  ungroup()%>%
  mutate(year_watershed=paste(Year, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  select(Year,Watershed, Transect, time_fire, rep_id)%>%
  distinct()%>%
  left_join(b_rich_trsct, by=c("Year","rep_id","Watershed","Transect"))

####analysis####
#convert to factors for analysis
b_rich_t_fire_tst$Watershed<-as.factor(b_rich_t_fire_tst$Watershed)
b_rich_t_fire_tst$time_fire<-as.factor(b_rich_t_fire_tst$time_fire)

b_rich_t_fire_tst_m<-lmer(richness~time_fire*Year+(1|Watershed), data=b_rich_t_fire_tst)
check_model(b_rich_t_fire_tst_m)
anova(b_rich_t_fire_tst_m)

b_even_t_fire_tst_m<-lmer(Evar~time_fire*Year+(1|Watershed), data=b_rich_t_fire_tst)
check_model(b_even_t_fire_tst_m)
anova(b_even_t_fire_tst_m)
ggplot(b_rich_t_fire_tst, aes(Year, richness, col=time_fire))+
  geom_boxplot()
ggplot(b_rich_t_fire_tst, aes(Year, Evar, col=time_fire))+
  geom_boxplot()
####visual####
b_rich_tsf<-interactionMeans(b_rich_t_fire_tst_m)%>%
  mutate(rich =`adjusted mean`,
         r_up=`adjusted mean`+`SE of link`,
         r_low=`adjusted mean`-`SE of link`)%>%
  group_by(time_fire)%>%
  summarise(rich=mean(rich),
            r_up=mean(r_up),
            r_low=mean(r_low))
b_even_tsf<-interactionMeans(b_even_t_fire_tst_m)%>%
  mutate(rich =`adjusted mean`,
         r_up=`adjusted mean`+`SE of link`,
         r_low=`adjusted mean`-`SE of link`)%>%
  group_by(time_fire)%>%
  summarise(rich=mean(rich),
            r_up=mean(r_up),
            r_low=mean(r_low))

b_t_fire_rich_avg_fig<-ggplot(b_rich_tsf,aes(time_fire, rich,col=time_fire))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=r_low,
                    ymax=r_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label=expression("Bird richness"))+
  xlab(label="Year since last fire")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
b_t_fire_even_avg_fig<-ggplot(b_even_tsf,aes(time_fire, rich,col=time_fire))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=r_low,
                    ymax=r_up),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label=expression("Bird evenness"))+
  xlab(label="Year since last fire")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
####combine####
b_div_timefire<-b_t_fire_rich_avg_fig/b_t_fire_even_avg_fig+ plot_layout(guides = "collect")&plot_annotation(tag_levels = "A")&theme(legend.position = "none")

###Bird NMDS####
b_t_fire_comp<-b_tcomm_df%>%
  mutate(year_watershed=paste(Year, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  ungroup()%>%
  select(-total_max, -total_count, -rep_id,-wsd_rep,-year_watershed)%>%
  distinct()%>%
  pivot_wider(names_from = SpeciesCode, values_from = rel_abund, values_fill = 0)
#split environmental and species data
b_burn_time_sp_data <- b_t_fire_comp %>%
  ungroup()%>%
  dplyr::select(-1:-6)
b_burn_time_env_data <- b_t_fire_comp%>%dplyr::select(1:6)

#get nmds1 and 2
b_burn_time_mds <- metaMDS(b_burn_time_sp_data, distance = "bray",k=2) 

#combine NMDS1 and 2 with factor columns and create centroids
b_burn_time_mds_scores<- data.frame(b_burn_time_env_data, scores(b_burn_time_mds, display="sites"))%>%
  group_by(Year,time_fire,Watershed)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))

#####visual-plotting centroid through time####
b_t_s<-ggplot(b_burn_time_mds_scores, aes(x=NMDS1_mean, y=NMDS2_mean, fill=time_fire, shape=Watershed))+
  geom_point(size=8, stroke=2)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  scale_shape_manual(values=c(21:24))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )

#compare by timesince fire include watershed and year
b_dist<-vegdist(b_burn_time_sp_data)
b_permanova<-adonis2(b_dist~b_burn_time_env_data$time_fire+b_burn_time_env_data$Watershed+as.factor(b_burn_time_env_data$Year), by="terms")
b_permanova
#multiple comparisons
b_pairwise<-pairwise.adonis2(b_dist~time_fire+Watershed+as.factor(Year),by="terms", data=b_burn_time_env_data)
b_pairwise
#only PBG0 vs ABG0 different also ABGO vs PBG2 

####groups based on grassland obligate####
bird_group<-b_tcomm_df%>%
  select(-total_count:-wsd_rep)%>%
  select(Year, FireGrzTrt, Watershed, Transect, SpeciesCode,total_max)%>%
  pivot_wider(names_from = "SpeciesCode", values_from = "total_max", values_fill = 0)%>%
  pivot_longer(5:56, names_to = "SpeciesCode", values_to = "abundance")%>%
  group_by(Year, FireGrzTrt, Watershed, Transect, SpeciesCode)%>%
  summarise(abundance=mean(abundance, na.rm=T))%>%
  left_join(PBG_bird_class, by="SpeciesCode")
bird_group$Grassland=as.factor(bird_group$Grassland)
bird_grp_data<-bird_group%>%
  group_by(Year,FireGrzTrt, Watershed, Transect, Grassland)%>%
  summarise(abund=sum(abundance))%>%
  group_by(Year,FireGrzTrt, Watershed,Transect)%>%
  mutate(t_abund=sum(abund),
         rel_abund=abund/t_abund)
bird_grp_data$Year<-as.factor(bird_grp_data$Year)
bird_grp_data$Watershed<-as.factor(bird_grp_data$Watershed)
bird_grp_data$FireGrzTrt<-as.factor(bird_grp_data$FireGrzTrt)
#####analysis####
#obligate grassland species
bird_grass<-lmer(rel_abund~FireGrzTrt*Year+(1|Watershed), data=bird_grp_data[bird_grp_data$Grassland=="TRUE",])
check_model(bird_grass)
check_normality(bird_grass)
anova(bird_grass)
testInteractions(bird_grass,pairwise="FireGrzTrt")
b_grass_data<-interactionMeans(bird_grass)%>%
  mutate(group="Grassland obligate")

#non-grassland obligate species
bird_non_grass<-lmer(rel_abund~FireGrzTrt*Year+(1|Watershed), data=bird_grp_data[bird_grp_data$Grassland=="FALSE",])
check_model(bird_non_grass)
anova(bird_non_grass)
testInteractions(bird_non_grass,pairwise="FireGrzTrt")
check_homogeneity(bird_non_grass)
check_normality(bird_non_grass)
b_nongrass_data<-interactionMeans(bird_non_grass)%>%
  mutate(group="Generalist")
#combine for visual
b_group_data<-b_grass_data%>%
  bind_rows(b_nongrass_data)%>%
  mutate(grp_mean=`adjusted mean`,
         upper=(`adjusted mean`+`SE of link`),
         lower=(`adjusted mean`-`SE of link`))

b_grp_data_ready<-b_group_data%>%
  group_by(group, FireGrzTrt)%>%
  summarise(grp_m=mean(grp_mean, na.rm=T),
            grp_upper=mean(upper, na.rm=T),
            grp_lower=mean(lower, na.rm=T))
b_group<-ggplot(b_grp_data_ready,aes(group, grp_m,col=FireGrzTrt))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=grp_lower,
                    ymax=grp_upper),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#009E73"))+
  ylab(label="Relative cover")+
  xlab(label="Bird group")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#####time since fire####
bird_grp_data_tf<-bird_group%>%
  mutate(year_watershed=paste(Year, Watershed, sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  group_by(Year,FireGrzTrt, Watershed, time_fire, Transect, Grassland)%>%
  summarise(abund=sum(abundance))%>%
  group_by(Year,FireGrzTrt, Watershed,Transect)%>%
  mutate(t_abund=sum(abund),
         rel_abund=abund/t_abund)
bird_grp_data_tf$Year<-as.factor(bird_grp_data_tf$Year)
bird_grp_data_tf$Watershed<-as.factor(bird_grp_data_tf$Watershed)
bird_grp_data_tf$time_fire<-as.factor(bird_grp_data_tf$time_fire)

#####analysis####
#obligate grassland species
bird_grass_tf<-lmer(rel_abund~time_fire*Year+(1|Transect), data=bird_grp_data_tf[bird_grp_data_tf$Grassland=="TRUE",])
check_model(bird_grass_tf)
check_normality(bird_grass_tf)
anova(bird_grass_tf)
testInteractions(bird_grass_tf,pairwise="time_fire")
b_grass_tf_data<-interactionMeans(bird_grass_tf)%>%
  mutate(group="Grassland obligate")

#non-grassland obligate species
bird_non_grass_tf<-lmer(rel_abund~time_fire*Year+(1|Transect), data=bird_grp_data_tf[bird_grp_data_tf$Grassland=="FALSE",])
check_model(bird_non_grass_tf)
anova(bird_non_grass_tf)
testInteractions(bird_non_grass_tf,pairwise="time_fire")

b_nongrass_tf_data<-interactionMeans(bird_non_grass_tf)%>%
  mutate(group="Generalist")
#combine for visual
b_group_data_tf<-b_grass_tf_data%>%
  bind_rows(b_nongrass_tf_data)%>%
  mutate(grp_mean=`adjusted mean`,
         upper=(`adjusted mean`+`SE of link`),
         lower=(`adjusted mean`-`SE of link`))

b_grp_data_tf_ready<-b_group_data_tf%>%
  group_by(group, time_fire)%>%
  summarise(grp_m=mean(grp_mean, na.rm=T),
            grp_upper=mean(upper, na.rm=T),
            grp_lower=mean(lower, na.rm=T))
b_group_tf<-ggplot(b_grp_data_tf_ready,aes(group, grp_m,col=time_fire))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=grp_lower,
                    ymax=grp_upper),width=0.0125)+
  scale_color_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab(label="Relative cover")+
  xlab(label="Bird group")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()   # Remove minor gridlines
  )
#All figures####
#combine the transect temporal change

#combine spatial change

#combine temporal change

#combine diversity and beta
#combine diversity and beta showing each year

#combine NMDS

#combine groups
