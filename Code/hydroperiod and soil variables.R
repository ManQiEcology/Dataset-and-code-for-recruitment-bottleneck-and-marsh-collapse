rm(list=ls())
library(here)
library(fs)
library(readxl)
library(hydrostats)
library(lubridate)
library(Rcpp)
library(ggplot2)
library(dplyr)
library(FSA)  
library(car)
library(DescTools)
library(multcomp)
library(Rmisc)  
library(agricolae)
library(rcompanion)
library(patchwork)
library(rstatix) # For Kruskal-Wallis and post-hoc tests

setwd(path_dir(here::here()))



##########################################################################################
#HYDROPERIOD CALCULATION: THIS CODE IS MODIFIED FROM CODE OF GROUNDWATER PAPER
##########################################################################################

Sys.setlocale("LC_TIME", "English") #set the default system language to be English
#1.load well info
well_ele<-read.csv("Data/Groundwater/Well info.txt",sep="", header = T)

#1.1 load groundwater data of Deal Island
DI_TC<-read_excel("Data/Groundwater/20138781_DI_G_TC.xlsx")
DI_H<-read_excel("Data/Groundwater/20574343_DI_G_H.xlsx")
DI_E<-read_excel("Data/Groundwater/20574344_DI_G_E.xlsx")
DI_P<-read_excel("Data/Groundwater/20292500_DI_G_P.xlsx")
#assign site name to sensor depth
DI_TC$DI_TC<-DI_TC$`Sensor depth (Meters)`
DI_H$DI_H<-DI_H$`Sensor depth (Meters)`
DI_E$DI_E<-DI_E$`Sensor depth (Meters)`
DI_P$DI_P<-DI_P$`Sensor depth (Meters)`
#merge groundwater data from different sites from Deal Island
DI<-NULL
DI<-merge(DI_TC,DI_H,by="Time, GMT-04:00")
DI<-merge(DI,DI_E,by="Time, GMT-04:00")
DI<-merge(DI,DI_P,by="Time, GMT-04:00")
#extract groundwater data from different sites of Deal Island and generate groundwater data file for Deal Island
DI_SD<-data.frame("DateTime"=DI$`Time, GMT-04:00`,TC=DI$DI_TC,H=DI$DI_H,E=DI$DI_E,P=DI$DI_P)
DI_SD$DateTime<-ymd_hms(DI_SD$DateTime,tz="US/Eastern") #define time with lubridate
#DI_SD above contains sensor depth data of Deal Island

#remove original DI dataset
rm(DI_E,DI_H,DI_P,DI_TC,DI)

#1.2 load groundwater of Farm Creek Marsh
FCM_HC<-read_excel("Data/Groundwater/20574349_FCM_G_HC.xlsx")
FCM_DC<-read_excel("Data/Groundwater/20574347_FCM_G_DC.xlsx")
FCM_PC<-read.csv("Data/Groundwater/Site4S-FCM_G_PC.csv",header = T)#the time is in EST
#Sensor depth of healthy and dieback
FCM_HC$FCM_HC<-FCM_HC$`Sensor depth (Meters)`
FCM_DC$FCM_DC<-FCM_DC$`Sensor depth (Meters)`
FCM<-NULL
FCM<-merge(FCM_HC,FCM_DC,by="Time, GMT-04:00")
FCM_SD<-data.frame("DateTime"=FCM$`Time, GMT-04:00`,H=FCM$FCM_HC,E=FCM$FCM_DC) #datetime in EST
FCM_SD$DateTime<-ymd_hms(FCM_SD$DateTime,tz="US/Eastern") #define time with lubridate
FCM_SD<-subset(FCM_SD,DateTime>min(DI_SD$DateTime,na.rm = T)&DateTime<max(DI_SD$DateTime,na.rm = T))
#FCM_SD above contains sensor depth data healthy (H) and dieback (E) from May 2019 to March 2020
#remove orginal dataset from FCM
rm(FCM_DC,FCM_HC,FCM)


#preliminary process of groundwater data from FCM pond
yr2mth<-function(x) {
  for (i in 1:9) {
    x<-gsub(paste("200",i,sep = ""),i,x)
  }
  for (i in 10:12) {
    x<-gsub(paste("20",i,sep = ""),i,x)
  }
  x
}
FCM_PC$DateTime<-yr2mth(FCM_PC$DateTime)
FCM_PC$DateTime<-mdy_hm(FCM_PC$DateTime,tz="US/Eastern")# datetime is US/Eastern '(UTC-5)'
FCM_PC$waterlevel <-FCM_PC$waterlevel*0.3048 # Clevel is the water level relative to MSL(NAVD88),turn foot to meter 1 foot =0.3048*meter, adjust water level (-0.15m) to get it consistent with the healthy and the dieback
FCM_PC<-subset(FCM_PC,DateTime>min(FCM_SD$DateTime,na.rm = T)&DateTime<max(FCM_SD$DateTime,na.rm = T),tz="US/Eastern") #keep pond dataseries within dieback and healthy patch
#############################################################################################
#2. Turn water depth to water level relative to MSL (NAVD1988)
Logger_depth<-well_ele$Logger_depth_2019_.cm./100#logger depth has a unit of cm, transfer to m
surface_elevation<-well_ele$Soil_surface_elevation.m._2019
logger_elevation<-well_ele$Well_top_elevation.m._2020-well_ele$Well_top_to_logger.cm./100 #consider well height and well top to logger doesb't change much, within 2cm from 2019 to 2020, thus use this to calculate water level
#for DI, calckulate water depth relatuve to Mean Sea Level (NAVD1988)
DI_WL_MSL<-data.frame("DateTime"=DI_SD$DateTime,
                      TC=DI_SD$TC-Logger_depth[1]+surface_elevation[1], 
                      H=DI_SD$H+surface_elevation[4]-Logger_depth[4], #marsh elevation might decreased by 4cm from 2019 May to 2020 March
                      E=DI_SD$E+surface_elevation[2]-Logger_depth[2], #marsh surface elevation and well elevation might decreased by 2cm from 2019 May to 2020 Mar
                      P=DI_SD$P+surface_elevation[3]-Logger_depth[3]) # marsh surface elevation and well elevation variation is around 1cm from 2019 May to 2020 March

FCM_WL_MSL<-data.frame("DateTime"=FCM_SD$DateTime,
                       H=FCM_SD$H+surface_elevation[7]-Logger_depth[7], #marsh elevation decrease by 2cm, but surface accreation increase by 2cm
                       E=FCM_SD$E+0.5*(surface_elevation[6]-Logger_depth[6]+logger_elevation[6]))# since marsh surface elevation and well elevation increased by 0.58 cm from 2019 May to 2020 March, we use average logger height

#load bishop tide water,Datum:NAVD88
bishop<-as.data.frame(read.csv("Data/Groundwater/bishop.csv",header = T))
bishop$DateTime<-paste(bishop$Date,bishop$Time..GMT.) #merge data and hm
bishop$DateTime<-ymd_hm(bishop$DateTime,tz="GMT") #define the date format with GMT time zone
bishop$DateTime<-with_tz(bishop$DateTime,tz="US/Eastern") # convert the time zone to US/Eastern
bishop<-data.frame("DateTime"=bishop$DateTime,TH=as.numeric(bishop$Verified..m.))
bishop<-subset(bishop,DateTime>min(FCM_SD$DateTime,na.rm = T)&DateTime<max(FCM_SD$DateTime,na.rm=T))
#correct DI-TC data at low tide with bishop data
bishoplus0.07<-data.frame(DateTime=bishop$DateTime,
                          TH=bishop$TH+0.07) #water level in DI_TC is 0.07higher than that in bishop, so use bishop+0.07 as replacement
DI_WL_MSL$TC[which(DI_WL_MSL$TC<=0)]<-NA #it is found that data is incorrect when water level of DI_TC is less than zero, thus replace this value with NA

for (i in which(is.na(DI_WL_MSL$TC)==T)) { #replace NA data with bishop+0.07
  x<-which(bishoplus0.07$DateTime==DI_WL_MSL$DateTime[i])
  if (identical(x, integer(0))) next
  else (DI_WL_MSL$TC[i]<-bishoplus0.07$TH[x])
}
#interpolate DI-TC 
model <- lm(TC ~ DateTime, data = DI_WL_MSL) 
a<-approx(DI_WL_MSL$DateTime, DI_WL_MSL$TC, xout=DI_WL_MSL$DateTime[which(is.na(DI_WL_MSL$TC)==T)])
DI_WL_MSL$TC[which(is.na(DI_WL_MSL$TC)==T)] <- a$y
#data in healthy and dieback patch of DI is incorrect at high water level
#######################################################################################
#5.calculate and plot soil saturation index SSI=0 means hydroperid, SSI=0.1, 0.2,0.3 means fraction of tidal period over which the soil remains fully saturated at depth of SSI
##########################################################################################
#5.1 calculate on site monitored water level relative to soil surface, above zero-flooded, below zero-air exposed
#for DI, generate DI_WL_SS
SSI<-0.0
DI_WL_SS<-data.frame("DateTime"=DI_WL_MSL$DateTime,
                     TC=DI_WL_MSL$TC-surface_elevation[1]+SSI,
                     H=DI_WL_MSL$H-surface_elevation[4]+SSI,
                     E=DI_WL_MSL$E-surface_elevation[2]+SSI,
                     P=DI_WL_MSL$P-surface_elevation[3]+SSI)
#for FCM, generate FCM_WL_SS_HE and FCM_WL_SS_P
FCM_WL_SS_HE<-data.frame("DateTime"=FCM_WL_MSL$DateTime,
                         H=FCM_WL_MSL$H-surface_elevation[7]+SSI,
                         E=FCM_WL_MSL$E-surface_elevation[6]+SSI)
FCM_WL_SS_P<-data.frame("DateTime"=FCM_PC$DateTime,P=FCM_PC$waterlevel-surface_elevation[9]+SSI )

#5.2 cacculate elevation drived water level relative to soil surface with tidal gauge records and elevation
bishop$DITC<-bishop$TH - surface_elevation[1]+SSI
bishop$DIH<-bishop$TH-surface_elevation[4]+SSI
bishop$DIE<-bishop$TH-surface_elevation[2]+SSI
bishop$DIP<-bishop$TH-surface_elevation[3]+SSI
bishop$FCMHC<-bishop$TH-surface_elevation[7]+SSI
bishop$FCMDC<-bishop$TH-surface_elevation[6]+SSI
bishop$FCMPC<-bishop$TH-surface_elevation[9]+SSI
#5.3. compare hydroperiod                                                                  #
##########################################################################################
hydrpd <-function(x){
  y<-sum(x>0,na.rm=T)/sum(1-is.na(x),na.rm=T)
  return(y)
}
##in situ reading
#DI
DI_hydrpd<-as.matrix(DI_WL_SS[3:5],ncol=3)
DI_hydrpd<-apply(DI_hydrpd,2,hydrpd)
#FCM_HE
FCM_HE_hydrpd<-as.matrix(FCM_WL_SS_HE[2:3])
FCM_HE_hydrpd<-apply(FCM_HE_hydrpd,2,hydrpd)
#FCM_P
FCM_P_hydrpd<-hydrpd(FCM_WL_SS_P[2])

#construct hydroperiod data frame
hydrpd<-data.frame(Hydroperiod=c(DI_hydrpd,FCM_HE_hydrpd,FCM_P_hydrpd),
                   Zone=factor(c("Neighbor patch","Dieback patch","Pond","Neighbor patch","Dieback patch", "Pond"),
                               levels=c("Pond","Dieback patch", "Neighbor patch")),
                   Site=factor(c(rep("Deal Island",3),rep("Farm Creek Marsh",3)),
                               levels=c("Farm Creek Marsh", "Deal Island")))

Juncus<-data.frame(Hydroperiod=NA,
                   Zone=factor(rep("Juncus patch", 2),
                               levels=c("Juncus patch")),
                   Site=factor(c("Deal Island","Farm Creek Marsh"),
                               levels=c("Farm Creek Marsh", "Deal Island")))
hydrpd<-rbind(hydrpd, Juncus)
str(hydrpd)
##########################################################################################
#REDOX
##########################################################################################

DI_salinity<-read_excel("Data/Pond plant transplanting result-2019/DI.xlsx", sheet = "Salinity")
FCM_salinity<- read_excel("Data/Pond plant transplanting result-2019/FCM.xlsx", sheet = "Salinity")
DI_Redox<-read_excel("Data/Pond plant transplanting result-2019/DI.xlsx", sheet = "Redox and Soil strength")
FCM_Redox<-read_excel("Data/Pond plant transplanting result-2019/FCM.xlsx", sheet = "Redox and Soil strength")
DI_salinity$Month<-factor(DI_salinity$Month,levels=unique(DI_salinity$Month))
FCM_salinity<-subset(FCM_salinity,FCM_salinity$Month!="Sep-20")
FCM_salinity$Month<-factor(FCM_salinity$Month,levels=unique(FCM_salinity$Month))

################################################################################
#1.Redox 
################################################################################
#1.1 Prepare the data
Data_DI<-DI_Redox 
#Data_DI<-subset(Data_DI,Zone!="Juncus patch")
Data_DI$Zone[which(Data_DI$Zone=="Healthy")]<-"Neighbor patch"
#Specify the order of factor levels for plots and Dunnett comparison
Data_DI = mutate(Data_DI, Zone = factor(Zone, 
                                        levels= c("Pond", "Dieback","Neighbor patch", "Juncus patch"))) #Specify the order of factor levels for plots and Dunnett comparison
Data_FCM<-FCM_Redox 
#Data_FCM<-subset(Data_FCM,Zone!="Juncus patch")
Data_FCM$Zone[which(Data_FCM$Zone=="Healthy")]<-"Neighbor patch"
Data_FCM = mutate(Data_FCM, Zone = factor(Zone, 
                                          levels= c("Pond", "Dieback","Neighbor patch", "Juncus patch"))) #Specify the order of factor levels for plots and Dunnett comparison
####################################################################################################
stts_display<-function(x) {
  print(summary(x))
  pairwise<-TukeyHSD(x)$`Site:Zone`
  pairwise<-as.data.frame(pairwise)
  pairwise$contrast<-row.names(pairwise)
  pairwise <- pairwise %>%
    tidyr::separate(contrast, into = c("Group1", "Group2"), sep = "-")
  # Extract pairwise comparisons within Farm Creek Marsh
  pairwise_farm_creek <- pairwise %>%
    filter(grepl("Farm Creek Marsh", Group1) & grepl("Farm Creek Marsh", Group2))
  
  # Extract pairwise comparisons within Deal Island
  pairwise_deal_island <- pairwise %>%
    filter(grepl("Deal Island", Group1) & grepl("Deal Island", Group2))
  print(pairwise_farm_creek[,4:5])
  print(pairwise_deal_island[,4:5])
}
#2.Statistic test

Redox_data<-rbind(Data_FCM,Data_DI)
Redox_data$Site<-factor(c(rep("Farm Creek Marsh", nrow(Data_FCM)),
                          rep("Deal Island",nrow(Data_DI))),
                        levels = c("Farm Creek Marsh","Deal Island"))
Redox_data<-subset(Redox_data,is.na(Redox_data$`Redox (mV)`)==0)
anova_model_redox <- aov(`Redox (mV)` ~ Site + Site: Zone,
                         data=Redox_data)
cld(emmeans(anova_model_redox, ~ Zone | Site), Letters = letters, adjust = "tukey")



#3. Assumption test
# Shapiro-Wilk test for normality
residuals_redox <- residuals(anova_model_redox)
shapiro_test <- shapiro.test(residuals_redox)
print(shapiro_test)

# Levene’s Test (preferred, more robust to non-normality)
library(car)
leveneTest(`Redox (mV)` ~ Site * Zone, data = Redox_data)

#log transform of redox
Redox_data$log_redox <- log(Redox_data$`Redox (mV)` + abs(min(Redox_data$`Redox (mV)`)) + 1)
anova_model_redox_log <- aov(log_redox ~ Site * Zone, data = Redox_data)

# Recheck assumptions
shapiro.test(residuals(anova_model_redox_log))
leveneTest(log_redox ~ Site * Zone, data = Redox_data)

library(ARTool)
# Run ART model
art_model_redox <- art(`Redox (mV)` ~ Site * Zone, data = Redox_data)

# Get ANOVA table
anova(art_model_redox)


# Post-hoc for main effects
art_site <- artlm(art_model_redox, "Site")
art_zone <- artlm(art_model_redox, "Zone")
art_interaction <- artlm(art_model_redox, "Site:Zone")

# Tukey-adjusted pairwise comparisons
pairs(emmeans(art_interaction, ~ Zone | Site), adjust = "tukey")

# Get letters for zones at each site
cld(emmeans(art_interaction, ~ Zone | Site), Letters = letters, adjust = "tukey")


#4. Graphing the Redox

Plot_Data_DI = summarySE(data=Data_DI,
                         "Redox (mV)",
                         groupvars="Zone",
                         conf.interval = 0.95,na.rm =T)
Plot_Data_FCM = summarySE(data=Data_FCM,
                  "Redox (mV)",
                  groupvars="Zone",
                  conf.interval = 0.95,na.rm =T)

Redox<-rbind(Plot_Data_DI,Plot_Data_FCM)
Redox$Site<-factor(c(rep("Deal Island",4),rep("Farm Creek Marsh",4)),
                   levels=c("Farm Creek Marsh", "Deal Island"))

################################################################################
#SOIL STRENGTH
################################################################################
#1 prepare the data

library(dplyr)
#Data<-subset(DI_Redox, Zone!="Pond") 
Data_DI = mutate(Data_DI, `Soil_depth (cm)` = factor(`Soil_depth (cm)`, 
                                          levels=c("20 cm","50 cm")))
Data_FCM = mutate(Data_FCM, `Soil_depth (cm)` = factor(`Soil_depth (cm)`, 
                                                     levels=c("20 cm","50 cm")))
Data_DI_20<-subset(Data_DI,`Soil_depth (cm)`=="20 cm")
Data_FCM_20<-subset(Data_FCM,`Soil_depth (cm)`=="20 cm")
Data_DI_50<-subset(Data_DI,`Soil_depth (cm)`="50 cm")
Data_FCM_50<-subset(Data_FCM,`Soil_depth (cm)`="50 cm")

SS_20_data<-rbind(Data_FCM_20, Data_DI_20)
SS_20_data$Site=factor(c(rep("Farm Creek Marsh",nrow(Data_FCM_20)),
                          rep("Deal Island", nrow(Data_DI_20))),
                          levels=c("Farm Creek Marsh", "Deal Island"))
SS_50_data<-rbind(Data_FCM_50,Data_DI_50)
SS_50_data$Site=factor(c(rep("Farm Creek Marsh",nrow(Data_FCM_50)),
                          rep("Deal Island", nrow(Data_DI_50))),
                        levels=c("Farm Creek Marsh", "Deal Island"))
SS_20_data<-subset(SS_20_data,is.na(`Soil strength (kPa)`)==0)
SS_50_data<-subset(SS_50_data,is.na(`Soil strength (kPa)`)==0)
#2 statistical test

anova_model_SS_50 <- aov(`Soil strength (kPa)` ~ Site + Site: Zone,
                   data=SS_50_data)
anova_model_SS_20 <- aov(`Soil strength (kPa)` ~ Site + Site: Zone,
                         data=SS_20_data)

#3. Assumption test
# Shapiro-Wilk test for normality
residuals_SS_50 <- residuals(anova_model_SS_50)
residuals_SS_20 <- residuals(anova_model_SS_20)
shapiro_test <- shapiro.test(residuals_SS_50)
print(shapiro_test)
shapiro_test <- shapiro.test(residuals_SS_20)
print(shapiro_test)

# Levene’s Test (preferred, more robust to non-normality)
leveneTest(`Soil strength (kPa)` ~ Site * Zone, data = SS_50_data)
leveneTest(`Soil strength (kPa)` ~ Site * Zone, data = SS_20_data)

#SS_20 met the normalty and homogeneity, but SS_50 does not meet any

#Do a log transformation of SS_50 and recheck assumption
SS_50_data$log_SS_50 <- log(SS_50_data$`Soil strength (kPa)`)
anova_model_SS_50_log <- aov(log_SS_50 ~ Site / Zone, data = SS_50_data)
# Recheck assumptions
shapiro.test(residuals(anova_model_SS_50_log))
leveneTest(log_SS_50 ~ Site * Zone, data = SS_50_data)
# QQ plot for visual inspection
qqnorm(residuals(anova_model_SS_50_log))
qqline(residuals(anova_model_SS_50_log), col = "red")

#assumption is basically satisfied after log transformation of SS_50
summary(anova_model_SS_50_log)
summary(anova_model_SS_20)

cld(emmeans(anova_model_SS_50_log, ~ Zone | Site), Letters = letters, adjust = "tukey")
cld(emmeans(anova_model_SS_20, ~ Zone | Site), Letters = letters, adjust = "tukey")

#3 plot the figure
Plot_SS_DI_20 = summarySE(data=Data_DI_20,
                          "Soil strength (kPa)",
                         groupvars="Zone",
                         conf.interval = 0.95,na.rm =T)
Plot_SS_FCM_20 = summarySE(data=Data_FCM_20,
                          "Soil strength (kPa)",
                          groupvars="Zone",
                          conf.interval = 0.95,na.rm =T)
Plot_SS_DI_50 = summarySE(data=Data_DI_50,
                          "Soil strength (kPa)",
                          groupvars="Zone",
                          conf.interval = 0.95,na.rm =T)
Plot_SS_FCM_50 = summarySE(data=Data_FCM_50,
                           "Soil strength (kPa)",
                           groupvars="Zone",
                           conf.interval = 0.95,na.rm =T)

SS_20<-rbind(Plot_SS_DI_20,Plot_SS_FCM_20)
SS_50<-rbind(Plot_SS_DI_50,Plot_SS_FCM_50)

SS_20$Site<-factor(c(rep("Deal Island",4),rep("Farm Creek Marsh",4)),
                   levels=c("Farm Creek Marsh", "Deal Island"))
SS_50$Site<-factor(c(rep("Deal Island",4),rep("Farm Creek Marsh",4)),
                   levels=c("Farm Creek Marsh", "Deal Island"))

SS_DI_20<-c("a","a","a","a") #for Deal Island
SS_DI_50<-c("a","a","a","a") #for Deal Island
SS_FCM_20<-c("ab","b","a","b") #for Deal Island
SS_FCM_50<-c("b","b","a","b") #for Deal Island




##########################################
#5. DI Soil salinity
#5.1 Specify the order of factor levels for plots and Dunnett comparison
library(dplyr)
Data_DI_salinity<-DI_salinity 
Data_DI_salinity =
  mutate(Data_DI_salinity, Month= factor(Month, levels=unique(Month))) #Specify the order of factor levels for plots and Dunnett comparison

Data_FCM_salinity<-FCM_salinity 
Data_FCM_salinity =
  mutate(Data_FCM_salinity, Month = factor(Month, levels=unique(Month))) #Specify the order of factor levels for plots and Dunnett comparison

Plot_Data_DI_salinity = summarySE(data=Data_DI_salinity,
                                  "Salinity_(ppt)",
                                  groupvars="Month",
                                  conf.interval = 0.95,na.rm =T)
Plot_Data_FCM_salinity= summarySE(data=Data_FCM_salinity,
                                  "Salinity_(ppt)",
                                  groupvars="Month",
                                  conf.interval = 0.95,na.rm =T)
Salinity<-rbind(Plot_Data_DI_salinity,Plot_Data_FCM_salinity)
Salinity$Site<-factor(c(rep("Deal Island",6),rep("Farm Creek Marsh",6)),
                      levels=c("Farm Creek Marsh", "Deal Island"))

#Comparision the difference of salinity between DI and FCM

DI_salinity_Vld <-subset(DI_salinity,is.na(`Salinity_(ppt)`)==0)
FCM_salinity_Vld <-subset(FCM_salinity,is.na(`Salinity_(ppt)`)==0)
Salinity_merge<-data.frame("Salinity_(ppt)"=c(DI_salinity_Vld$`Salinity_(ppt)`,
                                      FCM_salinity_Vld$`Salinity_(ppt)`),
                           Month=factor(c(DI_salinity_Vld$Month,
                                   FCM_salinity_Vld$Month)),
                           Site=factor(c(rep("Deal Island",length(DI_salinity_Vld$`Salinity_(ppt)`)),
                                  rep("Farm Creek Marsh", length(FCM_salinity_Vld$`Salinity_(ppt)`)))
                           ))

summarySE(Salinity_merge, "Salinity_.ppt.",
          groupvars=c("Month", "Site"))

lmm_model <- lmer(Salinity_.ppt. ~ Site + (1 | Month), data = Salinity_merge)
summary(lmm_model)

# Test for significance using ANOVA
anova(lmm_model)
# Extract residuals and make QQ plot or Shapiro-Wilk test
qqnorm(resid(lmm_model))
qqline(resid(lmm_model), col = "red")

shapiro.test(resid(lmm_model)) # Only works if n < 5000
plot(fitted(lmm_model), resid(lmm_model))  # Look for random scatter
abline(h = 0, col = "red")

library(DHARMa)
sim_res <- simulateResiduals(fittedModel = lmm_model, n = 1000)
plot(sim_res)          # Visual overall diagnostic
testUniformity(sim_res) # Tests distribution of residuals (normality proxy)
testDispersion(sim_res) # Tests over/under-dispersion (variance issues)
testResiduals(sim_res)  # Combines uniformity + dispersion tests



#Plot salinity
library(ggplot2)

ggplot(Salinity_merge, aes(x = Month, y = Salinity_.ppt., group = Site, color = Site)) +
  geom_line() +
  geom_point() +
  labs(x = "Month", y = "Soil Salinity (ppt)", color = "Site") +
  theme_minimal()
################################################################################
#ELEVATION
################################################################################
#load the elevation data of transplanting waypoints
Waypoint_Elv<-read_excel("Data/RTK elevation survey 2019/RTK data of DI and FCM.xls", 
                         sheet="Transplanting waypoints")
well_substrate<-read_excel("Data/RTK elevation survey 2019/RTK data of DI and FCM.xls", 
                           sheet="Groundwater well")
temp<-summarySE(Waypoint_Elv, "Elevation (m)",
          groupvars=c("Note", "Site", "Treatment"))
Waypoint_Elv_mean<-data.frame(
  Site = temp$Site,
  Note = temp$Note,
  Treatment=temp$Treatment,
  "Elevation (m)"=temp$`Elevation (m)`
)

# Ensure factors are properly set
Waypoint_Elv_mean$Site <- factor(Waypoint_Elv_mean$Site,
                                    levels=unique(c("Farm Creek Marsh","Deal Island")))
Waypoint_Elv_mean$Treatment <- factor(Waypoint_Elv_mean$Treatment,
                                         levels=c("Pond substrate",
                                                  "Pond control",
                                                  "Pond elevated",
                                                  "Dieback",
                                                  "Neighbor patch",
                                                  "Juncus patch"))



# Mixed-effects model
model <- lmer(Elevation..m. ~ Site * Treatment + (1 | Note), data = Waypoint_Elv_mean)

model_FCM<-lmer(Elevation..m. ~ Treatment + (1 | Note), 
               subset(Waypoint_Elv_mean,Site=="Farm Creek Marsh"))
model_DI<-lmer(Elevation..m. ~ Treatment + (1 | Note), 
               subset(Waypoint_Elv_mean,Site=="Deal Island"))
summary(model)

sim_res <- simulateResiduals(fittedModel = model, n = 1000)
plot(sim_res)          # Visual overall diagnostic
testUniformity(sim_res) # Tests distribution of residuals (normality proxy)
testDispersion(sim_res) # Tests over/under-dispersion (variance issues)
testResiduals(sim_res)  # Combines uniformity + dispersion tests
anova(model)
# Tukey-adjusted pairwise comparisons
pairs(emmeans(model, ~ Treatment | Site), adjust = "tukey")

# Get letters for zones at each site
cld(emmeans(model, ~ Treatment | Site), Letters = letters, adjust = "tukey")

Ele_DI_letter<-c("d","c","b","a","a","abc")
Ele_FCM_letter<-c("c","b","a","a","a","a")

letters_df <- data.frame(
  Treatment = rep(c("Pond substrate", "Pond control", "Pond elevated", 
                    "Dieback", "Neighbor patch", "Juncus patch"), 2),
  Site = factor(rep(c("Farm Creek Marsh", "Deal Island"), each = 6), 
                levels=c("Farm Creek Marsh", "Deal Island")),
  Letter = c(Ele_FCM_letter,Ele_DI_letter),
  y=c(0.21,0.32,0.38,0.33,0.31,0.38,
      0.14,0.23,0.38,0.26,0.26,0.24)
)


################################################################################
#PLOT THE FIGURES FOR ALL ENVIRONMENTAL PARAMETERS
################################################################################

#select the color for three zones


col_scheme <- c(
  "Pond substrate" = "#4682B4",  # A strong, appealing blue
  "Pond control" = "#87CEEB",   # Light sky blue
  "Pond" = "#4682B4",           # Same as "Pond control" for consistency
  "Pond elevated" = "#B0E0E6",  # Lighter, pale blue for elevated substrate
  "Dieback" = "#FDAE61",  
  "Dieback patch" = "#FDAE61", # Muted orange for better contrast
  "Neighbor patch" = "#66C2A5", # Light teal
  "Neighbor removal" = "#1B9E77", # Darker teal for emphasis
  "Neighbor control" = "#66C2A5", # Same as "Neighbor patch" for consistency
  "Juncus patch" = "#C2C2B8"   # Neutral light gray
)

#col_scheme <- c(
#  "Pond substrate" = "#4671D5",  # A strong, appealing blue
#  "Pond control" = "#6CAFEA",   # Light sky blue
#  "Pond" = "#6CAFEA",           # Same as "Pond control" for consistency
#  "Pond elevated" = "#A7D8F5",  # Lighter, pale blue for elevated substrate
#  "Dieback" = "#F2A65A",        # Muted orange for better contrast
#  "Dieback patch" = "#F2A65A",
#  "Neighbor patch" = "#88C586", # Light teal
#  "Neighbor removal" = "#4E9A42", # Darker teal for emphasis
#  "Neighbor control" = "#88C586", # Same as "Neighbor patch" for consistency
#  "Juncus patch" = "#C2C2B8"    # Neutral light gray
#)


custom_colors3<-col_scheme[c("Pond","Dieback","Neighbor patch")] 
custom_colors4<-col_scheme[c("Pond","Dieback","Neighbor patch","Juncus patch")] 

#Now we need to call the dataframes for ploting:
#hydrpd,Salinity,Redox,SS_20,SS_50

#elevation
plot1<-ggplot(Waypoint_Elv_mean, aes(x = Treatment, y = Elevation..m., fill = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  geom_text(
    data = letters_df,
    aes(x = Treatment, 
        y = y + 0.03, 
        label = Letter),  # Adjust y as needed for position
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = col_scheme[levels(Waypoint_Elv_mean$Treatment)]
  ) +
  theme_bw()  +
  labs(x="Marsh zone/transplant",y = "Elevation (m)"
  ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles


#hydrpd
hydrpd<-subset(hydrpd, Zone!="Juncus patch")
plot2<-ggplot(hydrpd,aes(x= Zone,y=`Hydroperiod`*100,fill=Zone))+
  geom_bar(stat="identity", position = "dodge", width = 0.7,
           show.legend = TRUE) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(hydrpd$Zone)]) +
  labs(x="Marsh zone",y = "Hydroperiod (%)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles


#redox
plot3<-ggplot(Redox,
       aes(x = Zone, y = `Redox (mV)`, fill = Zone)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(Redox$Zone)]) +
  ylim(-500,0)+
  # Add error bars with different colors for each site
  geom_errorbar(aes(ymax = `Redox (mV)` + se, ymin = `Redox (mV)` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  # Add text labels with different colors for each site
  geom_text(aes(label = c("", "b", "a", "b", "", "b", "a", "ab"), 
                y = `Redox (mV)` - se - 40), 
            position = position_dodge(width = 0.7), 
            vjust = 0, size = 3, show.legend = FALSE) +
  labs(x="Marsh zone", y = "Soil redox potential (mV)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles


plot4<-ggplot(SS_20,
       aes(x = Zone, y = `Soil strength (kPa)`, fill = Zone)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(SS_20$Zone)]) +
  ylim(0, 100)+
  # Add error bars with different colors for each site
  geom_errorbar(aes(ymax = `Soil strength (kPa)` + se, ymin = `Soil strength (kPa)` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  
  # Add text labels with different colors for each site
  geom_text(aes(label = c(SS_DI_20, SS_FCM_20), 
                y = `Soil strength (kPa)` + se +2), 
            position = position_dodge(width = 0.7), 
            vjust = 0, size = 3, show.legend = FALSE) +
  
  labs(x = "Marsh zone", y = "Soil strength_20 cm (kPa)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles

plot5<-ggplot(SS_50,
       aes(x = Zone, y = `Soil strength (kPa)`, fill = Zone)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(SS_50$Zone)]) +
  ylim(0, 100)+
  # Add error bars with different colors for each site
  geom_errorbar(aes(ymax = `Soil strength (kPa)` + se, ymin = `Soil strength (kPa)` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  
  # Add text labels with different colors for each site
  geom_text(aes(label = c(SS_DI_50, SS_FCM_50), 
                y = `Soil strength (kPa)` + se +2),
            position = position_dodge(width = 0.7), 
            vjust = 0, size = 3, show.legend = FALSE) +
  
  labs(x = "Marsh zone", y = "Soil strength_50 cm (kPa)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles




tiff("Result/Soil salinity.tiff", unit="in", width=3, height=2, res=600, pointsize=10)
#salinity
ggplot(Salinity, aes(x = Month, y = `Salinity_(ppt)`,colour =  Site, group = Site)) +
  geom_line(size = 0.5) +  # Draw lines for each site
  scale_color_manual(values = c("Farm Creek Marsh" = "lightseagreen", "Deal Island" = "lightseagreen")) +  # Use custom colors for lines
  #geom_point(aes(shape = Site), color = "black",size = 2, ) +  # Add points at each data value
  geom_errorbar(aes(ymin = `Salinity_(ppt)` - se, ymax = `Salinity_(ppt)` + se),
                width = 0.2, size = 0.2) +  # Add error bars
  ylim(0,30)+
  scale_shape_manual(values = c("Farm Creek Marsh" = 1, "Deal Island" = 16)) +  # Open and closed circles
  labs(x = "Month-Year", y = "Soil Salinity (ppt)", color = "Marsh site", shape = "Marsh site") +  # Axis and legend labels
  theme_bw() +  # Use a clean theme
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),  # Rotate x-axis labels
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        plot.title = element_text(size = rel(1.5), face = "bold", vjust = 1.5),
        axis.title = element_text(face = "plain"),
        axis.title.y = element_text(vjust = 1.8),
        axis.title.x = element_text(vjust = -0.5),
        panel.border = element_rect(color = "black"),
        legend.position = "none")  # Hide the legend

dev.off()




# Modify to remove x-axis labels for the upper panel
tiff("Result/Hydrology1.tiff", unit="in", width=5, height=5, res=600, pointsize=10)
final_plot <- (plot1 | plot2)
final_plot
dev.off()

tiff("Result/Environmental variables.tiff", unit="in", width=5, height=5, res=600, pointsize=10)
final_plot <- (plot4 | plot5 | plot3) 
final_plot
dev.off()

