#statistic test dieback
rm(list=ls())

library(here)
library(lme4)
library(lmerTest)
library(dplyr)
library(readxl)

setwd(here::here())
# Read data
abio <- read_excel("Data/Vegetation survey/Plant diversity source data.xlsx")
phy_data <- read_excel("Data/Vegetation survey/Pore water Salinity, pH, Redox, water depth.xlsx", sheet = "clean data")


# ===== Calculate dieback  =====
# calculate dieback area
abio <- abio %>%
  mutate(Dieback_area = `Open Water` + `Bare ground`)

# fit the linear model and conduct ANOVA
model = lm(Dieback_area ~ Point + Site + Point:Site,
           data=abio)

library(car)

Anova(model, type="II")  
#check assumption
shapiro.test(residuals(model)) #dieback area failed the normality test
leveneTest(Dieback_area ~ Point * Site, 
           data = abio) #dieback area passed the heterogeneity test
#log transform
abio$Dieback_area_log <- 
  log(abio$Dieback_area+1)
#refit the model
# fit the linear model and conduct ANOVA
model_log = lm(Dieback_area_log ~ Point + Site + Point:Site,
           data=abio)
#recheck assumption
shapiro.test(residuals(model_log)) #dieback area passed the normality test
leveneTest(Dieback_area_log ~ Point * Site, 
           data = abio) #dieback area passed the heterogeneity test
Anova(model_log, type="II")  
# ===== Water depth  =====
model = lm(`Water depth/cm`~ Point + Site + Point:Site,
           phy_data)
#check assumption
shapiro.test(residuals(model)) #water depth failed the normality test
leveneTest(Dieback_area ~ Point * Site, 
           data = abio) #water depth passed the heterogeneity test
#log transform
phy_data$water_depth_log <- 
  log(phy_data$`Water depth/cm`+abs(min(phy_data$`Water depth/cm`))+1)

#refit the model
# fit the linear model and conduct ANOVA
model_log = lm(water_depth_log ~ Point + Site + Point:Site,
               data=phy_data)
#recheck assumption
shapiro.test(residuals(model_log)) #water depth_log failed the normality test
leveneTest(Dieback_area_log ~ Point * Site, 
           data = abio) #water depth_log passed the heterogeneity test
#ART transformation
library(ARTool)
# Run ART model
phy_data$Site <- factor(phy_data$Site)
phy_data$Point <- factor(phy_data$Point)
art_model <- art(`Water depth/cm` ~ Site * Point, data = phy_data)
anova(art_model)

# Post-hoc for main effects
art_site <- artlm(art_model, "Site")
art_point <- artlm(art_model, "Point")
art_interaction <- artlm(art_model, "Site:Point")

# Get letters for zones at each site
cld(emmeans(art_interaction, ~ Point | Site), Letters = letters, adjust = "tukey")
