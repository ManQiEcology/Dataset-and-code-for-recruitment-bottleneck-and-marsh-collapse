rm(list=ls())
library(readxl)
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

col_scheme <- c(
  "Pond substrate" = "#4682B4",  # A strong, appealing blue
  "Pond control" = "#87CEEB",   # Light sky blue
  "Pond" = "#4682B4",           # Same as "Pond control" for consistency
  "Pond elevated" = "#B0E0E6",  # Lighter, pale blue for elevated substrate
  "Dieback" = "#FDAE61",  
  "Dieback patch" = "#FDAE61", # Muted orange for better contrast
  "Neighbor patch" = "#66C2A5", # Light teal
  "Neighbor removal" = "#1B9E77", # Darker teal for emphasis
  "Neighbor intact" = "#66C2A5", # Same as "Neighbor patch" for consistency
  "Juncus patch" = "#C2C2B8"   # Neutral light gray
)


setwd(here::here())

DI_plant<- read_excel("Data/Pond plant transplanting result-2019/DI.xlsx", 
                      sheet = "Plant growth",na = c("NA", "missing", ""))
FCM_plant<- read_excel("Data/Pond plant transplanting result-2019/FCM.xlsx", 
                       sheet = "Plant growth",na = c("NA", "missing", ""))
Biomass_DI<-read_excel("Data/Pond plant transplanting result-2019/DI.xlsx", 
                       sheet = "Biomass",na = c("NA", "missing", ""))
Biomass_FCM<- read_excel("Data/Pond plant transplanting result-2019/FCM.xlsx", 
                       sheet = "Biomass",na = c("NA", "missing", ""))
Flower_DI<-read_excel("Data/Pond plant transplanting result-2019/DI.xlsx", 
                    sheet = "Flowering",na = c("NA", "missing", ""))
Flower_FCM<-read_excel("Data/Pond plant transplanting result-2019/FCM.xlsx", 
                       sheet = "Flowering",na = c("NA", "missing", ""))
Survival_DI<-read_excel("Data/Pond plant transplanting result-2019/DI.xlsx", 
                        sheet = "Survival rate",na = c("NA", "missing", ""))
Survival_FCM<-read_excel("Data/Pond plant transplanting result-2019/FCM.xlsx", 
                        sheet = "Survival rate",na = c("NA", "missing", ""))

  
Shoot_Data<-rbind(FCM_plant[,1:9],DI_plant[,1:9])
Biomass_Data<-rbind(Biomass_FCM,Biomass_DI)
Flower_Data<-rbind(Flower_FCM,Flower_DI)
Survival_Data<-rbind(Survival_DI,Survival_FCM)


Shoot_Data$Site<-factor(c(rep("Farm Creek Marsh", nrow(FCM_plant)),
                          rep("Deal Island", nrow(DI_plant))),
                        levels = c("Farm Creek Marsh","Deal Island"))
Biomass_Data$Site<-factor(c(rep("Farm Creek Marsh", nrow(Biomass_FCM)),
                            rep("Deal Island", nrow(Biomass_DI))),
                          levels = c("Farm Creek Marsh","Deal Island"))
Flower_Data$Site<-factor(c(rep("Farm Creek Marsh", nrow(Flower_FCM)),
                            rep("Deal Island", nrow(Flower_DI))),
                          levels = c("Farm Creek Marsh","Deal Island"))
Shoot_Data$`Average height (cm)`<-apply(as.matrix(Shoot_Data[4:8]),1,FUN=mean,na.rm=T)

Shoot_Data$Month[which(Shoot_Data$Month=="Jul-21")]<-"Jun-21"
Flower_Data$Month[which(Flower_Data$Month=="Jul-21")]<-"Jun-21"
Shoot_Data$Month<-factor(Shoot_Data$Month,levels = unique(Shoot_Data$Month))
Flower_Data$Month<-factor(Flower_Data$Month,levels = unique(Flower_Data$Month))
Survival_Data$Month<-factor(Survival_Data$Month,levels = unique(Shoot_Data$Month))
Shoot_Data_21<-subset(Shoot_Data, Month=="Jun-21")
Flower_Data_21<-subset(Flower_Data, Month=="Jun-21")
Flower_Data_20<-subset(Flower_Data, Month=="Jun-20")

#_______________________-shoot height
Shoot_height_april<-subset(Shoot_Data, Month=="Apr-19" & Site=="Farm Creek Marsh")
anova_model_height<- aov(`Average height (cm)`~ Treatment,
                             data=Shoot_height_april)
summary(anova_model_height)
# Check assumptions
shapiro.test(residuals(anova_model_height)) #height_GR didn't pass
leveneTest(`Average height (cm)` ~ Site * Treatment, data = Shoot_height_april) #height_GR didn't pass
# Post-hoc for main effects
zone<-PostHocTest(anova_model_height, "Treatment")


#__________________shoot number
Shoot_no_april<-subset(Shoot_Data, Month=="Jun-20" & Site=="Farm Creek Marsh")
anova_model_shoot_no<- aov(`Shoot No.`~ Treatment,
                         data=Shoot_height_april)
summary(anova_model_height)
# Check assumptions
shapiro.test(residuals(anova_model_height)) #height_GR didn't pass
leveneTest(`Average height (cm)` ~ Site * Treatment, data = Shoot_height_april) #height_GR didn't pass
# Post-hoc for main effects
zone<-PostHocTest(anova_model_shoot_no, "Treatment")







#calculate growth rate of shoot height and number
Apr19<-subset(Shoot_Data, `Month`=="Apr-19") # get the plant data in April (inital month)
Jun21<-subset(Shoot_Data,`Month`=="Jun-20") # get the plant data in June 2020 (the first monitoring of 2020)
Shoot_GR<-data.frame("ID"=Apr19[["ID"]],
                     "Treatment"=Apr19[["Treatment"]],
                     "Site"=Apr19[["Site"]],
                     "shoot_no_GR"=(Jun21$`Shoot No.`- Apr19$`Shoot No.`)/Apr19$`Shoot No.`,
                     "Height_GR"=(Jun21$`Height (cm)`-Apr19$`Height (cm)`)/Apr19$`Height (cm)`)

factorize_treatment<-function(x) {
  x<-as.data.frame(x)
  x$Treatment<-factor(x$Treatment,
                      levels = c("Pond control",
                                 "Pond elevated",
                                 "Dieback",
                                 "Neighbor intact",
                                 "Neighbor removal",
                                 "Juncus patch"))
  return(x)
}

Shoot_Data<-factorize_treatment(Shoot_Data)
Shoot_GR<-factorize_treatment(Shoot_GR)
Biomass_Data<-factorize_treatment(Biomass_Data)
Flower_Data_21<-factorize_treatment(Flower_Data_21)
Survival_Data<-factorize_treatment(Survival_Data)

sum_height = summarySE(Shoot_Data,
                measurevar="Average height (cm)",
                groupvars=c("Site", "Treatment","Month"),na.rm = T)

sum_shoot_no =summarySE(Shoot_Data,
                        measurevar="Shoot No.",
                        groupvars=c("Site", "Treatment","Month"),na.rm = T)

sum_height_2021<-summarySE(Shoot_Data_21,
                           measurevar="Average height (cm)",
                           groupvars=c("Site", "Treatment"),na.rm = T)
sum_shoot_no_2021<-summarySE(Shoot_Data_21,
                             measurevar="Shoot No.",
                             groupvars=c("Site", "Treatment"),na.rm = T)

sum_biomass_alive = summarySE(Biomass_Data,
                            measurevar="Live biomass (g)",
                            groupvars=c("Site", "Treatment"),na.rm = T)

sum_biomass_dead =  summarySE(Biomass_Data,
                              measurevar="Dead biomass (g)",
                              groupvars=c("Site", "Treatment"),na.rm = T)

sum_flower_21 = summarySE(Flower_Data_21,
                       measurevar="Flowering",
                       groupvars=c("Site", "Treatment","Month"),na.rm = T)

sum_heigh_GR = summarySE(Shoot_GR,
                       measurevar="Height_GR",
                       groupvars=c("Site", "Treatment"),na.rm = T)
sum_shoot_no_GR = summarySE(Shoot_GR,
                            measurevar="shoot_no_GR",
                            groupvars=c("Site", "Treatment"),na.rm = T)


sum_height<-factorize_treatment(sum_height)
sum_shoot_no<-factorize_treatment (sum_shoot_no)
sum_height_2021<-factorize_treatment (sum_height_2021)
sum_shoot_no_2021<-factorize_treatment(sum_shoot_no_2021)
sum_biomass_alive<-factorize_treatment(sum_biomass_alive)
sum_biomass_dead<-factorize_treatment(sum_biomass_dead)
sum_flower_21<-factorize_treatment(sum_flower_21)
sum_heigh_GR <- factorize_treatment(sum_heigh_GR)
sum_shoot_no_GR <- factorize_treatment(sum_shoot_no_GR)
################################################################################
#STATISTICAL ANALYSIS
################################################################################
#Nested ANOVA
library(emmeans)
library(multcompView)
stts_display<-function(x) {
  print(summary(x))
  pairwise<-TukeyHSD(x)$`Site:Treatment`
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

anova_model_height_GR <- aov(`Height_GR` ~ Site + Site: Treatment,
                   data=Shoot_GR)
anova_model_shoot_no_GR <- aov(`shoot_no_GR` ~ Site + Site: Treatment,
                             data=Shoot_GR)
anova_model_flower <- aov(Flowering ~ Site + Site: Treatment,
                               data=Flower_Data_21)

# Check assumptions
shapiro.test(residuals(anova_model_height_GR)) #height_GR didn't pass
leveneTest(Height_GR ~ Site * Treatment, data = Shoot_GR) #height_GR didn't pass

shapiro.test(residuals(anova_model_shoot_no_GR)) 
leveneTest(`shoot_no_GR` ~ Site * Treatment, data = Shoot_GR) #shoot_no passed the assumption test
summary(anova_model_shoot_no_GR)
cld(emmeans(anova_model_shoot_no_GR,~ Treatment | Site), Letters = letters, 
    decreasing = TRUE, adjust = "sidak")

shapiro.test(residuals(anova_model_flower)) #flowering didn't pass 
leveneTest(Flowering ~ Site * Treatment, data = Flower_Data_21) #flowering didn't pass 

#log transformation for shoot_GR and flowering
Shoot_GR$log_hgr<-log(Shoot_GR$Height_GR + abs(min(Shoot_GR$Height_GR,na.rm = T)) + 1)
anova_model_height_GR_log <- aov(log_hgr ~ Site + Site: Treatment, data = Shoot_GR)

Flower_Data_21$Flowering_log<-log(Flower_Data_21$Flowering + 1)
anova_model_flower_log <- aov(Flowering_log ~ Site + Site: Treatment, data = Flower_Data_21)

# Recheck assumptions
shapiro.test(residuals(anova_model_height_GR_log)) #log transformed height_GR didn't pass
leveneTest(log_hgr ~ Site * Treatment, data = Shoot_GR) #log transformed height_GR didn't pass

shapiro.test(residuals(anova_model_flower_log)) #flowering_log didn't pass 
leveneTest(Flowering_log ~ Site * Treatment, data = Flower_Data_21) #flowering_log didn't pass 

# Choosing ART for height_GR and Flowering
Shoot_GR_clean <- na.omit(Shoot_GR[, c("Height_GR", "Site", "Treatment")])
art_model_Height_GR <- art(Height_GR ~ Site * Treatment, data = Shoot_GR_clean)
anova(art_model_Height_GR)

Flower_Data_21_clean <- na.omit(Flower_Data_21[, c("Flowering", "Site", "Treatment")])
art_model_flower <- art(Flowering ~ Site * Treatment, data = Flower_Data_21_clean)
anova(art_model_flower)


# Post-hoc for main effects
art_site <- artlm(art_model_Height_GR, "Site")
art_zone <- artlm(art_model_Height_GR, "Treatment")
art_interaction <- artlm(art_model_Height_GR, "Site:Treatment")

art_site <- artlm(art_model_flower, "Site")
art_zone <- artlm(art_model_flower, "Treatment")
art_interaction <- artlm(art_model_flower, "Site:Treatment")

# Tukey-adjusted pairwise comparisons
pairs(emmeans(art_interaction, ~ Treatment | Site), adjust = "sidak")
# Get letters for zones at each site
cld(emmeans(art_interaction, ~ Treatment | Site), Letters = letters, decreasing = T, adjust = "sidak")
cld(emmeans(art_zone, ~ Treatment ), Letters = letters, adjust = "sidak")
cld(emmeans(art_site, ~ Site ), Letters = letters, adjust = "sidak")



Letters_plant<-data.frame("Height_GR_FCM"=c("a",	"a",	"a",	"a",	"a",	"a"),
                          "Height_GR_DI"= c("a",	"a",	"a",	"a",	"a",	"a"),
                          "shoot_no_FCM"=c("d",	"abc",	"ab",	"c",	"a",	"bc"),
                          "shoot_no_DI"=c("c",	"a",	"ab",	"b",	"ab",	"bc"),
                          "Flower_FCM"=c("a",	"b",	"a",	"ab",	"a",	"ab"),
                          "Flower_DI"=c("bc",	"a",	"c",	"ab",	"bc",	"abc"))


################################################################################
#PLOTTING THE RESULT
################################################################################
#Shoot height and number
sum_height_june<-subset(sum_height,Month=="Apr-19"|Month=="Jun-19"|Month=="Jun-20"|Month=="Jun-21")
sum_shoot_no_june<-subset(sum_shoot_no,Month=="Apr-19"|Month=="Jun-19"|Month=="Jun-20"|Month=="Jun-21")
plot_height_series<-ggplot(sum_height_june,
       aes(x = Month, y = `Average height (cm)`, group = Treatment, color = Treatment)) +
  ylim(0, 120) +
  geom_line(size = 1) +  # Add lines for each Treatment
  #geom_point(size = 3) +  # Add points at each data value
  scale_color_manual(values = col_scheme[levels(sum_height$Treatment)]) +
  geom_errorbar(aes(ymax = `Average height (cm)` + se, 
                    ymin = `Average height (cm)` - se),
                width = 0.2, size = 0.5) +  # Error bars
  
  labs(
    x = "Time (MM-YY)",
    y = "Average Height (cm)",
    color = "Treatment"
  ) +
  facet_wrap(~ Site,scales = "free_x", ncol = 1) +
  theme_bw()+ # Optional: cleaner theme
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

plot_shoot_no_series<-ggplot(sum_shoot_no_june,
              aes(x = Month, y = `Shoot No.`, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +  # Add lines for each Treatment
  #geom_point(size = 3) +  # Add points at each data value
  scale_color_manual(values = col_scheme[levels(sum_shoot_no$Treatment)]) +
  geom_errorbar(aes(ymax = `Shoot No.` + se, 
                    ymin = `Shoot No.` - se),
                width = 0.2, size = 0.5) +  # Error bars
  
  labs(
    x = "Time (MM-YY)",
    y = "Shoot number",
    color = "Treatment"
  ) +
  facet_wrap(~ Site,scales = "free_x", ncol = 1) +
  theme_bw()+ # Optional: cleaner theme
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        #legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles

plot_survival_series<-ggplot(Survival_Data,
                             aes(x = Month, y = `Survival rate (%)`, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +  # Add lines for each Treatment
  scale_color_manual(values = col_scheme[levels(Survival_Data$Treatment)]) +
  
  labs(
    x = "Time (MM-YY)",
    y = "Survival rate (%)",
    color = "Treatment"
  ) +
  ylim(0,100)+
  facet_wrap(~ Site,scales = "free_x", nrow = 1) +
  theme_bw()+ # Optional: cleaner theme
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        #legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles

plot_survival<-ggplot(subset(Survival_Data,Month=="Jun-21"),
                             aes(x = Treatment, y = `Survival rate (%)`,  group_by= Treatment, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(Survival_Data$Treatment)]) +
  labs(
    x = "Treatment",
    y = "Survival rate in June 2021 (%) ",
    color = "Treatment"
  ) +
  ylim(0,100)+
  theme_bw()+ # Optional: cleaner theme
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        #legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles

plot_height_GR<-ggplot(sum_heigh_GR,
               aes(x =Treatment, y = `Height_GR`, group = Treatment, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(sum_heigh_GR$Treatment)]) +
  ylim(-1,0)+
  geom_errorbar(aes(ymax = `Height_GR` + se, ymin = `Height_GR` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  # Add text labels with different colors for each site
  geom_text(aes(label = c(Letters_plant$Height_GR_FCM,Letters_plant$Height_GR_DI), 
                y = `Height_GR` + se*`Height_GR`/abs(`Height_GR`) +
                  `Height_GR`/abs(`Height_GR`)*0.1), 
            position = position_dodge(width = 0.7), 
            vjust = 0, size = 2, show.legend = FALSE) +
  labs(
    x = "Treatment",
    y = "Relative increase of shoot height",
    color = "Treatment"
  ) +
  theme_bw()+ # Optional: cleaner theme
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

plot_shoot_no_GR<-ggplot(sum_shoot_no_GR ,
               aes(x =Treatment, y = `shoot_no_GR`, group = Treatment, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(sum_shoot_no_GR$Treatment)]) +
  ylim(-1.5,2)+
  geom_errorbar(aes(ymax = `shoot_no_GR` + se, ymin = `shoot_no_GR` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  geom_text(aes(label = c(Letters_plant$shoot_no_FCM,Letters_plant$shoot_no_DI), 
                y = ifelse(shoot_no_GR>=0, `shoot_no_GR` + se + 0.1,
                           `shoot_no_GR` - se - 0.3)), 
            position = position_dodge(width = 0.7), 
            vjust = 0, size = 2, show.legend = FALSE) +
  labs(
    x = "Treatment",
    y = "Relative increase of shoot number",
    color = "Treatment"
  ) +
  theme_bw()+ # Optional: cleaner theme
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



plot_flower_21<-ggplot(subset(sum_flower_21, Month=="Jun-21"),
               aes(x =Treatment, y = Flowering, group = Treatment, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(sum_flower_21$Treatment)]) +
  geom_errorbar(aes(ymax = Flowering + se, ymin = Flowering - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  ylim(0,25)+
  geom_text(aes(label = c(Letters_plant$Flower_FCM,Letters_plant$Flower_DI), 
                y = Flowering + se*Flowering/abs(Flowering) + 
                  Flowering/abs(Flowering)*0.5), 
            position = position_dodge(width = 0.7), 
            vjust = 0, size = 2, show.legend = FALSE) +
  labs(
    x = "Treatment",
    y = "Inflorescence number in June 2021",
    color = "Treatment"
  ) +
  theme_bw()+ # Optional: cleaner theme
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


tiff("Result/Shoot height and number.tiff", unit="in", width=7, height=5, res=600, pointsize=10)
final_plot <- (plot_height_series | plot_shoot_no_series)
final_plot
dev.off()

tiff("Result/survival rate_series.tiff", unit="in", width=5, height=3, res=600, pointsize=10)
plot_survival_series
dev.off()



tiff("Result/plant performance in 2021.tiff", unit="in", width=7, height=5, res=600, pointsize=10)
final_plot <- (plot_height_GR |plot_shoot_no_GR|  plot_flower_21|plot_survival)
final_plot
dev.off()





#EXTRA FIGURES
plot_height_21<-ggplot(sum_height_2021,
                       aes(x =Treatment, y = `Average height (cm)`, group = Treatment, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(ssum_height_2021$Treatment)]) +
  geom_errorbar(aes(ymax = `Average height (cm)` + se, ymin = `Average height (cm)` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  
  labs(
    x = "Treatment",
    y = "Shoot height (cm)",
    color = "Treatment"
  ) +
  theme_bw()+ # Optional: cleaner theme
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

plot_shoot_no_21<-ggplot(sum_shoot_no_2021,
                         aes(x =Treatment, y = `Shoot No.`, group = Treatment, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(sum_shoot_no_2021$Treatment)]) +
  geom_errorbar(aes(ymax = `Shoot No.` + se, ymin = `Shoot No.` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  
  labs(
    x = "Treatment",
    y = "Shoot number",
    color = "Treatment"
  ) +
  theme_bw()+ # Optional: cleaner theme
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

plot_biomass_alive<-ggplot(sum_biomass_alive,
                           aes(x =Treatment, y = `Live biomass (g)`, group = Treatment, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(sum_shoot_no$Treatment)]) +
  geom_errorbar(aes(ymax = `Live biomass (g)` + se, ymin = `Live biomass (g)` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  
  labs(
    x = "Treatment",
    y = "Live biomass (g)",
    color = "Treatment"
  ) +
  theme_bw()+ # Optional: cleaner theme
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

plot_biomass_dead<-ggplot(sum_biomass_dead,
                          aes(x =Treatment, y = `Dead biomass (g)`, group = Treatment, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Site, scales = "free_x", ncol = 1) +  # Arrange facets vertically
  scale_fill_manual(values = col_scheme[levels(sum_shoot_no$Treatment)]) +
  geom_errorbar(aes(ymax = `Dead biomass (g)` + se, ymin = `Dead biomass (g)` - se), 
                position = position_dodge(width = 0.7),
                width = 0.2, size = 0.3, show.legend = FALSE) +
  
  labs(
    x = "Treatment",
    y = "Dead biomass (g)",
    color = "Treatment"
  ) +
  theme_bw()+ # Optional: cleaner theme
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        #legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles

ggplot(sum_flower,
       aes(x = Month, y = Flowering, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +  # Add lines for each Treatment
  #geom_point(size = 3) +  # Add points at each data value
  scale_color_manual(values = col_scheme[levels(sum_flower$Treatment)]) +
  geom_errorbar(aes(ymax = Flowering + se, 
                    ymin = Flowering - se),
                width = 0.2, size = 0.5) +  # Error bars
  
  labs(
    x = "Month",
    y = "Average Height (cm)",
    color = "Treatment"
  ) +
  facet_wrap(~ Site,scales = "free_x", ncol = 1) +
  theme_bw()+ # Optional: cleaner theme
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(size = rel(1.5),
                                  face = "plain", vjust = 1.5), 
        axis.title = element_text(face = "plain"),
        #legend.position = "none",  # Hide the legend
        axis.title.y = element_text(vjust= 1.8),
        axis.text.x = element_text(angle = 90,vjust=0.5, hjust = 1),
        #axis.ticks.x = element_blank(), # Remove x-axis text
        strip.text = element_blank())   # Remove facet titles




