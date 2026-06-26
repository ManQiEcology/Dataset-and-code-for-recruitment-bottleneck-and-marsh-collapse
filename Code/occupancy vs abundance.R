rm(list=ls())

library(vegan)
library(readxl)
library(tidyverse)
library(patchwork)
library(here)
library(ggtext)
library(dplyr)

setwd((here::here()))

# Read data
abio <- read_excel("Data/Vegetation survey/Plant diversity source data.xlsx")
phy_data <- read_excel("Data/Vegetation survey/Pore water Salinity, pH, Redox, water depth.xlsx", sheet = "clean data")

# 90% Confidence Interval Function
ci90 <- function(x) {
  n <- sum(!is.na(x))
  se <- sd(x, na.rm = TRUE) / sqrt(n)
  t_val <- qt(0.95, df = n - 1)
  mean_x <- mean(x, na.rm = TRUE)
  c(mean = mean_x, lower = mean_x - t_val * se, upper = mean_x + t_val * se)
}

# ===== Environmental data (top panel) =====
# Water depth
water_ci <- phy_data %>%
  group_by(Point) %>%
  summarise(
    mean = ci90(`Water depth/cm`)[1],
    lower = ci90(`Water depth/cm`)[2],
    upper = ci90(`Water depth/cm`)[3]
  ) %>%
  mutate(Variable = "Water depth")

# Bare + Open Cover
abio <- abio %>%
  mutate(Dieback_area = `Open Water` + `Bare ground`)

bare_ci <- abio %>%
  group_by(Point) %>%
  summarise(
    mean = ci90(Dieback_area)[1],
    lower = ci90(Dieback_area)[2],
    upper = ci90(Dieback_area)[3]
  ) %>%
  mutate(Variable = "Dieback area")

# Rescale water depth to match cover %
scale_factor <- max(bare_ci$upper) / max(water_ci$upper)
water_ci_scaled <- water_ci %>%
  mutate(across(c(mean, lower, upper), ~ .x * scale_factor))

# Combine
env_ci_dual <- bind_rows(
  bare_ci,
  water_ci_scaled %>% mutate(Variable = "Water depth")
)

env_ci_dual$Point <- factor(env_ci_dual$Point, levels = c("A", "B", "C", "D", "E"), ordered = TRUE)

# Environmental panel
env_plot <- ggplot(env_ci_dual, aes(x = Point, y = mean, color = Variable, fill = Variable, group = Variable)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.1) +
  scale_color_manual(values = c("Dieback area" = "#D95F02", "Water depth" = "#1B9E77")) +
  scale_fill_manual(values = c("Dieback area" = "#D95F02", "Water depth" = "#1B9E77")) +
  xlab("Sampling point")+
  scale_y_continuous(
    name = "Dieback area (%)",
    sec.axis = sec_axis(~ . / scale_factor, 
                        name = "Water depth (cm)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5),  # reduces white space around plot
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = ggtext::element_markdown(size = 10),
    axis.title.y.left = element_text(color = "#D95F02", size = 10, angle = 90, vjust = -22),
    axis.title.y.right = element_text(color = "#1B9E77", size = 10, angle = 90),
    legend.position = "none"
  )


# ===== Species heatmap (bottom panel) =====
# Long-format abundance data
flood_intolerant_species <- c("Spartina patens",
                              "Distichlis spicata",
                              "Iva frutescens", 
                              "Phragmites australis")
# Reorder species by ecological importance or flood tolerance
flood_tolerance_order <- c(
  "Juncus roemerianus",
  "Schoenoplectus americanus",
  "Spartina alterniflora",
  "Bolboschoenus fluviatilis",
  "Salicornia sp.",
  "Pluchea odorata",
  "Setaria geniculate",
  "Distichlis spicata",
  "Spartina patens",
  "Iva frutescens",
  "Phragmites australis"
)
Site_order <- c(
  "S1","S2", "S4", "S3"
)



long_data <- abio %>%
  pivot_longer(cols = `Phragmites australis`:`Pluchea odorata`,
               names_to = "Species", values_to = "Abundance")

# Filter: keep species that appear in ≥2 quadrats and reach ≥10% abundance
species_filter <- long_data %>%
  group_by(Species) %>%
  summarise(
    occupancy = sum(Abundance > 0),
    max_abundance = max(Abundance, na.rm = TRUE)
  ) %>%
  filter(occupancy >= 2, max_abundance >= 10) %>%
  pull(Species)

# Filter data
filtered_long_data <- long_data %>%
  filter(Species %in% species_filter)

filtered_long_data <- filtered_long_data %>%
  mutate(FloodClass = if_else(
    Species %in% flood_intolerant_species,
    "Flood-Sen",
    "Flood-Tol"
  ))

total_quadrats_per_point_site <- 3 # 


occupancy_abundance_by_point_site_species <- filtered_long_data %>%
  group_by(Point, Site, Species) %>%
  summarise(
    Occupancy = 100 * sum(Abundance > 0) / total_quadrats_per_point_site,
    MeanAbundance = mean(Abundance[Abundance > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Occupancy > 0)

occupancy_abundance_by_point_site <- filtered_long_data %>%
  group_by(Point, Site, Transect, FloodClass) %>%
  summarise(
    Presence = ifelse(sum(Abundance > 0), 1, 0),
    Abundance = sum(Abundance[Abundance > 0], na.rm = TRUE),
    .groups = "drop"
  ) 
occupancy_abundance_by_point_site <- occupancy_abundance_by_point_site %>%
  group_by(Point, Site,FloodClass) %>%
  summarise(
    Occupy=sum(Presence>0),
    Occupancy = 100 * sum(Presence > 0) / total_quadrats_per_point_site,
    MeanAbundance = mean(Abundance[Abundance > 0], na.rm = TRUE),
    .groups = "drop"
  ) 

occupancy_abundance_by_point_site_species <- filtered_long_data %>%
  group_by(Point, Site, Species) %>%
  summarise(
    Occupancy = 100 * sum(Abundance > 0) / total_quadrats_per_point_site,
    MeanAbundance = mean(Abundance[Abundance > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Occupancy > 0)

occupancy_abundance_by_point_site_species$Species <- factor(
  occupancy_abundance_by_point_site_species$Species,
  levels = flood_tolerance_order
)

occupancy_abundance_by_point_site_species$Site <- factor(
  occupancy_abundance_by_point_site_species$Site,
  levels = Site_order
)
# Define colors to match panel b
flood_colors <- c("Flood-Sen" = "#66c2a5", "Flood-Tol" = "#fc8d62")

occupancy_abundance_by_point_site_species <- occupancy_abundance_by_point_site_species %>%
  mutate(FloodClass = if_else(
    Species %in% flood_intolerant_species,
    "Flood-Sen",
    "Flood-Tol"
  ))

# Create a new column with HTML-style colored species labels
occupancy_abundance_by_point_site_species <- occupancy_abundance_by_point_site_species %>%
  mutate(Species_colored = paste0(
    "<span style='color:", flood_colors[FloodClass], "'><i>", Species, "</i></span>"
  ))


# Get order of Species and apply it to Species_colored
species_order <- occupancy_abundance_by_point_site_species %>%
  distinct(Species, Species_colored) %>%
  arrange(match(Species, flood_tolerance_order))  # Preserve original order

occupancy_abundance_by_point_site_species$Species_colored <- factor(
  occupancy_abundance_by_point_site_species$Species_colored,
  levels = species_order$Species_colored
)



# ===== difference in occupancy between flood tolerant and flood intolerant species =====

occupancy_diff <- ggplot(occupancy_abundance_by_point_site, aes(x = Point, y = Occupancy, fill = FloodClass)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.7)) +
  labs(x = "Sampling point", y = "Occupancy (%)") +
  ylim(0, 100) +
  scale_fill_manual(
    values = c("Flood-Sen" = "#66c2a5", "Flood-Tol" = "#fc8d62"),
    labels = c(
      "Flood-Sen" = "<span style='color:#66c2a5'>Flood-Sen</span>",
      "Flood-Tol" = "<span style='color:#fc8d62'>Flood-Tol</span>"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    #legend.position = c(0.77,0.2),
    legend.position ="none",
    panel.grid.minor.y = element_blank(),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = ggtext::element_markdown(size = 10),  # Enable markdown in legend
    legend.title = element_blank()
  )
ggsave("Result/Occupancy difference between flood tolerance groups.tiff", 
       occupancy_diff, width = 4, height = 4, dpi = 600, compression = "lzw")



bubble_plot <- ggplot(occupancy_abundance_by_point_site_species,
                      aes(x = Point, y = Species_colored)) +
  geom_point(aes(size = MeanAbundance, color = Occupancy)) +
  scale_size_area(max_size = 4.5, name = "Mean Abundance (%)") +
  scale_color_viridis_c(option = "D", direction=-1, name = "Occupancy (%)") +
  facet_wrap(~ Site, nrow = 1) +
  labs(x = "Sampling Point", y = "Plant Species") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 9),  # key line for colored labels
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = c(-0.3, 1.62), 
    legend.direction = "vertical",          # make it stack vertically
    legend.box = "vertical",                # allow legends to sit above/below
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(0.4, "cm")      # space between color and size legends
  )

# ===== Final combined plot =====
# Stack env_plot and occupancy_diff side-by-side
top_row <- env_plot | occupancy_diff

# Combine with bubble_plot on top
final_plot <- top_row / bubble_plot +
  plot_layout(heights = c(1.5, 2.2)) +  # adjust height ratio
  plot_annotation(tag_levels = 'a')     # panel tags a, b, c
  final_plot
#Save
ggsave("Result/Occupancy_abundance by site.tiff", final_plot, width = 5, height = 4, dpi = 600, compression = "lzw")


# test the difference in occupancy across point
library(lme4)
library(lmerTest)
library(DHARMa)

occupancy_abundance_by_point_site$Point<-
  factor(occupancy_abundance_by_point_site$Point)
occupancy_abundance_by_point_site$FloodClass<-
  factor(occupancy_abundance_by_point_site$FloodClass)


m1 <- lmer(
  Occupancy ~ FloodClass * Point +
    (1|Site),
  data = occupancy_abundance_by_point_site
)
anova(m1)
summary(m1)

sim <- simulateResiduals(m1)

plot(sim)
testUniformity(sim)  # Tests distribution of residuals (normality proxy)
testDispersion(sim)  # Tests over/under-dispersion (variance issues)
testOutliers(sim)
testResiduals(sim)  # Combines uniformity + dispersion tests

model <- lm(
  Occupancy ~ FloodClass * Point,
  data = occupancy_abundance_by_point_site_species
)

##Test assumptions
shapiro_test <- shapiro.test(residuals(model)) #failed for the normality
print(shapiro_test)

## --- (2) Homogeneity of variance ---
levene_test <- leveneTest(Occupancy ~ FloodClass * Point,
                          data = occupancy_abundance_by_point_site_species) #variance is homogeneous
print(levene_test)

# --- 1. Log transformation (add constant to avoid log(0))
occupancy_abundance_by_point_site_species$log_Occupancy <- 
  log(occupancy_abundance_by_point_site_species$Occupancy+1)

# --- 2. Refit model
log_model <- lm(log_Occupancy ~ FloodClass * Point, 
                data = occupancy_abundance_by_point_site_species)

# --- 3. Recheck assumptions
shapiro.test(residuals(log_model)) #log_occupancy failed the normality test
leveneTest(log_Occupancy ~ FloodClass * Point, 
           data = occupancy_abundance_by_point_site_species)
qqnorm(resid(log_model)); qqline(resid(log_model), col = "red")


# Fit ART model
# Convert to factors (if not already)
occupancy_abundance_by_point_site_species$FloodClass <- 
  factor(occupancy_abundance_by_point_site_species$FloodClass)

occupancy_abundance_by_point_site_species$Point <- 
  factor(occupancy_abundance_by_point_site_species$Point)

# Now run ART
library(ARTool)
art_model <- art(Occupancy ~ FloodClass * Point, 
                 data = occupancy_abundance_by_point_site_species)
art_model <- art(Occupancy ~ FloodClass * Point, 
                 data = occupancy_abundance_by_point_site_species)
anova(art_model)   # Check main and interaction effects

# Post-hoc: if interaction is NOT significant
art_flood <- artlm(art_model, "FloodClass")
cld(emmeans(art_flood, ~ FloodClass), Letters = letters, adjust = "sidak")

art_point <- artlm(art_model, "Point")
cld(emmeans(art_point, ~ Point), Letters = letters, adjust = "sidak")

# Post-hoc: if interaction IS significant
art_interaction <- artlm(art_model, "FloodClass:Point")
cld(emmeans(art_interaction, ~ FloodClass | Point),
    Letters = letters, adjust = "sidak")
