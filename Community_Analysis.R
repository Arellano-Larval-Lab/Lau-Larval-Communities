########################################################
# Larval Community Analysis – Lau Basin Hydrothermal Vents
# Dataset: NEW_Meroplankton_Counts_10.27.25.csv
# Last edited: 2025-11-18
# Author: Vanessa Jimenez
########################################################

#---------------------------------------------------------
# 0. Load Packages
#---------------------------------------------------------
# Install once if needed
  # May not need all of these 
library(vegan)
library(ggplot2)
library(ggforce)
library(mvabund)
library(ape)
library(pairwiseAdonis)
library(dplyr)
library(tidyverse)
library(GGally)
library(ggsci)
library(ggpubr)

#----------------------------------------------------------
# 1. Load Data
#----------------------------------------------------------
setwd("/Users/Nessa/Desktop")
data <- read.csv('NEW_Meroplankton_Counts_10.27.25.csv')
str(data)

#-----------------------------------------------------------
# 2. Data Cleaning & Factor Setup
#-----------------------------------------------------------
# Remove Mata Tolu samples (not included in this analysis)
data <- data[data$RecoveryDive != 'J21416',]
# Remove larval traps 
data <- data[data$SampleMethod != 'Larval Trap',]

data$larvae.30m3 <- (data$Individuals.per.SampleID / data$Sample.Volume.m3) * 30

# Clean and set factor levels
data <- data %>%
  mutate(
    # Factors
    Larval.Type = trimws(Larval.Type),
    Site = factor(Site, levels = c('Kilo Moana', 'Tow Cam', 'Tahi Moana', 'ABE', 'Tui Malila', 'Mariner')),
    Larval.Type = factor(Larval.Type, levels = c('Actinula', 'Bivalve', 'Cyprid', 'Gastropod', 'Nauplius', 'Nectochaete', 'Ophiuroid', 'Parenchymula', 'Trochophore', 'Zoea', 'Unknown')),
    SampleLocation = factor(SampleLocation, levels = c('Within vents', 'Above vents', 'Vent periphery')),
    Ash.depth..cm. = factor(Ash.depth..cm., levels =c('0', '8-15', '12-23', '23-54', '80-150')),
    SampleMethod = factor(SampleMethod, levels = c('McLane', 'SyPRID - above vents', 'SyPRID - vent periphery')),
    RecoveryDive = factor(RecoveryDive),
    larvae.30m3 = as.numeric(larvae.30m3)
  )

str(data)

#--------------------------------------------------------------
# 3. Wide Format for Meroplankton Analysis
#--------------------------------------------------------------
data.wide <- data %>%
  aggregate(larvae.30m3 ~ Larval.Type +RecoveryDive + Site + SampleLocation + Ash.depth..cm., sum)
data.wide <- reshape(data = data.wide,
                     timevar = 'Larval.Type',
                     idvar = c('RecoveryDive', 'Site', 'SampleLocation', 'Ash.depth..cm.'),
                     direction = 'wide')
data.wide[is.na(data.wide)] <- 0
str(data.wide)


#----------------------------------------------------------------
# 4. Ordination & PERMANOVA
#----------------------------------------------------------------
# Bray-Curtis distance
bray_mero <- vegdist(data.wide[, 6:ncol(data.wide)], method = "bray")

# PERMANOVAs: test effect of Site, SampleLocation, Ash depth
set.seed(127)
perm_Site.SampleLocation <- adonis2(data.wide[, 6:ncol(data.wide)] ~ Site * SampleLocation, data = data.wide, method = "bray", permutations = 999)
print(perm_Site.SampleLocation)
# Site p = 0.016 *
# SampleLocation p = 0.001 ***
# Site:SampleLocatoion p = 0.261

perm_Ash.SampleLocation <- adonis2(data.wide[, 6:ncol(data.wide)] ~ Ash.depth..cm. * SampleLocation, data = data.wide, method = "bray", permutations = 999)
print(perm_Ash.SampleLocation)
# Ash depth p = 0.006 **
# SampleLocation p = 0.001***
# Ash depth : SampleLocation p = 0.075

perm_Ash.Site<- adonis2(data.wide[, 6:ncol(data.wide)] ~ Ash.depth..cm. * Site, data = data.wide, method = "bray", permutations = 999)
print(perm_Ash.Site)
# Ash depth p = 0.268
# Site p = 0.974
# no interaction

# Pairwise comparisons 
pairwise_SampleLocation <- pairwise.adonis(data.wide[, 6:ncol(data.wide)], data.wide$SampleLocation)
print(pairwise_SampleLocation)
# Within vents vs above vents p = 0.001
# Within vents vs vent periphery p = 0.001
# Above vents vs vent periphery p = 0.284

pairwise_Ash <- pairwise.adonis(data.wide[, 6:ncol(data.wide)], data.wide$Ash.depth..cm.)
print(pairwise_Ash)

pairwise <- pairwise.adonis(data.wide[, 6:ncol(data.wide)], data.wide$Site)
print(pairwise_Site)
# No significant differences 


# Beta dispersion
beta_Site <- betadisper(bray_mero, data.wide$Site)
anova(beta_Site) # p = 0.9784
plot(beta)

beta_Ash <- betadisper(bray_mero, data.wide$Ash.depth..cm.)
anova(beta_Ash) # p = 0.9941
plot(beta_Ash)

beta_SampleLocation <- betadisper(bray_mero, data.wide$SampleLocation)
anova(beta_SampleLocation) # p = 0.2497
plot(beta_SampleLocation)

# SIMPER analysis
simper_SampleLocation <- simper(data.wide[, 6:ncol(data.wide)], data.wide$SampleLocation, permutations = 999)
summary(simper_SampleLocation)

simper_Ash <- simper(data.wide[, 6:ncol(data.wide)], data.wide$Ash.depth..cm., permutations = 999)
summary(simper_Ash)

simper_Site <- simper(data.wide[, 6:ncol(data.wide)], data.wide$Site, permutations = 999)
summary(simper_Site)

#----------------------------------------------------------------
# 5. PCA
#----------------------------------------------------------------
pca_mer <- prcomp(data.wide[, 6:ncol(data.wide)], scale = TRUE)
summary(pca_mer)
biplot(pca_mer)

# Eigenvalue scree plot
prop_var <- eigenvals(pca_mer)/sum(eigenvals(pca_mer))
plot(prop_var, type = "l", ylab = "Proportion of Variance Explained", xlab = "PC")

pca_df <- data.frame(
  Site = data.wide$Site,
  SampleLocation = data.wide$SampleLocation,
  PC1 = pca_mer$x[,1],
  PC2 = pca_mer$x[,2]
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Site, shape = SampleLocation)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "PCA of Larval Communities")

#-----------------------------
# 6. PCoA
#-----------------------------
pcoa_mer <- pcoa(bray_mero)

# Hierarchical clustering into three groups 
hc <- hclust(bray_mero, method = "ward.D2")
plot(as.dendrogram(hc))
rect.hclust(hc, k = 3, border = "salmon")

clusters <- cutree(hc, k = 3)
clusters <- factor(clusters, levels = 1:3)

pcoa_df <- data.frame(
  Axis1 = pcoa_mer$vectors[,1],
  Axis2 = pcoa_mer$vectors[,2],
  Cluster = as.factor(clusters),
  Site = data.wide$Site,
  SampleLocation = data.wide$SampleLocation
)

pal_site <- c('#AAAA00','#44BB99','#99DDFF','#FFAABB','#EEDD88','#EE8866','#77AADD','#BBCC33', '#CC6677')

# Final PCOA Figure
ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Site, shape = SampleLocation)) +
  geom_point(size = 4) +
  geom_mark_ellipse(aes(fill = Cluster, group = Cluster), alpha = 0.15, show.legend = FALSE) +
  scale_color_manual(values = pal_site) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  labs(x = "Axis 1", y = "Axis 2") 

#-----------------------------
# 7. Diversity Indices
#-----------------------------
diversity_results <- data %>%
  filter(SampleMethod != "Larval Traps") %>%
  group_by(Site, Meters.Above.Bottom) %>%
  summarise(
    Species_richness = sum(c_across(starts_with("larvae.30m3")) > 0),
    Shannon_Index = diversity(c_across(starts_with("larvae.30m3"))),
    Simpson_Index = diversity(c_across(starts_with("larvae.30m3")), index = "simpson"),
    .groups = "drop"
  )

# Shannon PERMANOVA (example)
perm_div <- adonis2(Shannon_Index ~ Meters.Above.Bottom, data = diversity_results, permutations = 999)
print(perm_div)

# Boxplots
ggplot(diversity_results, aes(x = Site, y = Shannon_Index, fill = Site)) +
  geom_boxplot() +
  theme_minimal() +
  labs(y = "Shannon Index")

ggplot(diversity_results, aes(x = Site, y = Species
                              

#-----------------------------------------------------------------
# 8. Stacked Bar Plots: Total Meroplankton by Larval Type
#-----------------------------------------------------------------
pal <- c(
  '#AAAA00', '#44BB99', '#99DDFF', '#FFAABB', '#EEDD88',
  '#EE8866', '#77AADD', '#DDCCFF', '#BBCC33', '#117733',
  '#332288'
)

make_stacked_bar <- function(data, SampleLocation, 
                             yvar = "larvae.30m3", 
                             fillvar = "Larval.Type",
                             groupvar = "RecoveryDive") {
  
  ggplot(data,
         aes(x = Site,
             y = .data[[yvar]],
             fill = .data[[fillvar]],
             group = .data[[groupvar]])) +
    
    geom_bar(position = position_dodge(width = 0.9),
             stat = "identity") +
    
    labs(
      x = "",
      y = expression('Larvae per 30' ~ m^3),
      fill = "Larval Type",
      title = SampleLocation
    ) +
    
    theme_classic() +
    theme(
      text = element_text(size = 10),
      axis.text.x = element_text(color = "black")
    ) +
    
    scale_fill_manual(values = pal, drop = FALSE)
}

data_within  <- subset(data, SampleLocation == 'Within vents', drop = FALSE)
data_above  <- subset(data, SampleLocation == 'Above vents', drop = FALSE)
data_periph <- subset(data, SampleLocation == 'Vent periphery', drop = FALSE)

within <- make_stacked_bar(data = data_within,SampleLocation = "Within vents")
above <- make_stacked_bar(data = data_above,SampleLocation = "Above vents")
periph <- make_stacked_bar(data = data_periph,SampleLocation = "Vent periphery",yvar = "larvae.30m3",fillvar = "Larval.Type")

fig <- ggarrange(
  periph + theme(axis.text.x = element_blank(), axis.title.y = element_blank()),
  above  + theme(axis.text.x = element_blank(), axis.title.y = element_blank()),
  within  + theme(axis.title.y = element_blank()),
  ncol = 1, nrow = 3,
  common.legend = TRUE,
  legend = "right")

fig <- annotate_figure(
  fig,
  left = text_grob(expression('Larvae per 30' ~m^3), rot = 90))

fig


#-------------------------------------------------------------------------
# 9. Meroplankton Relative Abundance ***Can't seem to separate by RecoveryDive and keep relative abundances. 
#-------------------------------------------------------------------------
mero_long <- data %>%
  filter(SampleMethod != "Larval Traps") %>%
  group_by(Site, SampleLocation, Larval.Type, RecoveryDive) %>%
  summarise(larvae.30m3 = sum(larvae.30m3), .groups = "drop") %>%
  group_by(Site, SampleLocation, RecoveryDive) %>%
  mutate(RelAbund = larvae.30m3 / sum(larvae.30m3)) %>%
  ungroup() %>%
  mutate(
    Site2 = factor(Site, levels = c('Kilo Moana','Tow Cam','Tahi Moana','ABE','Tui Malila','Mariner'),
                   labels = c('KM','TC','TM','ABE','TuM','Mar'))
  )

ggplot(mero_long, aes(x = Site2, y = RelAbund, fill = Larval.Type, group = RecoveryDive)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(SampleLocation~.) +
  scale_fill_manual(values = pal) +
  theme_classic()


#---------------------------------------------------------------------------
# 7. Gastropod Analysis and Plots by Taxonomic Group
#---------------------------------------------------------------------------  
gastro_data <- subset(data, Larval.Type == 'Gastropod')

gastro_data <- gastro_data %>%
  mutate(
    Taxonomic.Group = trimws(Taxonomic.Group),
    Taxonomic.Group = factor(Taxonomic.Group, levels = c('Gastropoda', 'Alviniconcha', 'Bruceiella', 'Buccinoidea', 'Cerithium', 'Desbruyeresia', 'Elachisinidae', 'Lamellomphalus', 'Mitrella', 'Neoleptopsidae', 'Neomphalida', 'Peltospiridae', 'Raphitomidae', 'Shinkailepas', 'Symmetromphalus', 'Thermosipho')))
str(gastro_data)

pal <- c("#E69F00", "#F0E442","#56B4E9", "#A6CEE3", "#B2DF8A","#99CC99", "#66C2A5","#8DD3C7",
  "#FDBF6F","#FFCC99","#F4A582", "#FB9A99", "#CAB2D6", "#DDCCFF", "#FFD92F", "#CCEBC5")

gastro_plot <- function(data, SampleLocation, 
                             yvar = "larvae.30m3", 
                             fillvar = "Taxonomic.Group",
                             groupvar = "RecoveryDive") {
  
  ggplot(data,
         aes(x = Site,
             y = .data[[yvar]],
             fill = .data[[fillvar]],
             group = .data[[groupvar]])) +
    
    geom_bar(position = position_dodge(width = 0.9),
             stat = "identity") +
    
    labs(
      x = "",
      y = expression('Larvae per 30' ~ m^3),
      fill = "Gastropod Group",
      title = SampleLocation
    ) +
    
    theme_classic() +
    theme(
      text = element_text(size = 10),
      axis.text.x = element_text(color = "black")
    ) +
    
    scale_fill_manual(values = pal, drop = FALSE)
}


gastro_within  <- subset(gastro_data, SampleLocation == 'Within vents', drop = FALSE)
gastro_above  <- subset(gastro_data, SampleLocation == 'Above vents', drop = FALSE)
gastro_periph <- subset(gastro_data, SampleLocation == 'Vent periphery', drop = FALSE)

within.gastro <- gastro_plot(data = gastro_within,SampleLocation = "Within vents")

above.gastro <- gastro_plot(data = gastro_above,SampleLocation = "Above vents")

periph.gastro <- gastro_plot(data = gastro_periph,SampleLocation = "Vent periphery",yvar = "larvae.30m3",fillvar = "Taxonomic.Group")

fig.gastro <- ggarrange(
  periph.gastro + theme(axis.text.x = element_blank(), axis.title.y = element_blank()),
  above.gastro  + theme(axis.text.x = element_blank(), axis.title.y = element_blank()),
  within.gastro  + theme(axis.title.y = element_blank()),
  ncol = 1, nrow = 3,
  common.legend = TRUE,
  legend = "right"   # or "right"
)

fig.gastro <- annotate_figure(
  fig.gastro,
  left = text_grob(expression('Larvae per 30' ~m^3), rot = 90)
)
fig.gastro

# Gastro Larval Traps 

#------------------------------------------------------------------------
# 7. Bivalve Plots 
#-------------------------------------------------------------------------
bivalve_data <- subset(data, Larval.Type == 'Bivalve')
str(bivalve_data) 

bivalve_data <- bivalve_data %>%
  mutate(Taxonomic.Group = factor(Taxonomic.Group, levels = c('Bathymodiolus', 'Bivalvia', 'Galeommatoidae', 'Lasaeidae', 'Myidae', 'Nuculanida', 'Pectinidae', 'Thyasiridae', 'Vesicomyidae')))

pal <- c('#00cec9','#fab1a0', '#74b9ff','#ff7675', '#a29bfe','#ffeaa7', '#00b894','#e17055', "#B3DE69", '#fdcb6e', '#a49bfe', '#ff7575' )

# Function to generate stacked bar plots by Altitude
plot_bivalve <- function(data, SampleLocation,
                             yvar = "larvae.30m3", 
                             fillvar = "Taxonomic.Group",
                             groupvar = "RecoveryDive") {
  
  ggplot(data,
         aes(x = Site,
             y = .data[[yvar]],
             fill = .data[[fillvar]],
             group = .data[[groupvar]])) +
    
    geom_bar(position = position_dodge(width = 0.9),
             stat = "identity") +
    
    labs(
      x = "",
      y = expression('Larvae per 30' ~ m^3),
      fill = "Taxonomic Group",
      title = SampleLocation
    ) +
    
    theme_classic() +
    theme(
      text = element_text(size = 10),
      axis.text.x = element_text(color = "black")
    ) +
    
    scale_fill_manual(values = pal, drop = FALSE)
}


bivalve_within  <- subset(bivalve_data, SampleLocation == 'Within vents', drop = FALSE)
bivalve_above  <- subset(bivalve_data, SampleLocation == 'Above vents', drop = FALSE)
bivalve_periph <- subset(bivalve_data, SampleLocation == 'Vent periphery', drop = FALSE)

within.bivalve <- plot_bivalve(data = bivalve_within,SampleLocation = "Within vents")

above.bivalve <- plot_bivalve(data = bivalve_above,SampleLocation = "Above vents")

periph.bivalve <- plot_bivalve(data = bivalve_periph,SampleLocation = "Vent periphery",yvar = "larvae.30m3",fillvar = "Taxonomic.Group")

fig.bivalve <- ggarrange(
  periph.bivalve + theme(axis.text.x = element_blank(), axis.title.y = element_blank()),
  above.bivalve  + theme(axis.text.x = element_blank(), axis.title.y = element_blank()),
  within.bivalve  + theme(axis.title.y = element_blank()),
  ncol = 1, nrow = 3,
  common.legend = TRUE,
  legend = "right"   # or "right"
)

fig.bivalve <- annotate_figure(
  fig.bivalve,
  left = text_grob(expression('Larvae per 30' ~m^3), rot = 90)
)
fig.bivalve

# Bivalve Relative Abundances

# Bivalve Larval Traps
