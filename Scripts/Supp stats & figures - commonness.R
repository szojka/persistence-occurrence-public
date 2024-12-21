##############################################
# How common were are species?
##############################################

# To make histograms of how many of the plots have the species in each abundance category.

library(tidyverse)
library(stringr)
library(car)
library(lme4)
library(glmmTMB)
library(visreg)
library(patchwork)
library(ggeffects)
library(DHARMa)

source("Scripts/Source fitness data.R")

############################

# relevant dataframe: green_all
# view(green_all) # abundance = ab_cat, plot = plot, species = species

# save as different dataframe so I can change level names without messing up other graphs:
df_common <- filter(green_all, treatment %in% 'B')

df_common$species <- case_match(df_common$species, 
                                'brohor' ~ "Bromus",
                                'miccal' ~ "Micropus",
                                'plaere' ~ 'Plantago',
                                'vulmic' ~ "Festuca")

# df_common$ab_cat <- case_match(df_common$ab_cat, 
#                                 0 ~ 0,
#                                 1 ~ 1,
#                                 2 ~ 10,
#                                 3 ~ 100)
#                                 

common_fig <- ggplot(df_common) +
  geom_histogram(aes(x = ab_cat, fill = species),
                 position = "dodge",
                 binwidth = 0.3,
                 color = "grey43") +
  labs(x = "Abundance category", y = "Number of blocks", fill = "Species", tag = "B") +
  scale_fill_manual(values = c("slategray1", "skyblue", "lightslateblue", "slateblue4")) + 
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.position = c(0.85,0.8)) 
common_fig
# 
# pdf(
#   'Figures/common_fig.pdf',
#   width = 5,
#   height = 5
# )
# common_fig
# dev.off()
# 
# png(
#   'Figures/common_fig.png',
#   width = 5,
#   height = 5,
#   units = "in",
#   res = 600
# )
# common_fig
# dev.off()


####################################
# Proportion occupancy in each plot

df_occ_prop <- df_common %>%
  group_by(species) %>%
  mutate(tot_n = n()) %>%
  group_by(species, occurrence) %>%
  mutate(prop_occ = n()/tot_n) %>%
  select(prop_occ,species,occurrence) %>%
  filter(occurrence == 1) %>%
  distinct()

# order species
df_occ_prop$species <- factor(df_occ_prop$species, levels = c("Festuca", "Micropus","Plantago","Bromus"))

occ_fig <- ggplot(df_occ_prop, aes(x = species, y = prop_occ, group = species)) +
  geom_point(size = 3,
             color = "skyblue3") +
  theme_classic() +
  labs(x = "Species rank", y = "Proportion of blocks occupied", group = "", tag = "A") +
  theme(text = element_text(size = 16)) +
  ylim(0,1) 
occ_fig

pdf(
  'Figures/common_fig.pdf',
  width = 10,
  height = 5
)
occ_fig + common_fig
dev.off()

png(
  'Figures/common_fig.png',
  width = 10,
  height = 5,
  units = "in",
  res = 600
)
occ_fig + common_fig
dev.off()
