
#----------------------------------------------------------------------.
# PERSISTENCE VS OCCURRENCE VS ABIOTIC PERSISTENCE ACROSS GRADIENT ####
#----------------------------------------------------------------------.

library(tidyverse)
library(stringr)
library(car)
library(lme4)
library(glmmTMB)
library(visreg)
library(patchwork)
library(ggeffects)
library(DHARMa)
library(emmeans)

source("Scripts/Source fitness data.R")

# DESCRIPTION:
# code for the data wrangling, modelling, and visualizing community level 
# responses across the productivity gradient.

# NOTE:
# we fit abiotic persistence, biotic persistence, and occurrence all together 
# so that the shared factor--biotic persistence--retained the same model structure.

##################
## Prepare data ##
##################

# note for figure 2, comparing occurrence and persistence, persistence in this figure is measured with neighbors, as this is indicative of natural conditions.

# combine O, PA, PB into one variable
pb <- fig_dat %>% 
  filter(treatment %in% c('B') & type %in% c('P')) %>%
  mutate(type = "PB") %>%
  select(-treatment)
pa <- fig_dat %>% 
  filter(treatment %in% c('A') & type %in% c('P')) %>%
  mutate(type = "PA") %>%
  select(-treatment)
ob <- fig_dat %>% 
  filter(treatment %in% c('B') & type %in% c('O')) %>%
  mutate(type = "OB") %>%
  select(-treatment)
fig_dat_fin <- rbind(pb,pa,ob)

hist(fig_dat_fin$species_no)

fig_dat_fin <- select(fig_dat_fin, grid, site, block, type, species_no, green_index_scaled, replicates)

# make response variable proportions
fig_dat_fin <- fig_dat_fin %>%
  mutate(species_prop = species_no/replicates)

# make type a factor
fig_dat_fin$type <- as.factor(fig_dat_fin$type)

# Sample sizes
fig_dat_fin %>%
  filter(replicates>0) %>%
  mutate(frequency = n()) %>%
  select(frequency) %>%
  distinct() #1285

# N block
fig_dat_fin %>%
  select(block) %>%
  distinct(block) %>%
  mutate(frequency = n()) %>%
  select(frequency) %>%
  distinct() # 435 blocks

#####################
## Model selection ##
#####################

# Notes:
# In order to use a vector of proportions as the response variable with glmer (., family = binomial), 
# you need to set the number of trials that led to each proportion using the weights argument. 
# For example, using the cbpp data from the lme4 package: 
# glmer (incidence / size ~ period + (1 | herd), weights = size, family = binomial, data = cbpp)
# glmmTMB converges, glmer did not converge (isSingular)
# replicates in fig_dat let's us know how many species are viable data within each plot

# random effects are equivalent to ... + (1|site/grid/block)
m2 <-
  glmmTMB(species_prop ~ type * green_index_scaled + (1|block) + (1|site) + (1|site:grid),
        data= fig_dat_fin,
        family = binomial(),
        weights = replicates) # the weights are the number of trials used to generate each proportion

m2q <-  glmmTMB(species_prop ~ type * poly(green_index_scaled,2,raw=TRUE) + (1|block) + (1|site) + (1|site:grid), 
              data = fig_dat_fin,
              family = binomial(),
              weights = replicates)
  

AIC(m2, m2q)
# df      BIC
# m2   9 3073.549
# m2q 12 3091.261

BIC(m2, m2q)
# df      BIC
# m2   9 3073.549
# m2q 12 3091.261

#####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(m2)
testDispersion(m2) 
plot(fitted(m2), residuals(m2))
hist(residuals(m2)) 
loc_sim_ouput <- simulateResiduals(m2)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# m2 and m2q no sign issues

#############################
## Interpretation of coefs ##
#############################

Anova(m2, type = 3) #type:green_index 36.9901  1  1.187e-09 ***
# Anova is formal test, can usually just report statistics from summary (estimate +/- SE, p-value)
# best way to code in ANOVA is car package - type 3 ANOVA

# results show us the relationsip between species number and vegetation index depends on treatment (persistence or occurrence)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
#Response: species_prop
# Chisq Df Pr(>Chisq)    
# (Intercept)               4.8826  1    0.02713 *  
#   type                    397.8593  2  < 2.2e-16 ***
#   green_index_scaled       57.7232  1  3.017e-14 ***
#   type:green_index_scaled  84.9378  2  < 2.2e-16 ***



summary(m2) # type O is baseline (occurrence)
# 
#Conditional model:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                  0.20686    0.09362   2.210   0.0271 *  
#   typePA                    -0.91405    0.07914 -11.550  < 2e-16 ***
#   typePB                    -1.65271    0.08407 -19.660  < 2e-16 ***
#   green_index_scaled         0.56852    0.07483   7.598 3.02e-14 ***
#   typePA:green_index_scaled -0.39752    0.08769  -4.533 5.81e-06 ***
#   typePB:green_index_scaled -0.81888    0.08889  -9.213  < 2e-16 ***

#############################
## Estimate contrasts ##
#############################

emtrends(m2, ~ type | green_index_scaled, "green_index_scaled", component="cond")
# green_index_scaled = 0.00401, degree = linear:
# type    green_index_scaled.trend  SE  df  asymp.LCL asymp.UCL
# OB                      0.569 0.0748 Inf    0.4219    0.7152
# PA                      0.171 0.0763 Inf    0.0215    0.3205
# PB                     -0.250 0.0785 Inf   -0.4043   -0.0964

em <- emmeans(m2, ~type | green_index_scaled, at = list(green_index_scaled = c(-4, -1, 2))) # specify green_index_scaled values
pairs(em)
# green_index_scaled = -4:
#   contrast estimate    SE  df z.ratio p.value
# OB - PA    -0.676 0.364 Inf  -1.856  0.1517
# OB - PB    -1.623 0.363 Inf  -4.467  <.0001
# PA - PB    -0.947 0.374 Inf  -2.532  0.0304
# 
# green_index_scaled = -1:
#   contrast estimate    SE  df z.ratio p.value
# OB - PA     0.517 0.122 Inf   4.245  0.0001
# OB - PB     0.834 0.121 Inf   6.904  <.0001
# PA - PB     0.317 0.125 Inf   2.529  0.0307
# 
# green_index_scaled =  2:
#   contrast estimate    SE  df z.ratio p.value
# OB - PA     1.709 0.188 Inf   9.094  <.0001
# OB - PB     3.290 0.199 Inf  16.571  <.0001
# PA - PB     1.581 0.202 Inf   7.812  <.0001

#contrast(em, "trt.vs.ctrl", ref = "OB", type = "response", adjust = "mvt")

###############
## Visualize ##
###############

# All together three lines:

vis_prod <- ggpredict(m2, 
                      terms = c("green_index_scaled[all]", "type"), 
                      type = "fe", allow.new.levels=TRUE); plot(vis_prod)

length(vis_prod$x) # 1152/3 = 384
colortreatpred <- c(rep("tan4", times = 384), rep("maroon4", times = 384), rep("turquoise1", times = 384))

# all lines together:
ggplot() +
  geom_jitter(
    data = fig_dat_fin,
    aes(x = green_index_scaled, y = species_prop, color = type),
    alpha = 0.3,
    width = 0.01,
    height = 0.1,
    size = 1.5,
  #  color = colortreatraw
  )  +
  scale_color_manual(values =c("tan4","maroon4","turquoise1"), 
                     labels = c("Occurrence","Persistence without neighbors", "Persistence with neighbors")) + # up here so that it only deals with this layer
  geom_line(data = vis_prod, aes(x = x, y = predicted, group = group, color = group),
           # color = colortreatpred,
            linewidth = 1
  ) +
  geom_ribbon(data = vis_prod, aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high),
    alpha = 0.2,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() + # should always be above other themes
  theme(legend.title = element_blank(),
        legend.position = c(0.2,0.85)) + 
  labs(y = "Proportion of species", x = "Vegetation index (G-R)") 

######################################
# separate three lines into two panels:
######################################

temp <- filter(vis_prod, group %in% c('OB', "PB"))
length(temp$x) # 768/2 = 384
#colors <- c("tan4","maroon4","springgreen")
colortreatpred <- c(rep("tan4", times = 384), rep("turquoise1", times = 384))


# a = occurrence vs persistence (ob and pb)
a <- 
  ggplot() +
  geom_jitter(
    data = filter(fig_dat_fin, type %in% c('OB', "PB")),
    aes(x = green_index_scaled, y = species_prop, color = type),
    alpha = 0.2,
    width = 0.01,
    height = 0.05,
    size = 1.5,
    #  color = colortreatraw
  )  +
  scale_color_manual(values =c("tan4","turquoise1"), 
                     labels = c("Occurrence", "Persistence")) + # up here so that it only deals with this layer
  geom_line(data = filter(vis_prod, group %in% c('OB', "PB")), aes(x = x, y = predicted, group = group, color = group),
            # color = colortreatpred,
            linewidth = 1
  ) +
  geom_ribbon(data = filter(vis_prod, group %in% c('OB', "PB")), aes(
    x = x,
    y = predicted,
    group=group,
    ymin = conf.low,
    ymax = conf.high),
    alpha = 0.25,
    show.legend = F,
   fill = colortreatpred
  ) +
  theme_classic() + # should always be above other themes
  theme(legend.title = element_blank(),
        legend.position = c(0.2,0.85),
        text = element_text(size = 16)) + 
  labs(y = "Proportion of species", x = "Vegetation index (G-R)") 


# b = persistence aboiotic vs biotic (pa vs pb)
temp <- filter(vis_prod, group %in% c('PA', "PB"))
length(temp$x) # 768/2 = 384
colortreatpred <- c(rep("maroon4", times = 384), rep("turquoise1", times = 384))


# b = pa vs pb (comparing persistence in abiotic versus biotic)
b <- 
  ggplot() +
  geom_jitter(
    data = filter(fig_dat_fin, type %in% c('PA', "PB")),
    aes(x = green_index_scaled, y = species_prop, color = type),
    alpha = 0.2,
    width = 0.01,
    height = 0.05,
    size = 1.5,
    #  color = colortreatraw
  )  +
  scale_color_manual(values =c("maroon4","turquoise1"), 
                     labels = c("Persistence without neighbors", "Persistence with neighbors")) + # up here so that it only deals with this layer
  geom_line(data = filter(vis_prod, group %in% c('PA', "PB")), aes(x = x, y = predicted, group = group, color = group),
            # color = colortreatpred,
            linewidth = 1
  ) +
  geom_ribbon(data = filter(vis_prod, group %in% c('PA', "PB")), aes(
    x = x,
    y = predicted,
    group=group,
    ymin = conf.low,
    ymax = conf.high),
    alpha = 0.25,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() + # should always be above other themes
  theme(legend.title = element_blank(),
        legend.position = c(0.4,0.85),
        text = element_text(size = 16)) + 
  labs(y = "Proportion of species", x = "Vegetation index (G-R)") 

#-------------------------------

pdf(
  'Figures/fig_persist_vs_occur.pdf',
  width = 6,
  height = 4
)
a
dev.off()

png(
  'Figures/fig_persist_vs_occur.png',
  width = 6,
  height = 4,
  units = 'in',
  res = 600
)
a
dev.off()

#----------------------------

pdf(
  'Figures/fig_biotic_interactions.pdf',
  width = 6,
  height = 4
)
b
dev.off()

png(
  'Figures/fig_biotic_interactions.png',
  width = 6,
  height = 4,
  units = 'in',
  res = 600
)
b
dev.off()
