#------------------------------------------------------------------------
# SPECIES SPECIFIC OCCURRENCE PATTERNS
#------------------------------------------------------------------------

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

# DESCRIPTION:
# for each species (Plantago, Micropus, Bromus, and Vulpia (aka Festuca)) model occurrence patterns across the productivity gradient

# INSTRUCTIONS:
# run before "Stats & figures - spp seed production.R"

# NOTES:
# occurrence is a binary so using logistic binomial to fit

#######################
# PLANTAGO
#######################

occ_plantago <- filter(green_all, species %in% 'plaere' & treatment %in% 'B')
hist(occ_plantago$occurrence) # no zi needed for binomials
mpo <- glmmTMB(occurrence ~ green_index_scaled + (1|site) + (1|grid:site), data = occ_plantago, family = binomial(link = 'logit'))

summary(mpo)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.2081     0.6765   0.308    0.758
# green_index_scaled   0.2592     0.2602   0.996    0.319
length(occ_plantago$grid)

vis_occ_p <-
  ggpredict(mpo,
            terms = c("green_index_scaled[all]"),
            type = "fe",allow.new.levels=TRUE)

f <- ggplot() +
  geom_jitter(
    data = occ_plantago,
    aes(green_index_scaled, occurrence),
    alpha = 0.1,
    size = 1.5,
    color = "midnightblue",
    height = 0.1
  ) +
  geom_line(data = vis_occ_p, aes(x = x, y = predicted),
            color = "midnightblue",
            linewidth = 1) +
  geom_ribbon(
    data = vis_occ_p,
    aes(
      x = x,
      y = predicted,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = "midnightblue"
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "Vegetation index (G-R)") +
  ylim(0, 1) + 
  annotate("text", x = -4, y = 1, label = "F", size = 7)
f

#######################
# MICROPUS
#######################

occ_micropus <- filter(green_all, species %in% 'miccal'& treatment %in% 'B')
mmo <- glmmTMB(occurrence ~ green_index_scaled + (1|site) + (1|grid:site), data = occ_micropus, family = binomial(link = 'logit'))

summary(mmo)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.6644     0.4370   1.520    0.128
# green_index_scaled   0.2936     0.2271   1.293    0.196
length(occ_micropus$grid)

vis_occ_m <-
  ggpredict(mmo,
            terms = c("green_index_scaled[all]"),
            type = "fe", allow.new.levels=TRUE)

e <- ggplot() +
  geom_jitter(
    data = occ_micropus,
    aes(green_index_scaled, occurrence),
    alpha = 0.1,
    size = 1.5,
    color = "midnightblue",
    height = 0.1
  ) +
  geom_line(data = vis_occ_m, aes(x = x, y = predicted),
            color = "midnightblue",
            linewidth = 1) +
  geom_ribbon(
    data = vis_occ_m,
    aes(
      x = x,
      y = predicted,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = "midnightblue"
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "Occurrence", x = "Vegetation index (G-R)") +
  ylim(0, 1)+ 
  annotate("text", x = -4, y = 1, label = "E", size = 7)
e

#######################
# BROMUS
#######################

occ_bromus <- filter(green_all, species %in% 'brohor' & treatment %in% 'B')
mbo <- glmmTMB(occurrence ~ green_index_scaled + (1|site) + (1|grid:site), data = occ_bromus, family = binomial(link = 'logit'))


summary(mbo)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -0.7113     0.9627  -0.739 0.459988    
# green_index_scaled   1.0350     0.3046   3.398 0.000678 ***

length(occ_bromus$grid)

vis_occ_b <-
  ggpredict(mbo,
            terms = c("green_index_scaled[all]"),
            type = "fe",allow.new.levels=TRUE)

g <- ggplot() +
  geom_jitter(
    data = occ_bromus,
    aes(green_index_scaled, occurrence),
    alpha = 0.1,
    size = 1.5,
    color = "midnightblue",
    height = 0.1
  ) +
  geom_line(data = vis_occ_b, aes(x = x, y = predicted),
            color = "midnightblue",
            linewidth = 1) +
  geom_ribbon(
    data = vis_occ_b,
    aes(
      x = x,
      y = predicted,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = "midnightblue"
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "Vegetation index (G-R)") +
  ylim(0, 1) + 
  annotate("text", x = -4, y = 1, label = "G", size = 7)
g

#######################
# VULPIA
#######################

occ_vulpia <- filter(green_all, species %in% 'vulmic'& treatment %in% 'B')
mvo <- glmmTMB(occurrence ~ green_index_scaled + (1|site) + (1|grid:site), data = occ_vulpia, family = binomial(link = 'logit'))

summary(mvo)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)          0.6022     0.4003   1.504   0.1325   
# green_index_scaled   0.5662     0.2186   2.590   0.0096 **

length(occ_vulpia$tag)

vis_occ_v <-
  ggpredict(mvo,
            terms = c("green_index_scaled[all]"),
            type = "fe",allow.new.levels=TRUE)

h <- ggplot() +
  geom_jitter(
    data = occ_vulpia,
    aes(green_index_scaled, occurrence),
    alpha = 0.1,
    size = 1.5,
    color = "midnightblue",
    height = 0.1
  ) +
  geom_line(data = vis_occ_v, aes(x = x, y = predicted),
            color = "midnightblue",
            linewidth = 1) +
  geom_ribbon(
    data = vis_occ_v,
    aes(
      x = x,
      y = predicted,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = "midnightblue"
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "Vegetation index (G-R)") +
  ylim(0, 1) + 
  annotate("text", x = -4, y = 1, label = "H", size = 7)
h
