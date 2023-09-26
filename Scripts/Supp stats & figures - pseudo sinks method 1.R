
#################################
# Supplemental look at density 
# dependence and pseudo-sinks
#################################

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

###############
## METHOD 1 ##
###############

# First, look at seed production ~ abundance category * greeness 
# Visualize looking at different slices of greeness.

green_all$ab_cat_shared <- as.factor(green_all$ab_cat_shared)
# filter for the unmanipulated plots (Biotic)
temp_dat <- filter(green_all, treatment %in% 'B')

# check variables
temp_dat$ab_cat_shared <- as.factor(temp_dat$ab_cat_shared)
temp_dat$species <- as.factor(temp_dat$species)

# change names for plotting
temp_dat$ab_cat_shared <-
  recode_factor(temp_dat$ab_cat_shared, 
                "0" = "0 conspecifics",
                "1" = "1-10 conspecifics",
                "2" = "11-100 conspecifics",
                "3" = "> 101 conspecifics")

################################################################################
# Plantago


#########
# Models
#########

mddp.1 <- glmmTMB(seed_half ~ ab_cat_shared*green_index_scaled + 
                    (1 | site) + (1| grid:site), 
                  data = filter(temp_dat, species %in% c("plaere")),
                  ziformula = ~1,
                  family = poisson(link="log")) 
mddp.2 <- glmmTMB(seed_half ~ poly(green_index_scaled,2, raw = TRUE)*ab_cat_shared + 
                    (1 | site) + (1| grid:site), 
                  data = filter(temp_dat, species %in% c("plaere")),
                  ziformula = ~1,
                  family = poisson(link="log"))
mddp.3 <- glmmTMB(seed_half ~ ab_cat_shared*green_index_scaled + 
                    (1 | site) + (1| grid:site), 
                  ziformula = ~.,
                  data = filter(temp_dat, species %in% c("plaere")),
                  family = poisson(link="log")) 
mddp.4 <- glmmTMB(seed_half ~ poly(green_index_scaled,2, raw = TRUE)*ab_cat_shared + 
                    (1 | site) + (1| grid:site), 
                  ziformula = ~.,
                  data = filter(temp_dat, species %in% c("plaere")),
                  family = poisson(link="log"))

AIC(mddp.1,mddp.2,mddp.3,mddp.4)
# df      AIC
# mddp.1 11 905.2504
# mddp.2 15 912.1818
# mddp.3 20 903.0918 # linear w. zi slope
# mddp.4 28 914.7320

#############
# Check fit
#############

# using DHARMa:
testZeroInflation(mddp.3)
testDispersion(mddp.3) 
plot(fitted(mddp.3), residuals(mddp.3))
hist(residuals(mddp.3)) 
loc_sim_ouput <- simulateResiduals(mddp.3)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)
# all good

########
# Stats
########

summary(mddp.3)
Anova(mddp.3, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: seed_half
# Chisq Df Pr(>Chisq)    
#   (Intercept)                       0.0029  1     0.9571    
#   ab_cat_shared                    57.8229  3  1.715e-12 ***
#   green_index_scaled                1.1764  1     0.2781    
#   ab_cat_shared:green_index_scaled  3.1088  3     0.3752  

# investigate which abundance categories are difference:
emm <- emmeans(mddp.3, ~ ab_cat_shared)
simp <- pairs(emm, simple = "each")
simp
# contrast                                    estimate    SE  df z.ratio p.value
# 0 conspecifics - (1-10 conspecifics)         -0.4319 0.210 Inf  -2.054  0.1685
# 0 conspecifics - (11-100 conspecifics)        1.1180 0.291 Inf   3.840  0.0007 # here
# 0 conspecifics - > 101 conspecifics           0.0606 0.732 Inf   0.083  0.9998
# (1-10 conspecifics) - (11-100 conspecifics)   1.5499 0.208 Inf   7.451  <.0001 # here
# (1-10 conspecifics) - > 101 conspecifics      0.4926 0.711 Inf   0.692  0.9001
# (11-100 conspecifics) - > 101 conspecifics   -1.0573 0.734 Inf  -1.440  0.4742

############
# visualize
############

vis_dd_p <-
  ggpredict(mddp.3,
            terms = c("green_index_scaled[all]","ab_cat_shared"),
            type = "fe.zi", allow.new.levels=TRUE) 
# fe.zi means seed production

ddp <- ggplot() +
  geom_line(data = vis_dd_p, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1
            # color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_dd_p,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high,
      fill = group
    ),
    alpha = 0.1,
    show.legend = F
    #fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 3) +
  ggtitle("Plantago") + 
  labs(x = "Vegetation index (G-R)", color = "Abundance") +
  annotate("text", x = -4, y = 3, label = "B", size = 6)
ddp


################################################################################
# Micropus


#########
# Models
#########

mddm.1 <- glmmTMB(seed_half ~ ab_cat_shared*green_index_scaled + 
                    (1 | site) + (1| grid:site), 
                  data = filter(temp_dat, species %in% c("miccal")),
                  ziformula = ~1,
                  family = poisson(link="log")) 
mddm.2 <- glmmTMB(seed_half ~ poly(green_index_scaled,2, raw = TRUE)*ab_cat_shared + 
                    (1 | site) + (1| grid:site), 
                  data = filter(temp_dat, species %in% c("miccal")),
                  ziformula = ~1,
                  family = poisson(link="log"))
mddm.3 <- glmmTMB(seed_half ~ ab_cat_shared*green_index_scaled + 
                    (1 | site) + (1| grid:site), 
                  ziformula = ~.,
                  data = filter(temp_dat, species %in% c("miccal")),
                  family = poisson(link="log")) 
mddm.4 <- glmmTMB(seed_half ~ poly(green_index_scaled,2, raw = TRUE)*ab_cat_shared + 
                    (1 | site) + (1| grid:site), 
                  ziformula = ~.,
                  data = filter(temp_dat, species %in% c("miccal")),
                  family = poisson(link="log"))

AIC(mddm.1,mddm.2,mddm.3,mddm.4)
# df      AIC
# mddm.1 11 119.4972 # linear w. zi intercept
# mddm.2  9       NA
# mddm.3 20       NA
# mddm.4 16       NA

#############
# Check fit
#############

# using DHARMa:
testZeroInflation(mddm.1)
testDispersion(mddm.1) 
plot(fitted(mddm.1), residuals(mddm.1))
hist(residuals(mddm.1)) 
loc_sim_ouput <- simulateResiduals(mddm.1)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)
# all good

########
# Stats
########

summary(mddm.1)
Anova(mddm.1, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Chisq Df Pr(>Chisq)
# (Intercept)                      1.0967  1     0.2950
# ab_cat_shared                    1.0088  3     0.7991
# green_index_scaled               0.8646  1     0.3524
# ab_cat_shared:green_index_scaled 0.0079  3     0.9998   

##############
# visualize
##############

vis_dd_m <-
  ggpredict(mddm.1,
            terms = c("green_index_scaled[all]","ab_cat_shared"),
            type = "fe.zi", allow.new.levels=TRUE) 
# fe.zi means seed production

plot(vis_dd_m)

ddm <- ggplot() +
  geom_line(data = vis_dd_m, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1
            # color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_dd_m,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high,
      fill = group
    ),
    alpha = 0.1,
    show.legend = F
    #fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.6,0.7),
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 3) +
  ggtitle("Micropus") + 
  labs(x = "Vegetation index (G-R)", y= "Seed production", color = "Abundance") +
  annotate("text", x = -4, y = 3, label = "A", size = 6)
ddm

################################################################################
# Bromus


#########
# Models
#########


# simplified random effects to fix singular fit
mddb.1 <- glmmTMB(seed_half ~ ab_cat_shared*green_index_scaled + 
                    (1 | site),# + (1| grid:site), 
                  data = filter(temp_dat, species %in% c("brohor")),
                  ziformula = ~1,
                  family = poisson(link="log")) 
mddb.2 <- glmmTMB(seed_half ~ poly(green_index_scaled,2, raw = TRUE)*ab_cat_shared + 
                    (1 | site),# + (1| grid:site), 
                  data = filter(temp_dat, species %in% c("brohor")),
                  ziformula = ~1,
                  family = poisson(link="log"))
mddb.3 <- glmmTMB(seed_half ~ ab_cat_shared*green_index_scaled + # doesn't run - infinite values
                    (1 | site),# + (1| grid:site), 
                  ziformula = ~.,
                  data = filter(temp_dat, species %in% c("brohor")),
                  family = poisson(link="log")) 
mddb.4 <- glmmTMB(seed_half ~ poly(green_index_scaled,2, raw = TRUE)*ab_cat_shared + 
                    (1 | site),# + (1| grid:site), 
                  ziformula = ~.,
                  data = filter(temp_dat, species %in% c("brohor")),
                  family = poisson(link="log"))

AIC(mddb.1,mddb.2,mddb.3,mddb.4)
# df      AIC
# mddb.1 10 554.5423 # linear w. zi as intercept
# mddb.2 14 560.4602
# mddb.4 26       NA

#############
# Check fit
#############

# using DHARMa:
testZeroInflation(mddb.1)
testDispersion(mddb.1) 
plot(fitted(mddb.1), residuals(mddb.1))
hist(residuals(mddb.1)) 
loc_sim_ouput <- simulateResiduals(mddb.1)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)
# all good

########
# Stats
########

summary(mddb.1)
Anova(mddb.1, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Chisq Df Pr(>Chisq)
# (Intercept)                      0.0010  1     0.9747
# ab_cat_shared                    0.9215  3     0.8202
# green_index_scaled               0.0012  1     0.9725
# ab_cat_shared:green_index_scaled 2.4594  3     0.4827

##############
# visualize
##############

vis_dd_b <-
  ggpredict(mddb.1,
            terms = c("green_index_scaled[all]","ab_cat_shared"),
            type = "fe.zi", allow.new.levels=TRUE) 
# fe.zi means seed production

plot(vis_dd_b)

ddb <- ggplot() +
  geom_line(data = vis_dd_b, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1
            # color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_dd_b,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high,
      fill = group
    ),
    alpha = 0.1,
    show.legend = F
    #fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 5) +
  ggtitle("Bromus") + 
  labs(x = "Vegetation index (G-R)", y= "", color = "Abundance") +
  annotate("text", x = -4, y = 5, label = "C", size = 6)
ddb

################################################################################
# Festuca


#########
# Models
#########
#rank deficient so need to remove interactions https://stackoverflow.com/questions/38766155/rank-deficiency-warning-mixed-model-lmer

mddv.1 <- glmmTMB(seed_half ~ ab_cat_shared + green_index_scaled + 
                    (1 | site) + (1| grid:site), 
                  data = filter(temp_dat, species %in% c("vulmic")),
                  ziformula = ~1,
                  family = poisson(link="log")) 
mddv.2 <- glmmTMB(seed_half ~ poly(green_index_scaled,2, raw = TRUE) + ab_cat_shared + 
                    (1 | site) + (1| grid:site), 
                  data = filter(temp_dat, species %in% c("vulmic")),
                  ziformula = ~1,
                  family = poisson(link="log"))
mddv.3 <- glmmTMB(seed_half ~ ab_cat_shared + green_index_scaled + 
                    (1 | site) + (1| grid:site), 
                  ziformula = ~.,
                  data = filter(temp_dat, species %in% c("vulmic")),
                  family = poisson(link="log")) 
mddv.4 <- glmmTMB(seed_half ~ poly(green_index_scaled,2, raw = TRUE) + ab_cat_shared + 
                    (1 | site) + (1| grid:site), 
                  ziformula = ~.,
                  data = filter(temp_dat, species %in% c("vulmic")),
                  family = poisson(link="log"))

AIC(mddv.1,mddv.2,mddv.3, mddv.4)
# df      AIC
# mddv.1  8 899.4781 # linear w. zi as intercept
# mddv.2  9 899.6512
# mddv.3 14 902.8260
# mddv.4 16 903.9124

#############
# Check fit
#############

# using DHARMa:
testZeroInflation(mddv.1)
testDispersion(mddv.1) 
plot(fitted(mddv.1), residuals(mddv.1))
hist(residuals(mddv.1)) 
loc_sim_ouput <- simulateResiduals(mddv.1)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)
# all good

########
# Stats
########

summary(mddv.1)
Anova(mddv.1, type = 3)

# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: seed_half
# Chisq Df Pr(>Chisq)    
#   (Intercept)        33.1628  1  8.476e-09 ***
#   ab_cat_shared       1.0836  3    0.78104    
#   green_index_scaled  3.8222  1    0.05058 . 

# no effect of abundance cat

##############
# visualize
##############

vis_dd_v <-
  ggpredict(mddv.1,
            terms = c("green_index_scaled[all]","ab_cat_shared"),
            type = "fe.zi", allow.new.levels=TRUE) 
# fe.zi means seed production

plot(vis_dd_v)

ddv <- ggplot() +
  geom_line(data = vis_dd_v, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1
            # color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_dd_v,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high,
      fill = group
    ),
    alpha = 0.1,
    show.legend = F
    #fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 20) +
  ggtitle("Festuca") + 
  labs(x = "Vegetation index (G-R)", y= "", color = "Abundance") +
  annotate("text", x = -4, y = 20, label = "D", size = 6)
ddv

#########################
# Greenness as category:
#########################

# note weirdness: festuca 1-10 conspecifics at prod = 2 is showing significant deviation, that is not shown in the continuous version. I think it is because we are coercing a polynomial fit into categories, which is wrong.

vis_dd_p <-
  ggpredict(mddp.3,
            terms = c("ab_cat_shared","green_index_scaled[-4,-1,2]"),
            type = "fe.zi", allow.new.levels=TRUE) 

vis_dd_m <-
  ggpredict(mddm.1,
            terms = c("ab_cat_shared","green_index_scaled[-4,-1,2]"),
            type = "fe.zi", allow.new.levels=TRUE) 

vis_dd_b <-
  ggpredict(mddb.1,
            terms = c("ab_cat_shared","green_index_scaled[-4,-1,2]"),
            type = "fe.zi", allow.new.levels=TRUE) 


vis_dd_v <-
  ggpredict(mddv.1,
            terms = c("ab_cat_shared","green_index_scaled[-4,-1,2]"),
            type = "fe.zi", allow.new.levels=TRUE) 

###############################
# PLOT ALL
###############################

p <- plot(vis_dd_p)
ddp <- p +
  theme_classic() +
  theme(legend.position = c(0.8,0.8),
        text = element_text(size = 16),
        axis.text.x  = element_blank()) + 
  labs(x = "Vegetative index (G-R)", y = "Seed production") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 5) +
  ggtitle("Plantago") + 
  labs(x = "", y= "", color = "Productivity") +  
  annotate("text", x = 0, y = 5, label = "B", size = 6)
ddp

p <- plot(vis_dd_m)
ddm <- p +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 16),
        axis.text.x  = element_blank()) +
  labs(x = "Vegetative index (G-R)", y = "Seed production") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 5) +
  ggtitle("Micropus") + 
  labs(x = "", y= "Seed Production", color = "Productivity") +  
  annotate("text", x = 0, y = 5, label = "A", size = 6)
ddm

p <- plot(vis_dd_b)
ddb <- p +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 16),
        axis.text.x  = element_text(angle = 45, hjust = 1)) + 
  labs(x = "Vegetative index (G-R)", y = "Seed production") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 5) +
  ggtitle("Bromus") + 
  labs(x = "", y= "Seed Production", color = "Productivity") +  
  annotate("text", x = 0, y = 5, label = "C", size = 6)
ddb


p <- plot(vis_dd_v)
ddv <- p +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 16),
        axis.text.x  = element_text(angle = 45, hjust = 1)) + 
  labs(x = "Vegetative index (G-R)", y = "Seed production") +
  geom_hline(yintercept = 1,
                                  linetype = "dashed",
                                  color = "black") +
  ylim(0, 5) +
  ggtitle("Festuca") + 
  labs(x = "", y= "", color = "Productivity") +  
  annotate("text", x = 0, y = 5, label = "D", size = 6)
ddv

png(
  'Figures/supp_fig_dd_cat.png',
  width = 8,
  height = 8,
  units = "in",
  res = 600
)
(ddm + ddp) / (ddb + ddv)
dev.off()

pdf(
  'Figures/supp_fig_dd_cat.pdf',
  width = 8,
  height = 8
)
(ddm + ddp) / (ddb + ddv)
dev.off()
