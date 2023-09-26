
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

source("Scripts/Source fitness data.R")

###############
## METHOD 2 ##
###############

# Looking at occ+per (dispersal categories) changing across abundance, with productivity accounted for.

##################
## Prepare data ##
##################

dd_temp1 <- green_all %>%
  dplyr::select(site, grid, plot, species, treatment, persistence, occurrence, green_index_scaled, contingency, replicates,ab_cat_shared) %>%
  distinct()
dd_temp1$contingency <-
  recode_factor(dd_temp1$contingency, "SS_n" = "SS")
dd_temp1$contingency <-
  recode_factor(dd_temp1$contingency, "SS_y" = "SS")
dd_temp1 <- dd_temp1 %>%
  filter(contingency %in% c('DL', "ME", 'SS')) %>%
   group_by(plot, site, grid, species, contingency, treatment) %>%
   mutate(cont_freq = n()) %>%
  dplyr::filter(., treatment %in% 'B') %>%
  ungroup() %>%
  dplyr::select(plot, green_index_scaled, contingency, cont_freq, site, grid, replicates, ab_cat_shared, species) %>%
  distinct() %>%
  pivot_wider(., names_from = 'contingency', values_from = 'cont_freq') %>%
  dplyr::select(-DL, -SS) # so that I don't have to use multinomial model, binomial data instead!
dd_temp1
dd_temp1[is.na(dd_temp1)] <- 0 # NA to 0

# check min and max of ME
min(dd_temp1$ME) # 0 good
max(dd_temp1$ME) # 1 good

# check variables
dd_temp1$ab_cat_shared <- as.factor(dd_temp1$ab_cat_shared)
dd_temp1$species <- as.factor(dd_temp1$species)

# Check: no sinks (persistence no, occurrence yes) should be in abundance cat 0
list_key2 <- dd_temp1 %>%
  filter(ab_cat_shared %in% '0' & ME > 0) %>%
  select(plot) %>%
  distinct()
list_key2
# good, it is zero

# change names for plotting
dd_temp1$ab_cat_shared <-
  recode_factor(dd_temp1$ab_cat_shared, 
                "0" = "0 conspecifics",
                "1" = "1-10 conspecifics",
                "2" = "11-100 conspecifics",
                "3" = "> 101 conspecifics")

##########
# MODELS
##########

################################################################################
# Bromus

mdd2.1b <- glmmTMB(ME ~ ab_cat_shared*green_index_scaled + 
                    (1 | site) + (1| grid:site), 
                  data = filter(dd_temp1, species %in% "brohor"),
                  family = binomial())
# mdd2.2b <- glmmTMB(ME ~ poly(green_index_scaled,2, raw = TRUE)*ab_cat_shared + 
#                     (1 | site) + (1| grid:site),
#                   data =  filter(dd_temp1, species %in% "brohor"),
#                   family = binomial())
# AIC(mdd2.1b)
# df      AIC
# mdd2.1b 10 153.3375
# mdd2.2b 14       NA

# Check rest of fit
# using DHARMa:
testZeroInflation(mdd2.1b)
testDispersion(mdd2.1b) 
plot(fitted(mdd2.1b), residuals(mdd2.1b))
hist(residuals(mdd2.1b)) 
loc_sim_ouput <- simulateResiduals(mdd2.1b)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)
# good

summary(mdd2.1b)
Anova(mdd2.1b, type = 3)

# Response: ME
# Chisq Df Pr(>Chisq)    
# (Intercept)                      0.0000  1     0.9991
# ab_cat_shared                    5.9495  3     0.1141
# green_index_scaled               0.0000  1     0.9999
# ab_cat_shared:green_index_scaled 2.1830  3     0.5353


################################################################################
# Micropus

# convergence - very small eigenvalue,
# soln: reduce random effect structure, remove interactions

mdd2.1m <- glmmTMB(ME ~ ab_cat_shared + green_index_scaled + 
                     (1 | site), #+ (1| grid:site), 
                   data = filter(dd_temp1, species %in% "miccal"),
                   family = binomial())
mdd2.2m <- glmmTMB(ME ~ poly(green_index_scaled,2, raw = TRUE) + ab_cat_shared + 
                     (1 | site), #+ (1| grid:site),
                   data =  filter(dd_temp1, species %in% "miccal"),
                   family = binomial())
AIC(mdd2.2m)
# df      AIC
# mdd2.1m  6 66.23596 false convergence
# mdd2.2m  7 68.06758 winner

# Check rest of fit
# using DHARMa:
testZeroInflation(mdd2.2m)
testDispersion(mdd2.2m) 
plot(fitted(mdd2.2m), residuals(mdd2.2m))
hist(residuals(mdd2.2m)) 
loc_sim_ouput <- simulateResiduals(mdd2.2m)
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

summary(mdd2.2m)
Anova(mdd2.2m)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: ME
# Chisq Df Pr(>Chisq)
# poly(green_index_scaled, 2, raw = TRUE) 0.7683  2      0.681
# ab_cat_shared                           0.0000  3      1.000


################################################################################
# Plantago

# mdd2.1p <- glmmTMB(ME ~ ab_cat_shared*green_index_scaled + 
#                      (1 | site), #+ (1| grid:site), 
#                    data = filter(dd_temp1, species %in% "plaere"),
#                    family = binomial())
mdd2.2p <- glmmTMB(ME ~ poly(green_index_scaled,2, raw = TRUE)*ab_cat_shared + 
                     (1 | site), #+ (1| grid:site),
                   data =  filter(dd_temp1, species %in% "plaere"),
                   family = binomial())
AIC(mdd2.2p)
# df      AIC
# dd2.1p 11 NA hessian error
# mdd2.2p 13 269.5955 winner

# Check rest of fit
# using DHARMa:
testZeroInflation(mdd2.2p)
testDispersion(mdd2.2p) 
plot(fitted(mdd2.2p), residuals(mdd2.2p))
hist(residuals(mdd2.2p)) 
loc_sim_ouput <- simulateResiduals(mdd2.2p)
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

summary(mdd2.2p)

Anova(mdd2.2p, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: ME
# Chisq Df Pr(>Chisq)
# (Intercept)                                           0.0000  1     0.9984
# poly(green_index_scaled, 2, raw = TRUE)               0.0000  2     1.0000
# ab_cat_shared                                         4.0439  3     0.2568
# poly(green_index_scaled, 2, raw = TRUE):ab_cat_shared 0.3283  6     0.9993


##############################################################################
# Festuca

#rank deficient: solution is to remove interaction, deviations in fit mean I have to reduce random effect structure as well. 

mdd2.1v <- glmmTMB(ME ~ ab_cat_shared + green_index_scaled + 
                     (1 | site),# + (1| grid:site), 
                   data = filter(dd_temp1, species %in% "vulmic"),
                   family = binomial())
mdd2.2v <- glmmTMB(ME ~ poly(green_index_scaled,2, raw = TRUE) + ab_cat_shared + 
                     (1 | site),# + (1| grid:site),
                   data =  filter(dd_temp1, species %in% "vulmic"),
                   family = binomial())
AIC(mdd2.1v,mdd2.2v)
# df      AIC
# mdd2.1v  6 313.1887
# mdd2.2v  7 311.4738 winner


# Check rest of fit
# using DHARMa:
testZeroInflation(mdd2.2v)
testDispersion(mdd2.2v) 
plot(fitted(mdd2.2v), residuals(mdd2.2v))
hist(residuals(mdd2.2v)) 
loc_sim_ouput <- simulateResiduals(mdd2.2v)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)
# bad deviations - try to reduce random effect structure, fixed it!

summary(mdd2.2v)

Anova(mdd2.2v, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: ME
# Chisq Df Pr(>Chisq)
# (Intercept)                             0.0000  1     0.9974
# poly(green_index_scaled, 2, raw = TRUE) 3.4516  2     0.1780
# ab_cat_shared                           0.1176  3     0.9896


###############################
# PLOT ALL
###############################

vis_dd_p <-
  ggpredict(mdd2.2p,
            terms = c("ab_cat_shared[1-10 conspecifics,11-100 conspecifics,> 101 conspecifics]", # remove 0 category
                      "green_index_scaled[-4, -1, 2]"),
            type = "fe", allow.new.levels=TRUE) 

vis_dd_m <-
  ggpredict(mdd2.2m,
            terms = c("ab_cat_shared[1-10 conspecifics,11-100 conspecifics,> 101 conspecifics]", # remove 0 category
                      "green_index_scaled[-4, -1, 2]"),
            type = "fe", allow.new.levels=TRUE) 

vis_dd_b <-
  ggpredict(mdd2.1b,
            terms = c("ab_cat_shared[1-10 conspecifics,11-100 conspecifics,> 101 conspecifics]", # remove 0 category
                      "green_index_scaled[-4, -1, 2]"),
            type = "fe", allow.new.levels=TRUE) 

vis_dd_v <-
  ggpredict(mdd2.2v,
            terms = c("ab_cat_shared[1-10 conspecifics,11-100 conspecifics,> 101 conspecifics]", # remove 0 category
                      "green_index_scaled[-4, -1, 2]"),
            type = "fe", allow.new.levels=TRUE) 

p <- plot(vis_dd_p)
sink_supp_p <- p +
  theme_classic() +
  theme(legend.position ='right',
        text = element_text(size = 16),
        axis.text.x  = element_blank()) + 
  ylim(0, 1) +
  ggtitle("Plantago") + 
  labs(x = "", y= "", color = "Productivity") +  
  annotate("text", x = 0, y = 1, label = "B", size = 6)
sink_supp_p

p <- plot(vis_dd_m)
sink_supp_m <- p +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 16),
        axis.text.x  = element_blank()) +
  ylim(0, 1) +
  ggtitle("Micropus") + 
  labs(x = "", y= "Proportion of sinks", color = "Productivity") +  
  annotate("text", x = 0, y = 1, label = "A", size = 6)
sink_supp_m

p <- plot(vis_dd_b)
sink_supp_b <- p +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 16),
        axis.text.x  = element_text(angle = 45, hjust = 1)) + 
  labs(x = "", y = "Proportion of sinks", color = "Productivity") +
  ylim(0, 1) +
  ggtitle("Bromus") + 
  annotate("text", x = 0, y = 1, label = "C", size = 6)
sink_supp_b


p <- plot(vis_dd_v)
sink_supp_v <- p +
  theme_classic() +
  theme(legend.position = 'none',
        text = element_text(size = 16),
        axis.text.x  = element_text(angle = 45, hjust = 1)) + 
  ylim(0, 1) +
  ggtitle("Festuca") + 
  labs(x = "", y= "", color = "Productivity") +  
  annotate("text", x = 0, y = 1, label = "D", size = 6)
sink_supp_v


png(
  'Figures/supp_fig_sinks_abundance.png',
  width = 9,
  height = 8,
  units = "in",
  res = 600
)
(sink_supp_m + sink_supp_p) / (sink_supp_b + sink_supp_v) + plot_layout(guides = 'collect') & theme(legend.position = 'right')
dev.off()

pdf(
  'Figures/supp_fig_sinks_abundance.pdf',
  width = 9,
  height = 8
)
(sink_supp_m + sink_supp_p) + (sink_supp_b + sink_supp_v)
dev.off()

###############################
# HOW COMMON ARE SINKS BY SPP
###############################

# instances of sinks
table <- dd_temp1 %>%
  group_by(species, ab_cat_shared) %>% 
  filter(ME == 1) %>%
  mutate(freq = n()) %>%
  select(species, ab_cat_shared, freq) %>%
  distinct() %>%
  arrange(ab_cat_shared, species)
view(table)

# total sample size for each species within abundances
table2 <- dd_temp1 %>%
  group_by(species, ab_cat_shared) %>% 
  mutate(freq = n()) %>%
  select(species, ab_cat_shared, freq) %>%
  distinct() %>%
  arrange(ab_cat_shared, species)
view(table2)










