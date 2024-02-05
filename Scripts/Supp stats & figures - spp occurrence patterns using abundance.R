#------------------------------------------------------------------------
# SPECIES SPECIFIC OCCURRENCE PATTERNS USING ABUNDANCE
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
# for each species (Plantago, Micropus, Bromus, and Vulpia (aka Festuca)) model occurrence patterns across the productivity gradient using numerical abundance instead of binary which is in main text figure/models

# NOTES:
# occurrence is a based on abundance category but treat as numerical counts, thus using poisson

# first change ab_cat_shared from 0-3 to 0->101
 green_all_supp <- green_all
 
 # calculate how often each species and block ends up in each abundance category
 
 green_all_supp %>%
   select(species, block, ab_cat_shared) %>%
   group_by(species) %>%
   summarize(abundance_total = n())
 # 1 brohor              708
 # 2 miccal              740
 # 3 plaere              724
 # 4 vulmic              759
 
 green_all_supp %>%
   select(species, block, ab_cat_shared) %>%
   group_by(species, ab_cat_shared) %>%
   summarize(abundance_name = n())
# brohor              0            428
# brohor              1            167
#  (428+167)/708 = 0.8403955
# miccal              0            269
# miccal              1            410
#  (269+410)/740 = 0.9175676
# plaere              0            345
# plaere              1            295
#  (345+295)/724 = 0.8839779
# vulmic              0            275
# vulmic              1            440
# (275+440)/759 = 0.942029

 green_all_supp %>%
   select(species, block, ab_cat_shared) %>%
   group_by(ab_cat_shared) %>%
   summarize(abundance_name = n()) 
# (1317+1312)/(1317+1312+256+46) = 0.8969635
 
 #--------------------------------------------------------------------
 # PLANTAGO ####
 #--------------------------------------------------------------------

 ######################
 #### Prepare data ####
 ######################
 
occ_plantago <- filter(green_all_supp, species %in% 'plaere' & treatment %in% 'B')
hist(occ_plantago$ab_cat_shared) 

#######################
####    Model      ####
#######################

#####linear  
# zi=1
mp1 <- glmmTMB(
  ab_cat_shared ~ green_index_scaled + (1 | site) + (1| grid:site),
  data = occ_plantago,
  family = poisson(link="log"),
  ziformula=~1
)

# zi=~.
mp2 <- glmmTMB(
  ab_cat_shared ~ green_index_scaled + (1 | site) + (1| grid:site),
  data = occ_plantago,
  family = poisson(link="log"),
  ziformula=~.
)
#####quadratic 
# zi =1
mp3 <- glmmTMB(
  ab_cat_shared ~ poly(green_index_scaled,2, raw=TRUE) +  (1 | site) + (1| grid:site),
  data = occ_plantago,
  family = poisson(link="log"),
  ziformula=~1
)
# zi =~.
mp4 <- glmmTMB(
  ab_cat_shared ~ poly(green_index_scaled,2, raw=TRUE) +  (1 | site) + (1| grid:site),
  data = occ_plantago,
  family = poisson(link="log"),
  ziformula=~.
)

BIC(mp1, mp2, mp3, mp4) # could be better to use if overfitting is possible
# df      BIC
# mp1  5 716.1003 win
# mp2  8       NA
# mp3  6 722.0389
# mp4 10       NA

AIC(mp1, mp2, mp3, mp4)
# df      AIC
# mp1  5 696.2696 win
# mp2  8       NA
# mp3  6 698.2420
# mp4 10       NA

mpo <- mp1
  
#####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(mpo)
testDispersion(mpo) 
plot(fitted(mpo), residuals(mpo))
hist(residuals(mpo)) 
loc_sim_ouput <- simulateResiduals(mpo); plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# KS dev but fine

#############################
## Interpretation of coefs ##
#############################

summary(mpo)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)         -0.7920     0.3054  -2.593  0.00951 **
#   green_index_scaled   0.2599     0.1183   2.197  0.02800 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -21.97    5444.01  -0.004    0.997
length(occ_plantago$grid) # 390

###############
## Visualize ##
###############

# for zi plots
vis_occ_p_zi <-
  ggpredict(mpo,
            terms = c("green_index_scaled[all]"),
            type = "zi_prob"); plot(vis_occ_p_zi)

# for conditional plot:
vis_occ_p <-
  ggpredict(mpo,
            terms = c("green_index_scaled[all]"),
            type = "fe.zi", allow.new.levels=TRUE); plot(vis_occ_p)

a <- ggplot() +
  geom_jitter(
    data = occ_plantago,
    aes(green_index_scaled, ab_cat_shared),
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
  labs(y = "Predicted abundance", x = "Vegetation index (G-R)") +
  ylim(0, 3) + 
  scale_y_continuous(breaks = c(0,1,2,3), labels =c("0","1-10","11-100",">101") ) + 
  ggtitle("Plantago") +
  annotate("text", x = -4, y = 3, label = "A", size = 7)
a

e <- ggplot() +
  geom_jitter(
    data = filter(occ_plantago, ab_cat_shared %in% 0),
    aes(green_index_scaled, ab_cat_shared),
    alpha = 0.1,
    size = 1.5,
    color = "midnightblue",
    height = 0.1
  ) +
  geom_line(data = vis_occ_p_zi, aes(x = x, y = predicted),
            color = "midnightblue",
            linewidth = 1) +
  geom_ribbon(
    data = vis_occ_p_zi,
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
  labs(y = "Predicted absent", x = "Vegetation index (G-R)") +
  ylim(0, 3) + 
  annotate("text", x = -4, y = 3, label = "E", size = 7)
e

#---------------------------------------
# MICROPUS ####
#---------------------------------------

##################
## Prepare data ##
##################

occ_micropus <- filter(green_all_supp, species %in% 'miccal' & treatment %in% 'B')
hist(occ_micropus$ab_cat_shared) 

###################
##    Model      ##
###################

#####linear  
# zi=1
mm1 <- glmmTMB(
  ab_cat_shared ~ green_index_scaled +  (1 | site) + (1| grid:site),
  data = occ_micropus,
  family = poisson(link="log"),
  ziformula=~1
)
# # zi=~.
# mm2 <- glmmTMB(
#   ab_cat_shared ~ green_index_scaled +  (1 | site) + (1| grid:site),
#   data = occ_micropus,
#   family = poisson(link="log"),
#   ziformula=~.
# )
#####quadratic 
# zi =1
mm3 <- glmmTMB(
  ab_cat_shared ~ poly(green_index_scaled,2, raw=TRUE) +  (1 | site) + (1| grid:site),
  data = occ_micropus,
  family = poisson(link="log"),
  ziformula=~1
)
# zi =~.
# mm4 <- glmmTMB(
#   ab_cat_shared ~ poly(green_index_scaled,2, raw=TRUE) +  (1 | site) + (1| grid:site),
#   data = occ_micropus,
#   family = poisson(link="log"),
#   ziformula=~.
# )

BIC(mm1,mm3)
# df      BIC
# mm1  5 798.0972 win
# mm3  6 802.1230
AIC(mm1,mm3)
# df      AIC
# mm1  5 778.1274 win
# mm3  6 778.1593

mmo <- mm1

#####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(mmo)
testDispersion(mmo) 
plot(fitted(mmo), residuals(mmo))
hist(residuals(mmo)) 
loc_sim_ouput <- simulateResiduals(mmo); plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# KS significant and dispersion

#############################
## Interpretation of coefs ##
#############################
summary(mmo)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)         -0.4925     0.1849  -2.664  0.00771 **
#   green_index_scaled   0.1636     0.1118   1.463  0.14336   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -21.95    4804.05  -0.005    0.996 
length(occ_micropus$grid) # 401

###############
## Visualize ##
###############

# for zi plots
vis_occ_m_zi <-
  ggpredict(mmo,
            terms = c("green_index_scaled[all]"),
            type = "zi_prob"); plot(vis_occ_m_zi)

# for conditional plot:
vis_occ_m <-
  ggpredict(mmo,
            terms = c("green_index_scaled[all]"),
            type = "fe.zi", allow.new.levels=TRUE); plot(vis_occ_m)

b <- ggplot() +
  geom_jitter(
    data = occ_micropus,
    aes(green_index_scaled, ab_cat_shared),
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
  labs(y = "", x = "Vegetation index (G-R)") +
  ylim(0, 3) + 
  scale_y_continuous(breaks = c(0,1,2,3), labels =c("0","1-10","11-100",">101") ) + 
  ggtitle("Micropus") +
  annotate("text", x = -4, y = 3, label = "B", size = 7)
b

f <- ggplot() +
  geom_jitter(
    data = filter(occ_micropus, ab_cat_shared %in% 0),
    aes(green_index_scaled, ab_cat_shared),
    alpha = 0.1,
    size = 1.5,
    color = "midnightblue",
    height = 0.1
  ) +
  geom_line(data = vis_occ_m_zi, aes(x = x, y = predicted),
            color = "midnightblue",
            linewidth = 1) +
  geom_ribbon(
    data = vis_occ_m_zi,
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
  ylim(0, 3) + 
  annotate("text", x = -4, y = 3, label = "F", size = 7)
f

#--------------------------------------------------------------------
# BROMUS
#--------------------------------------------------------------------

##################
## Prepare data ##
##################

occ_bromus <- filter(green_all_supp, species %in% 'brohor' & treatment %in% 'B')
hist(occ_bromus$ab_cat_shared) 

############
## Model  ##
############

#####linear 
# zi=1
mb1 <- glmmTMB(
  ab_cat_shared ~ green_index_scaled +  (1 | site) + (1| grid:site),
  data = occ_bromus,
  family = poisson(link='log'),
  ziformula=~1
) 
# zi=~.
# mb2 <- glmmTMB(
#   ab_cat_shared ~ green_index_scaled +  (1 | site) + (1| grid:site),
#   data = occ_bromus,
#   family = poisson(link='log'),
#   ziformula=~.
# )
#####quadratic 
# zi =1
mb3 <- glmmTMB(
  ab_cat_shared ~ poly(green_index_scaled,2, raw=TRUE) +  (1 | site) + (1| grid:site),
  data = occ_bromus,
  family = poisson(link='log'),
  ziformula=~1
) 
# zi =~.
mb4 <- glmmTMB(
  ab_cat_shared ~ poly(green_index_scaled,2, raw=TRUE) +  (1 | site) + (1| grid:site),
  data = occ_bromus,
  family = poisson(link='log'),
  ziformula=~.
)

BIC(mb1,mb3,mb4)
# df      BIC
# mb1  5 664.2377 win
# mb3  6 666.5755
# mb4 10       NA
AIC(mb1,mb3,mb4)
# df      AIC
# mb1  5 644.5898
# mb3  6 642.9980 win
# mb4 10       NA

# disagree, use BIC for simplicity

mbo <- mb1
  
#####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(mbo)
testDispersion(mbo) 
plot(fitted(mbo), residuals(mbo))
hist(residuals(mbo)) 
loc_sim_ouput <- simulateResiduals(mbo)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# KS deviation but fine

#############################
## Interpretation of coefs ##
#############################

summary(mbo)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)         -1.2882     0.5350  -2.408  0.01604 * 
#   green_index_scaled   0.4324     0.1324   3.266  0.00109 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -21.44    5026.92  -0.004    0.997
length(occ_bromus$grid) # 376

###############
## Visualize ##
###############

# for zi plots
vis_occ_b_zi <-
  ggpredict(mbo,
            terms = c("green_index_scaled[all]"),
            type = "zi_prob",allow.new.levels=TRUE); plot(vis_occ_b_zi)

# for conditional plot:
vis_occ_b <-
  ggpredict(mbo,
            terms = c("green_index_scaled[all]"),
            type = "fe",allow.new.levels=TRUE); plot(vis_occ_b)

c <- ggplot() +
  geom_jitter(
    data = occ_bromus,
    aes(green_index_scaled, ab_cat_shared),
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
  ylim(0, 3) + 
  scale_y_continuous(breaks = c(0,1,2,3), labels =c("0","1-10","11-100",">101") ) + 
  ggtitle("Bromus") +
  annotate("text", x = -4, y = 3, label = "C", size = 7)
c

j <- ggplot() +
  geom_jitter(
    data = filter(occ_bromus, ab_cat_shared %in% 0),
    aes(green_index_scaled, ab_cat_shared),
    alpha = 0.1,
    size = 1.5,
    color = "midnightblue",
    height = 0.1
  ) +
  geom_line(data = vis_occ_b_zi, aes(x = x, y = predicted),
            color = "midnightblue",
            linewidth = 1) +
  geom_ribbon(
    data = vis_occ_b_zi,
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
  ylim(0, 3) + 
  annotate("text", x = -4, y = 3, label = "J", size = 7)
j

#--------------------------------------------------------------------
# VULPIA (now known by FESTUCA)
#--------------------------------------------------------------------

##################
## Prepare data ##
##################

occ_vulpia <- filter(green_all_supp, species %in% 'vulmic' & treatment %in% 'B')
hist(occ_vulpia$ab_cat_shared) 

###########
## Model ##
###########

#####linear 
# zi=1
mv1 <- glmmTMB(
  ab_cat_shared ~ green_index_scaled +  (1 | site) + (1| grid:site),
  data = occ_vulpia,
  family = poisson(link="log"),
  ziformula=~1
)
# zi=~.
mv2 <- glmmTMB(
  ab_cat_shared ~ green_index_scaled +  (1 | site) + (1| grid:site),
  data = occ_vulpia,
  family = poisson(link="log"),
  ziformula=~.
) 

#####quadratic 
# zi =1
mv3 <- glmmTMB(
  ab_cat_shared ~ poly(green_index_scaled,2, raw=TRUE) +  (1 | site) + (1| grid:site),
  data = occ_vulpia,
  family = poisson(link="log"),
  ziformula=~1
)
# zi =~.
mv4 <- glmmTMB(
  ab_cat_shared ~ poly(green_index_scaled,2, raw=TRUE) +  (1 | site) + (1| grid:site),
  data = occ_vulpia,
  family = poisson(link="log"),
  ziformula=~.
)

BIC(mv1,mv2,mv3,mv4)
# df      BIC
# mv1  5 786.3237 win
# mv2  8       NA
# mv3  6 789.3174
# mv4 10       NA
AIC(mv1,mv2,mv3,mv4)
# df      AIC
# mv1  5 766.4293
# mv2  8       NA
# mv3  6 765.4441 win
# mv4 10       NA

# disagree, go with BIC for simplicity
mvo <- mv1
  
####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(mvo)
testDispersion(mvo) 
plot(fitted(mvo), residuals(mvo))
hist(residuals(mvo)) 
loc_sim_ouput <- simulateResiduals(mvo)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# KS and dispersion issues

#############################
## Interpretation of coefs ##
#############################

summary(mvo)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)         -0.4752     0.1601  -2.968   0.0030 **
#   green_index_scaled   0.2576     0.1168   2.206   0.0274 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -22.07    5321.65  -0.004    0.997
length(occ_vulpia$grid) #395

###############
## Visualize ##
###############

# for zi plots
vis_occ_v_zi <-
  ggpredict(mvo,
            terms = c("green_index_scaled[all]"),
            type = "zi_prob",allow.new.levels=TRUE); plot(vis_occ_v_zi)

# for conditional plot:
vis_occ_v <-
  ggpredict(mvo,
            terms = c("green_index_scaled[all]"),
            type = "fe",allow.new.levels=TRUE); plot(vis_occ_v)

d <- ggplot() +
  geom_jitter(
    data = occ_vulpia,
    aes(green_index_scaled, ab_cat_shared),
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
  ylim(0, 3) + 
  scale_y_continuous(breaks = c(0,1,2,3), labels =c("0","1-10","11-100",">101") ) + 
  ggtitle("Festuca") +
  annotate("text", x = -4, y = 3, label = "D", size = 7)
d

h <- ggplot() +
  geom_jitter(
    data = filter(occ_vulpia, ab_cat_shared %in% 0),
    aes(green_index_scaled, ab_cat_shared),
    alpha = 0.1,
    size = 1.5,
    color = "midnightblue",
    height = 0.1
  ) +
  geom_line(data = vis_occ_v_zi, aes(x = x, y = predicted),
            color = "midnightblue",
            linewidth = 1) +
  geom_ribbon(
    data = vis_occ_v_zi,
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
  ylim(0, 3) + 
  annotate("text", x = -4, y = 3, label = "H", size = 7)
h

############################
# PULL TOGETHER AND SAVE
############################

pdf("Figures/supp_occ_w_abundance.pdf",
    width = 13,
    height = 5)
(a + b + c + d + plot_layout(ncol = 4))
dev.off()


png("Figures/supp_occ_w_abundance.png",
    width = 13,
    height = 5,
    units = "in",
    res = 600)
(a + b + c + d + plot_layout(ncol = 4))
dev.off()




