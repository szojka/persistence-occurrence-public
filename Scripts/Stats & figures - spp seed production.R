
#------------------------------------------------------------.
# SPECIES-SPECIFIC SEED PRODUCTION ACROSS GRADIENT ####
#------------------------------------------------------------.

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

#transformations used below are coded in Source
source("Scripts/Source fitness data.R")
source("Scripts/Stats & figures - spp occurrence patterns.R")

# Description of script: 
# for each species (Plantago, Micropus, Bromus, and Vulpia (aka Festuca)), we... 
#1. test linear vs quadratic models with either zero-inflation as an intercept or allowed to vary across the gradient. 
  # biologically, these different forms of zero inflation indicate if germination success is a constant across the gradient, or is better fit as a function across the gradient.
#2. We then use BIC and AIC to find the best model. We use the DHARMa package to ensure this selected model meets assumptions. 
#3. We then interpret the model using summary() and Anova() from the car package. 
#4. Finally we visualize the model outputs over the raw data points. 

#--------------------------------------------------------------------
# PLANTAGO ####
#--------------------------------------------------------------------

######################
#### Prepare data ####
######################

seed_plantago <- filter(green_all, species %in% "plaere")
# # find total sample size
# N = 724
# find non-zero sample size
seed_plantago %>%
  filter(seed>0) %>%
  mutate(frequency = n()) %>%
  select(frequency) %>%
  distinct()
# N = 293

ggplot(seed_plantago, aes(x = seed_half)) +
  geom_histogram(aes(fill = treatment),
                 position = "identity", alpha = 0.5) # ziformula=~1

#######################
####    Model      ####
#######################

#####linear  
# zi=1
mp1 <- glmmTMB(
  seed_half ~ green_index_scaled * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_plantago,
  family = poisson(link="log"),
  ziformula=~1
)

# zi=~.
mp2 <- glmmTMB(
  seed_half ~ green_index_scaled * treatment + (1 | block) + (1 | site) + (1| grid:site),
  data = seed_plantago,
  family = poisson(link="log"),
  ziformula=~.
)
#####quadratic 
# zi =1
mp3 <- glmmTMB(
  seed_half ~ poly(green_index_scaled,2, raw=TRUE) * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_plantago,
  family = poisson(link="log"),
  ziformula=~1
)
# zi =~.
mp4 <- glmmTMB(
  seed_half ~ poly(green_index_scaled,2, raw=TRUE) * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_plantago,
  family = poisson(link="log"),
  ziformula=~.
)

BIC(mp1, mp2, mp3, mp4) # could be better to use if overfitting is possible
# df      BIC
# mp1  8 2177.637 Winner
# mp2 14 2191.503
# mp3 10 2187.842
# mp4 18 2205.665

AIC(mp1, mp2, mp3, mp4)
# df      AIC
# mp1  8 2140.959
# mp2 14 2127.316
# mp3 10 2141.994
# mp4 18 2123.139 winner

# AIC and BIC disagree
# for now trust AIC, and BIC works better with small N. My N is 724
# going ahead with mp4

mp <- mp4

#####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(mp)
testDispersion(mp) 
plot(fitted(mp), residuals(mp))
hist(residuals(mp)) 
loc_sim_ouput <- simulateResiduals(mp)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# no issues except resisutatl vs fitted quantile plots. Not a deal breaker.

#############################
## Interpretation of coefs ##
#############################

summary(mp)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                                          1.13768    0.37413   3.041  0.00236 ** 
#   poly(green_index_scaled, 2, raw = TRUE)1             0.09236    0.11727   0.788  0.43093    
#   poly(green_index_scaled, 2, raw = TRUE)2            -0.09145    0.05536  -1.652  0.09854 .  
#   treatmentB                                          -1.39503    0.10788 -12.931  < 2e-16 ***
#   poly(green_index_scaled, 2, raw = TRUE)1:treatmentB -0.10149    0.12059  -0.842  0.40002    
#   poly(green_index_scaled, 2, raw = TRUE)2:treatmentB  0.10389    0.06213   1.672  0.09453 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                                          0.08428    0.28405   0.297 0.766693    
#   poly(green_index_scaled, 2, raw = TRUE)1            -0.57001    0.22574  -2.525 0.011569 *  
#   poly(green_index_scaled, 2, raw = TRUE)2            -0.35122    0.14634  -2.400 0.016394 *  
#   treatmentB                                           0.05096    0.26950   0.189 0.850013    
#   poly(green_index_scaled, 2, raw = TRUE)1:treatmentB  1.16840    0.30938   3.777 0.000159 ***
#   poly(green_index_scaled, 2, raw = TRUE)2:treatmentB -0.10912    0.26272  -0.415 0.677891    

# interested in the main effects vs zero inflation components:

Anova(mp,type=3,component="cond")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: seed_half
# Chisq Df Pr(>Chisq)    
#   (Intercept)                                         9.2469  1   0.002359 ** 
#   poly(green_index_scaled, 2, raw = TRUE)             3.5154  2   0.172442    
#   treatment                                         167.2165  1  < 2.2e-16 ***
#   poly(green_index_scaled, 2, raw = TRUE):treatment   4.9614  2   0.083683 .  

Anova(mp,type=3,component="zi")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: seed_half
# Chisq Df Pr(>Chisq)    
#   (Intercept)                                        0.0880  1  0.7666934    
#   poly(green_index_scaled, 2, raw = TRUE)           10.8950  2  0.0043072 ** 
#   treatment                                          0.0358  1  0.8500129    
#   poly(green_index_scaled, 2, raw = TRUE):treatment 14.6546  2  0.0006574 ***

###############
## emtrends ##
###############

#set max.degree = 2 if quadratic in the model
emtrends(mp, ~ treatment | green_index_scaled, "green_index_scaled", component="cond", max.degree = 2)

# green_index_scaled = -0.00483, degree = linear:
#   treatment green_index_scaled.trend     SE  df lower.CL upper.CL
# A                          0.09379 0.1172 706   -0.136   0.3239
# B                         -0.00932 0.1381 706   -0.281   0.2619
# 
# green_index_scaled = -0.00483, degree = quadratic:
#   treatment green_index_scaled.trend     SE  df lower.CL upper.CL
# A                         -0.09145 0.0554 706   -0.200   0.0172
# B                          0.01243 0.0613 706   -0.108   0.1327
# 
# Confidence level used: 0.95 


emtrends(mp, ~ treatment | green_index_scaled, "green_index_scaled", component="zi",max.degree = 2)
# 
# 
# green_index_scaled = -0.00483, degree = linear:
#   treatment green_index_scaled.trend    SE  df lower.CL upper.CL
# A                           -0.565 0.225 706  -1.0072  -0.1218
# B                            0.606 0.296 706   0.0241   1.1871
# 
# green_index_scaled = -0.00483, degree = quadratic:
#   treatment green_index_scaled.trend    SE  df lower.CL upper.CL
# A                           -0.351 0.146 706  -0.6385  -0.0639
# B                           -0.460 0.249 706  -0.9495   0.0288
# 
# Confidence level used: 0.95 

###############
## Visualize ##
###############

# these ggeffects all succesfully calculated CIs. 

vis_p_all <-
  ggpredict(mp,
            terms = c("green_index_scaled[all]","treatment"),
            type = "zi_prob"); plot(vis_p_all)

# for final plot:
vis_spp_p <-
  ggpredict(mp,
            terms = c("green_index_scaled[all]", "treatment"),
            type = "fe.zi", allow.new.levels=TRUE); plot(vis_spp_p)

# add color as column for plotting
seed_plantago$colortreat <- NA
seed_plantago$colortreat[seed_plantago$treatment == "A"]<- "steelblue3"
seed_plantago$colortreat[seed_plantago$treatment == "B"]<- "red4"
seed_plantago$colortreat <- as.factor(seed_plantago$colortreat)
colortreatraw <- c(seed_plantago$colortreat)
length(vis_spp_p$x) # 710/2 = 355
colortreatpred <- c(rep("steelblue3", times = 355), rep("red4", times = 355))

w <- ggplot() +
  geom_jitter(
    data = seed_plantago,
    aes(green_index_scaled, seed_half, color = treatment),
    alpha = 0.1,
    size = 1.5,
    color = colortreatraw
  ) +
  geom_line(data = vis_spp_p, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            color =colortreatpred
  ) +
  geom_ribbon(
    data = vis_spp_p,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
      ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
    ) +
  theme_classic() + # should always be above other themes
  theme(legend.title = element_blank(),
        legend.position = 'none',
        text = element_text(size = 16)) + 
  labs(y = "", x = "") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 15) +
  ggtitle("Plantago")+ 
  annotate("text", x = -4, y = 15, label = "B", size = 7)
w

#--------------------------------------------------------------------
# MICROPUS 
#--------------------------------------------------------------------

##################
## Prepare data ##
##################

seed_micropus <- filter(green_all, species %in% "miccal")
# N = 740 
seed_micropus %>%
  filter(seed>0) %>%
  mutate(frequency = n()) %>%
  select(frequency) %>%
  distinct()
# N = 34
ggplot(seed_micropus, aes(x = seed)) +
  geom_histogram(aes(fill = treatment),
                 position = "identity", alpha = 0.5) # ziformula=~treatment bc abiotic fit better without zi

###################
##    Model      ##
###################

#####linear  
# zi=1
mm1 <- glmmTMB(
  seed_half ~ green_index_scaled * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_micropus,
  family = poisson(link="log"),
  ziformula=~1
)
# zi=~.
mm2 <- glmmTMB(
  seed_half ~ green_index_scaled * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_micropus,
  family = poisson(link="log"),
  ziformula=~.
)
#####quadratic 
# zi =1
mm3 <- glmmTMB(
  seed_half ~ poly(green_index_scaled,2, raw=TRUE) * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_micropus,
  family = poisson(link="log"),
  ziformula=~1
)
# zi =~.
mm4 <- glmmTMB(
  seed_half ~ poly(green_index_scaled,2, raw=TRUE) * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_micropus,
  family = poisson(link="log"),
  ziformula=~.
)

BIC(mm1,mm2,mm3,mm4)
# df      BIC
# mm1  8 353.4620 winner
# mm2 14 388.1800
# mm3 10 361.9736
# mm4 18 402.6239
AIC(mm1,mm2,mm3,mm4)
# df      AIC
# mm1  8 316.6088
# mm2 14 323.6869
# mm3 10 315.9071 winner
# mm4 18 319.7042

# They don't agree, and BIC is better for small sample so going with that

mm <- mm1

#####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(mm)
testDispersion(mm) 
plot(fitted(mm), residuals(mm))
hist(residuals(mm)) 
loc_sim_ouput <- simulateResiduals(mm)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# no probs here

###############
## emtrends ##
###############

#set max.degree = 2 if quadratic in the model
emtrends(mm, ~ treatment | green_index_scaled, "green_index_scaled", component="cond",max.degree = 1)
# not undefined columns are selected...? Emtrends was mistaking 'T' for a variable when I refered to it as TRUE. Simple mistake.

# green_index_scaled = 0.0365:
#   treatment green_index_scaled.trend    SE  df lower.CL upper.CL
# A                            0.512 0.611 732  -0.6864     1.71
# B                            0.998 0.519 732  -0.0203     2.02
# 
# Confidence level used: 0.95 

# zi = 1, not relevent:
#emtrends(mm, ~ treatment | green_index_scaled, "green_index_scaled", component="zi",max.degree = 1)


#############################
## Interpretation of coefs ##
#############################

Anova(mm,type=3,component="cond")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: seed_half
# Chisq Df Pr(>Chisq)  
# (Intercept)                  3.4564  1    0.06301 .
# green_index_scaled           0.7038  1    0.40151  
# treatment                    6.3185  1    0.01195 *
# green_index_scaled:treatment 0.4219  1    0.51601 

Anova(mm,type=3,component="zi")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: seed_half
# Chisq Df Pr(>Chisq)    
# (Intercept) 115.07  1  < 2.2e-16 ***

summary(mm)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                     0.8322     0.4476   1.859   0.0630 .
# green_index_scaled              0.5122     0.6106   0.839   0.4015  
# treatmentB                     -1.6065     0.6391  -2.514   0.0119 *
# green_index_scaled:treatmentB   0.4858     0.7480   0.649   0.5160  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    2.896      0.270   10.73   <2e-16 ***

###############
## Visualize ##
###############

# confidence intervals all work

vis_m_all <-
  ggpredict(mm,
            terms = c("green_index_scaled[all]","treatment"),
            type = "zi_prob"); plot(vis_m_all)

# for final plot:
vis_spp_m <-
  ggpredict(mm,
            terms = c("green_index_scaled[all]", "treatment"),
            type = "fe.zi", allow.new.levels=TRUE); plot(vis_spp_m)

# add color as column for plotting
seed_micropus$colortreat <- NA
seed_micropus$colortreat[seed_micropus$treatment == "A"]<- "steelblue3"
  seed_micropus$colortreat[seed_micropus$treatment == "B"]<- "red4"
    seed_micropus$colortreat <- as.factor(seed_micropus$colortreat)
    colortreatraw <- c(seed_micropus$colortreat)
   length(vis_spp_m$x)# 744/2
    colortreatpred <- c(rep("steelblue3", times = 372), rep("red4", times = 372))

x <- ggplot() +
      geom_jitter(
        data = seed_micropus,
        aes(green_index_scaled, seed_half, color = treatment),
        alpha = 0.1,
        size = 1.5,
        show.legend = F,
        color = colortreatraw
      ) +
      geom_line(data = vis_spp_m, aes(x = x, y = predicted, group = group, color = group),
                linewidth = 1,
                #color = colortreatpred,
                show.legend = T
      ) +
  scale_color_manual(values =c("steelblue3","red4"), 
                     labels = c("Without neighbors","With neighbors")) +
      geom_ribbon(
        data = vis_spp_m,
        aes(
          x = x,
          y = predicted,
          group=group,
          ymin = conf.low,
          ymax = conf.high
        ),
        alpha = 0.1,
        show.legend = F,
        fill = colortreatpred
      ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.49,0.7),
        text = element_text(size = 16)) + 
  labs(y = "Population growth", x = "") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 15) +
  ggtitle("Micropus") + 
  annotate("text", x = -4, y = 15, label = "A", size = 7)
x 

#--------------------------------------------------------------------
# BROMUS
#--------------------------------------------------------------------

##################
## Prepare data ##
##################

seed_bromus <- filter(green_all, species %in% "brohor")
# N = 708
seed_bromus %>%
  filter(seed>0) %>%
  mutate(frequency = n()) %>%
  select(frequency) %>%
  distinct()
# N = 272
# look at this for zi formula
ggplot(seed_bromus, aes(x = seed_half)) +
  geom_histogram(aes(fill = treatment),
                 position = "identity", alpha = 0.5) # ziformula=~1

#####################################
## Model selection seed production ##
#####################################

#####linear 
# zi=1
mb1 <- glmmTMB(
  seed_half ~ green_index_scaled * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_bromus,
  family = nbinom1(link='log'),
  ziformula=~1
) # CONVERGENCE ISSUE non-postitve-definite hessian matrix 
# zi=~.
mb2 <- glmmTMB(
  seed_half ~ green_index_scaled * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_bromus,
  family = nbinom1(link='log'),
  ziformula=~.
)
#####quadratic 
# zi =1
mb3 <- glmmTMB(
  seed_half ~ poly(green_index_scaled,2, raw=TRUE) * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_bromus,
  family = nbinom1(link='log'),
  ziformula=~1
) # CONVERGENCE ISSUE nonpostitive dfinite hessian matrix AND false convergence
# zi =~.
mb4 <- glmmTMB(
  seed_half ~ poly(green_index_scaled,2, raw=TRUE) * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_bromus,
  family = nbinom1(link='log'),
  ziformula=~.
)

BIC(mb1,mb2,mb3,mb4)
# df      BIC
# mb1  9       NA
# mb2 15 1581.987 winner
# mb3 11       NA
# mb4 19 1607.157
AIC(mb1,mb2,mb3,mb4)
# df      AIC
# mb1  9       NA
# mb2 15 1513.550 winner
# mb3 11       NA
# mb4 19 1520.471

mb <- mb2

#####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(mb)
testDispersion(mb) 
plot(fitted(mb), residuals(mb))
hist(residuals(mb)) 
loc_sim_ouput <- simulateResiduals(mb)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# no worries.

###############
## emtrends ##
###############

#set max.degree = 2 if quadratic in the model
emtrends(mb, ~ treatment | green_index_scaled, "green_index_scaled", component="cond",max.degree = 1)
# green_index_scaled = 0.0165:
#   treatment green_index_scaled.trend    SE  df lower.CL upper.CL
# A                            0.638 0.151 693    0.342   0.9338
# B                           -0.366 0.169 693   -0.697  -0.0343
# 
# Confidence level used: 0.95 

emtrends(mb, ~ treatment | green_index_scaled, "green_index_scaled", component="zi",max.degree = 1)
# 
# green_index_scaled = 0.0165:
#   treatment green_index_scaled.trend    SE  df lower.CL upper.CL
# A                            0.575 0.462 693   -0.333    1.482
# B                           -0.579 0.757 693   -2.067    0.908
# 
# Confidence level used: 0.95 


#############################
## Interpretation of coefs ##
#############################

Anova(mb,type=3,component="cond")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: seed_half
# Chisq Df Pr(>Chisq)    
#   (Intercept)                   6.8246  1   0.008991 ** 
#   green_index_scaled           17.9391  1  2.281e-05 ***
#   treatment                    91.0498  1  < 2.2e-16 ***
#   green_index_scaled:treatment 30.9259  1  2.681e-08 ***

Anova(mb,type=3,component="zi")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: seed_half
# Chisq Df Pr(>Chisq)  
# (Intercept)                  1.2262  1     0.2682
# green_index_scaled           1.5456  1     0.2138
# treatment                    2.0773  1     0.1495
# green_index_scaled:treatment 1.7560  1     0.1851

summary(mb)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                     0.4781     0.1830   2.612  0.00899 ** 
#   green_index_scaled              0.6380     0.1506   4.235 2.28e-05 ***
#   treatmentB                     -1.8398     0.1928  -9.542  < 2e-16 ***
#   green_index_scaled:treatmentB  -1.0038     0.1805  -5.561 2.68e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                    -0.5633     0.5087  -1.107    0.268
# green_index_scaled              0.5748     0.4623   1.243    0.214
# treatmentB                     -2.6494     1.8383  -1.441    0.150
# green_index_scaled:treatmentB  -1.1541     0.8709  -1.325    0.185

###############
## Visualize ##
###############

vis_b_all <-
  ggpredict(mb,
            terms = c("green_index_scaled[all]","treatment"),
            type = "zi_prob"); plot(vis_b_all)

# for final plot:
vis_spp_b <-
  ggpredict(mb,
            terms = c("green_index_scaled[all]", "treatment"),
            type = "fe.zi", allow.new.levels=TRUE); plot(vis_spp_b)


# add color as column for plotting
seed_bromus$colortreat <- NA
seed_bromus$colortreat[seed_bromus$treatment == "A"]<- "steelblue3"
  seed_bromus$colortreat[seed_bromus$treatment == "B"]<- "red4"
    seed_bromus$colortreat <- as.factor(seed_bromus$colortreat)
    colortreatraw <- c(seed_bromus$colortreat)
    length(vis_spp_b$x)# 700/2 # length of vis_spp_b
    colortreatpred <- c(rep("steelblue3", times = 351), rep("red4", times = 351))
    
y <- ggplot() +
      geom_jitter(
        data = seed_bromus,
        aes(green_index_scaled, seed_half, color = treatment),
        alpha = 0.1,
        size = 1.5,
        color = colortreatraw
      ) +
      geom_line(data = vis_spp_b, aes(x = x, y = predicted, group = group, color = group),
                linewidth = 1,
                color = colortreatpred
      ) +
      geom_ribbon(
        data = vis_spp_b,
        aes(
          x = x,
          y = predicted,
          group=group,
          ymin = conf.low,
          ymax = conf.high
        ),
        alpha = 0.1,
        show.legend = F,
        fill = colortreatpred
      ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 15) +
  ggtitle("Bromus")+ 
  annotate("text", x = -4, y = 15, label = "C", size = 7)
y

#--------------------------------------------------------------------
# VULPIA
#--------------------------------------------------------------------

##################
## Prepare data ##
##################

seed_vulpia <- filter(green_all, species %in% "vulmic")
# N = 759
seed_vulpia %>%
  filter(seed>0) %>%
  mutate(frequency = n()) %>%
  select(frequency) %>%
  distinct()
# N = 308
# look at this for zi formula
ggplot(seed_vulpia, aes(x = seed_half)) +
  geom_histogram(aes(fill = treatment),
                 position = "identity", alpha = 0.5) # ziformula=~1

#####################################
## Model selection seed production ##
#####################################

#####linear 
# zi=1
mv1 <- glmmTMB(
  seed_half ~ green_index_scaled * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_vulpia,
  family = poisson(link="log"),
  ziformula=~1
)
# zi=~.
mv2 <- glmmTMB(
  seed_half ~ green_index_scaled * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_vulpia,
  family = poisson(link="log"),
  ziformula=~.
) # CONVERGENCE ISSUE non-positive hessian matrix
#####quadratic 
# zi =1
mv3 <- glmmTMB(
  seed_half ~ poly(green_index_scaled,2, raw=TRUE) * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_vulpia,
  family = poisson(link="log"),
  ziformula=~1
)
# zi =~.
mv4 <- glmmTMB(
  seed_half ~ poly(green_index_scaled,2, raw=TRUE) * treatment + (1 | block) +  (1 | site) + (1| grid:site),
  data = seed_vulpia,
  family = poisson(link="log"),
  ziformula=~.
)

BIC(mv1,mv2,mv3,mv4)
# df      BIC
# mv1  8 2382.351 winner
# mv2 14       NA
# mv3 10 2390.531
# mv4 18 2407.640
AIC(mv1,mv2,mv3,mv4)
# df      AIC
# mv1  8 2345.295
# mv2 14       NA
# mv3 10 2344.211
# mv4 18 2324.264 winner

# AIC and BIC disagree, Vulpia doens't suffer from overfitting, so going with AIC

mv <- mv4
  
#####################
## Check model fit ##
#####################

# using DHARMa:
testZeroInflation(mv)
testDispersion(mv) 
plot(fitted(mv), residuals(mv))
hist(residuals(mv)) 
loc_sim_ouput <- simulateResiduals(mv)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

# no problems at all

###############
## emtrends ##
###############

#set max.degree = 2 if quadratic in the model
emtrends(mv, ~ treatment | green_index_scaled, "green_index_scaled", component="cond",max.degree = 2)

# green_index_scaled = -0.0464, degree = linear:
#   treatment green_index_scaled.trend     SE  df lower.CL upper.CL
# A                           0.0926 0.1014 741  -0.1064  0.29163
# B                           0.2571 0.1220 741   0.0176  0.49664
# 
# green_index_scaled = -0.0464, degree = quadratic:
#   treatment green_index_scaled.trend     SE  df lower.CL upper.CL
# A                          -0.1169 0.0560 741  -0.2268 -0.00694
# B                          -0.0914 0.0689 741  -0.2267  0.04381
# 
# Confidence level used: 0.95 


emtrends(mv, ~ treatment | green_index_scaled, "green_index_scaled", component="zi",max.degree = 2)

# green_index_scaled = -0.0464, degree = linear:
#   treatment green_index_scaled.trend     SE  df lower.CL upper.CL
# A                          -0.4653 0.1645 741   -0.788   -0.142
# B                           0.1462 0.1822 741   -0.211    0.504
# 
# green_index_scaled = -0.0464, degree = quadratic:
#   treatment green_index_scaled.trend     SE  df lower.CL upper.CL
# A                           0.0085 0.0962 741   -0.180    0.197
# B                          -0.1364 0.1287 741   -0.389    0.116
# 
# Confidence level used: 0.95 

#############################
## Interpretation of coefs ##
#############################

Anova(mv,type=3,component="cond")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: seed_half
# Chisq Df Pr(>Chisq)    
#   (Intercept)                                       198.0443  1    < 2e-16 ***
#   poly(green_index_scaled, 2, raw = TRUE)             5.5216  2    0.06324 .  
#   treatment                                          78.3598  1    < 2e-16 ***
#   poly(green_index_scaled, 2, raw = TRUE):treatment   1.9020  2    0.38636    

Anova(mv,type=3,component="zi")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: seed_half
# Chisq Df Pr(>Chisq)    
#   (Intercept)                                        0.1383  1   0.710010    
#   poly(green_index_scaled, 2, raw = TRUE)            8.7883  2   0.012349 *  
#   treatment                                         13.6098  1   0.000225 ***
#   poly(green_index_scaled, 2, raw = TRUE):treatment  9.7338  2   0.007697 ** 


summary(mv)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                                          1.78835    0.12708  14.073   <2e-16 ***
#   poly(green_index_scaled, 2, raw = TRUE)1             0.08110    0.10195   0.796   0.4263    
#   poly(green_index_scaled, 2, raw = TRUE)2            -0.11686    0.05599  -2.087   0.0369 *  
#   treatmentB                                          -0.97021    0.10960  -8.852   <2e-16 ***
#   poly(green_index_scaled, 2, raw = TRUE)1:treatmentB  0.16698    0.12138   1.376   0.1689    
#   poly(green_index_scaled, 2, raw = TRUE)2:treatmentB  0.02544    0.07271   0.350   0.7264    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                                          0.062715   0.168661   0.372 0.710010    
#   poly(green_index_scaled, 2, raw = TRUE)1            -0.464498   0.167266  -2.777 0.005486 ** 
#   poly(green_index_scaled, 2, raw = TRUE)2             0.008501   0.096242   0.088 0.929614    
#   treatmentB                                           0.772039   0.209273   3.689 0.000225 ***
#   poly(green_index_scaled, 2, raw = TRUE)1:treatmentB  0.597293   0.211877   2.819 0.004817 ** 
#   poly(green_index_scaled, 2, raw = TRUE)2:treatmentB -0.144852   0.149233  -0.971 0.331726 

###############
## Visualize ##
###############

# CIs working for ggeffects 

vis_v_all <-
  ggpredict(mv,
            terms = c("green_index_scaled[all]","treatment"),
            type = "zi_prob"); plot(vis_v_all)

# for final plot:
vis_spp_v <-
  ggpredict(mv,
            terms = c("green_index_scaled[all]", "treatment"),
            type = "fe.zi", allow.new.levels=TRUE); plot(vis_spp_v)


# add color as column for plotting
seed_vulpia$colortreat <- NA
seed_vulpia$colortreat[seed_vulpia$treatment == "A"]<- "steelblue3"
seed_vulpia$colortreat[seed_vulpia$treatment == "B"]<- "red4"
seed_vulpia$colortreat <- as.factor(seed_vulpia$colortreat)
colortreatraw <- c(seed_vulpia$colortreat)
length(vis_spp_v$x)  #746/2 # length of vis_spp_v
colortreatpred <- c(rep("steelblue3", times = 373), rep("red4", times = 373))

z <- ggplot() +
      geom_jitter(
        data = seed_vulpia,
        aes(green_index_scaled, seed_half, color = treatment),
        alpha = 0.1,
        size = 1.5,
        color = colortreatraw
      ) +
      geom_line(data = vis_spp_v, aes(x = x, y = predicted, group = group, color = group),
                linewidth = 1,
                color = colortreatpred
      ) +
      geom_ribbon(
        data = vis_spp_v,
        aes(
          x = x,
          y = predicted,
          group=group,
          ymin = conf.low,
          ymax = conf.high
        ),
        alpha = 0.1,
        show.legend = F,
        fill = colortreatpred
      ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  ylim(0, 15) +
  ggtitle("Festuca") + 
  annotate("text", x = -4, y = 15, label = "D", size = 7)
z


#-----------------------------------------------------------------------------.

# FULL FIGURE

library(patchwork)
pdf("Figures/spp_specific.pdf",
    width = 12,
    height = 6)
(x+w+y+z+e+f+g+h+plot_layout(ncol = 4)) # made micropus first
dev.off()


png("Figures/spp_specific.png",
    width = 12,
    height = 6,
    units = "in",
    res = 600)
(x+w+y+z+e+f+g+h+plot_layout(ncol = 4)) # made micropus first
dev.off()

#--------------------------------------------------------------------------------.

# ZERO-INFLATION SUPP FIGURE

colortreatpred <- c(rep("steelblue3", times = 372), rep("red4", times = 372))
zm <- ggplot() +
  geom_line(data = vis_m_all, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            #color = colortreatpred, 
            show.legend = T
  ) +
  scale_color_manual(values =c("steelblue3","red4"), 
                     labels = c("Without neighbors","With neighbors")) +
  geom_ribbon(
    data = vis_m_all,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.49,0.7),
        text = element_text(size = 16)) +
  labs(y = "Proportion of zeros", x = "") +
  ylim(0, 1) +
  ggtitle("Micropus") + 
  annotate("text", x = -4, y = 1, label = "A", size = 7)
zm

colortreatpred <- c(rep("steelblue3", times = 355), rep("red4", times = 355))
zp <- ggplot() +
  geom_line(data = vis_p_all, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_p_all,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  ylim(0, 1) +
  ggtitle("Plantago") + 
  annotate("text", x = -4, y = 1, label = "B", size = 7)
zp

colortreatpred <- c(rep("steelblue3", times = 351), rep("red4", times = 351))
zb <- ggplot() +
  geom_line(data = vis_b_all, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_b_all,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  ylim(0, 1) +
  ggtitle("Bromus") + 
  annotate("text", x = -4, y = 1, label = "C", size = 7)
zb

colortreatpred <- c(rep("steelblue3", times = 373), rep("red4", times = 373))
zv <- ggplot() +
  geom_line(data = vis_v_all, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_v_all,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "") +
  ylim(0, 1) +
  ggtitle("Festuca") + 
  annotate("text", x = -4, y = 1, label = "D", size = 7)
zv

#-----------------------------------------------------------------------------

#SEED PRODUCTION ALONE SUPP FIGURE (fe)


# enter new vis predict:

vis_p <-
  ggpredict(mp,
            terms = c("green_index_scaled[all]", "treatment"),
            type = "fe", allow.new.levels=TRUE)

vis_m <-
  ggpredict(mm,
            terms = c("green_index_scaled[all]", "treatment"),
            type = "fe", allow.new.levels=TRUE)

vis_b <-
  ggpredict(mb,
            terms = c("green_index_scaled[all]", "treatment"),
            type = "fe", allow.new.levels=TRUE)

vis_v <-
  ggpredict(mv,
            terms = c("green_index_scaled[all]", "treatment"),
            type = "fe", allow.new.levels=TRUE)


# For plotting, set max of confidence interval to the maximum of our y-axis 
vis_m$conf.low[vis_m$conf.low > 15] <- 15
vis_m$conf.high[vis_m$conf.high > 15] <- 15

colortreatpred <- c(rep("steelblue3", times = 372), rep("red4", times = 372))
fzm <- ggplot() +
  geom_line(data = vis_m, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            #color = colortreatpred
            show.legend = T
  ) +
  scale_color_manual(values =c("steelblue3","red4"), 
                     labels = c("Without neighbors","With neighbors")) +
  geom_ribbon(
    data = vis_m,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'none',
        text = element_text(size = 16)) +
  labs(y = "Seed production", x = "Vegetation index (G-R)") +
  ylim(0, 15) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  annotate("text", x = -4, y = 15, label = "E", size = 7)
fzm

colortreatpred <- c(rep("steelblue3", times = 355), rep("red4", times = 355))
fzp <- ggplot() +
  geom_line(data = vis_p, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_p,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "Vegetation index (G-R)") +
  ylim(0, 15) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  annotate("text", x = -4, y = 15, label = "F", size = 7)
fzp

colortreatpred <- c(rep("steelblue3", times = 351), rep("red4", times = 351))
fzb <- ggplot() +
  geom_line(data = vis_b, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_b,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "Vegetation index (G-R)") +
  ylim(0, 15) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  annotate("text", x = -4, y = 15, label = "G", size = 7)
fzb

colortreatpred <- c(rep("steelblue3", times = 373), rep("red4", times = 373))
fzv <- ggplot() +
  geom_line(data = vis_v, aes(x = x, y = predicted, group = group, color = group),
            linewidth = 1,
            color = colortreatpred
  ) +
  geom_ribbon(
    data = vis_v,
    aes(
      x = x,
      y = predicted,
      group=group,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.1,
    show.legend = F,
    fill = colortreatpred
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        text = element_text(size = 16)) +
  labs(y = "", x = "Vegetation index (G-R)") +
  ylim(0, 15) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "black") +
  annotate("text", x = -4, y = 15, label = "H", size = 7)
fzv


pdf("Figures/zero-inflation.pdf",
    width = 12,
    height = 6)
(zm + zp + zb + zv + fzm + fzp + fzb + fzv + plot_layout(ncol = 4))
dev.off()


png("Figures/zero-inflation.png",
    width = 12,
    height = 6,
    units = "in",
    res = 600)
(zm + zp + zb + zv + fzm + fzp + fzb + fzv + plot_layout(ncol = 4))
dev.off()

