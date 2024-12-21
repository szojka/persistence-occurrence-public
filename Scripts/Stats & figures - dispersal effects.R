
#------------------------------------------------------------.
# DISPERSAL EFFECTS ON POPULATIONS ####
#------------------------------------------------------------.

# DESCRIPTION:
# data wrangling to model a multinomial response of species sorting, sink populations, or dispersal limitation, as dependent outcomes across the productivity gradient. Visualize these modelled relationships.

library(tidyverse)
library(stringr)
library(car)
library(lme4)
library(glmmTMB)
library(visreg)
library(patchwork)
library(ggeffects)
library(DHARMa)
library(mclogit)
library(optimx)

source("Scripts/Source fitness data.R")

# Necessary functions:

# confidence intervals for proportion
# level = confidence level, here 0.975 = total or 95% CI
# n = sample size
# p = proportion
ci_prop <- function(level = 0.975, n, p) qt(level,df=n-1)*sqrt(p*(1-p)/n) 

# Take outputs that are alphas and turn them into probabilities spanning 0->1 
# j = coefficient counter
# J = number of coefficients -1, same as number of possible categories - 1
# a = output from multinomial model (pairwise comparisons)
prob_trans <- function(j, J=2) { # have to have J = 2 because indexing of alpha vector 
  exp(a[j])/(1 + sum(sapply(1:J, function(i) exp(a[i])))) 
}

# normal confidence intervals
# se = standard error
# level = confidence level, here 0.975 = total or 95% CI
# n = sample size
ci_norm <- function(level = 0.975, n, se)  qt(level,df=n-1)*se

##################
## Prepare data ##
##################

green_all_DLMESS <- green_all %>%
  dplyr::select(site, grid, block, species, treatment, persistence, occurrence, green_index_scaled, contingency, replicates) %>%
  distinct()
green_all_DLMESS$contingency <-
  recode_factor(green_all_DLMESS$contingency, "SS_n" = "SS")
green_all_DLMESS$contingency <-
  recode_factor(green_all_DLMESS$contingency, "SS_y" = "SS")
green_all_DLMESS <- green_all_DLMESS %>%
  filter(contingency %in% c('DL', "ME", 'SS')) %>%
  group_by(block, site, grid, contingency, treatment) %>%
  mutate(cont_freq = n()) %>%
  mutate(cont_prop = cont_freq / replicates) %>%
  dplyr::select(green_index_scaled, contingency, cont_prop, site, grid, replicates) %>%
  distinct() %>%
  pivot_wider(., names_from = 'contingency', values_from = 'cont_prop')

# NAs  in SS ME or DL are zeros
green_all_DLMESS[is.na(green_all_DLMESS)] <- 0

y_mat <- as.matrix(green_all_DLMESS[,c(7:9)])

##########################################
# need to calculate raw proportions for figure
green_all_fig3_DLMESS <- green_all
green_all_fig3_DLMESS$contingency <-
  recode_factor(green_all_fig3_DLMESS$contingency, "SS_n" = "SS")
green_all_fig3_DLMESS$contingency <-
  recode_factor(green_all_fig3_DLMESS$contingency, "SS_y" = "SS")
green_all_fig3_DLMESS <- green_all_fig3_DLMESS %>%
  filter(contingency %in% c('DL', "ME", 'SS')) %>%
  group_by(green_index_scaled, treatment) %>%
  mutate(total_freq = n()) %>%
  group_by(green_index_scaled, treatment, contingency) %>%
  mutate(cont_freq = n()) %>%
  mutate(cont_prop = cont_freq / total_freq) %>%
  dplyr::select(green_index_scaled, contingency, cont_prop, block, site, grid) %>%
  distinct() %>%
  pivot_wider(., names_from = 'contingency', values_from = 'cont_prop')
# add zeros
green_all_fig3_DLMESS$DL[is.na(green_all_fig3_DLMESS$DL)] <- 0
green_all_fig3_DLMESS$ME[is.na(green_all_fig3_DLMESS$ME)] <- 0
green_all_fig3_DLMESS$SS[is.na(green_all_fig3_DLMESS$SS)] <- 0
green_all_fig3_DLMESS <- green_all_fig3_DLMESS %>%
  pivot_longer(
    .,
    cols = c("DL", "ME", "SS"),
    names_to = 'contingency',
    values_to = 'cont_prop'
  )


##################
## Fit a model 
##################

#accounting for potential correlation among observations with random effects

mfit <- mblogit(y_mat ~ green_index_scaled*treatment,
                random = c(~ 1|site/grid, ~1|block), 
                weights = replicates, 
                data = green_all_DLMESS,
                control = mmclogit.control(maxit = 500)) # converged on last iteration

summary(mfit)

# Call:
#   mblogit(formula = y_mat ~ green_index_scaled * treatment, data = green_all_DLMESS, 
#           random = c(~1 | site/grid, ~1 | block), weights = replicates, 
#           control = mmclogit.control(maxit = 500))
# 
# Equation for ME vs SS:
#   Estimate Std. Error z value Pr(>|z|)   
#   (Intercept)                   -0.42453    0.19750  -2.150  0.03159 * 
#   green_index_scaled             0.13993    0.09110   1.536  0.12456   
#   treatmentB                     0.27205    0.08470   3.212  0.00132 **
#   green_index_scaled:treatmentB  0.17615    0.09562   1.842  0.06545 . 
# 
# Equation for DL vs SS:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                    -1.3385     0.1508  -8.876  < 2e-16 ***
#   green_index_scaled              0.1094     0.1032   1.060   0.2892    
#   treatmentB                     -0.5747     0.1316  -4.366 1.26e-05 ***
#   green_index_scaled:treatmentB  -0.2656     0.1163  -2.284   0.0224 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###################
## Visualize pt1 ##
###################

# no longer necessary as predict automatically uses dataframe inputed to model
new <- green_all_DLMESS %>%
  ungroup() %>%
  dplyr::select(treatment, green_index_scaled, site, block, grid) %>%
  arrange(green_index_scaled) # organize by green_index
length(new$treatment) # 846

# predict

# evidence predict should work on mclogit: https://www.elff.eu/software/mclogit/manual/predict/

# predict without considering random effects for smooth lines
p <- predict(mfit, type = 'response', conditional = FALSE) # predict.mmblogit initiated with se.fit = T

# Use predict to calculate SE incorporating random effect variance
p.se <- predict(mfit, type = 'response', se.fit = TRUE, conditional = TRUE)

# check predictions around zero
temp1 <- cbind(p, new$green_index_scaled) # things jump around a lot ~ 0 so trying conditional

# organize predictions
pred1 <- data.frame(p)
pred1 <- pivot_longer(pred1, cols = 1:3, names_to = 'response.level', values_to = 'predicted')

pred2 <- data.frame(p.se$se.fit)
pred2 <- pivot_longer(pred2, cols = 1:3, names_to = 'response.level', values_to = 'std.error')
pred2 <- pred2[, -1]
pred <- cbind(pred1, pred2)

# make sure category dataframe is structure long to be compatable with predictions data
new_rep <- cbind(new, temp1[,-4])
new_rep <- pivot_longer(new_rep, cols = 6:8) 
new_rep <- new_rep[,c(-6,-7)]

# new dataframe with predicted fit and se responses by categories

predict_new <- cbind(pred, new_rep)

# make confidence intervals from SE

n_treat <- green_all_DLMESS %>%
  group_by(treatment) %>%
  mutate(n = n()) %>%
  dplyr::select(n) %>%
  distinct() %>%
  print() # get frequency of A vs B from raw data

predict_new <- left_join(predict_new, n_treat, by = 'treatment')
predict_new$ci <- NA
predict_new$ci <- ci_norm(n = predict_new$n, se = predict_new$std.error)

###########################
# Find confidence intervals
###########################
{
predict_new$high <- NA
predict_new$high <-  predict_new$predicted +  predict_new$ci
predict_new$low <- NA
predict_new$low <-  predict_new$predicted -  predict_new$ci

# set up colors for ggplot
length(predict_new$response.level)/3
colortreatpred <- c(rep("red4", times = 846), 
                    rep("lavenderblush4", times = 846), 
                    rep("steelblue1", times = 846))

predict_new$colors <- NA
predict_new$colors[predict_new$response.level == "SS"] <- "lavenderblush4"
predict_new$colors[predict_new$response.level == "ME"] <- "red4"
predict_new$colors[predict_new$response.level == "DL"] <- "steelblue1"
color_sting <- c(predict_new$colors)

predict_new$treatment  <- recode_factor(predict_new$treatment , "A" = "Without neighbors")
predict_new$treatment  <- recode_factor(predict_new$treatment , "B" = "With neighbors")
levels(predict_new$treatment)
green_all_fig3_DLMESS$treatment  <- recode_factor(green_all_fig3_DLMESS$treatment , "A" = "Without neighbors")
green_all_fig3_DLMESS$treatment  <- recode_factor(green_all_fig3_DLMESS$treatment , "B" = "With neighbors")
levels(green_all_fig3_DLMESS$treatment)


# change factor order so that Species sorting is in middle
predict_new$response.level <- as.factor(predict_new$response.level)
levels(predict_new$response.level) # was DL ME SS
predict_new$response.level <- factor(predict_new$response.level, levels=c("DL", "SS", "ME"))
levels(predict_new$response.level) # now DL SS ME

green_all_fig3_DLMESS$contingency <- factor(green_all_fig3_DLMESS$contingency, levels=c("DL", "SS", "ME"))
levels(green_all_fig3_DLMESS$contingency) # now DL SS ME
}

###############
# summary stats for table

predict_new_new <- predict_new %>%
  select(response.level, treatment, predicted, high, low, green_index_scaled)

view(predict_new_new)

#################
# Visualize pt 2
#################

disp <- 
ggplot() +
  geom_jitter(
    data = green_all_fig3_DLMESS, # need dataframe with raw proportions
    aes(x = green_index_scaled, y = cont_prop, group = contingency, color = contingency), # shapes = disp lim = down arrow  (shape = 6), sink = uparrow (shape = 2), ss = square (shape = 0)
    alpha = 0.2,
    width = 0.01,
    height = 0.001,
    size = 1.5
  ) +
  scale_color_manual(
    values = c("red4", "lavenderblush4","steelblue1"),
    labels = c("Dispersal limitation", "Species sorting","Sink populations"),
  ) +
  scale_fill_manual(
    values = c( "lavenderblush4", "steelblue1","red4"),
    labels = c("Species sorting","Dispersal limitation", "Sink populations"),
    guide = "none"
  ) +
  geom_line(data = predict_new, mapping = aes(x = green_index_scaled, 
                                              y = predicted, 
                                              group = response.level, 
                                              color = response.level),
            show.legend = TRUE,
            linewidth = 1) + 
  geom_ribbon(data = predict_new, mapping = aes(x = green_index_scaled,
                                                y = predicted,
                                                fill = colors,
                                                group=response.level,
                                                ymin = predicted-ci,
                                                ymax = predicted+ci),
              alpha = 0.35,
              show.legend = FALSE
  ) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'top',
        text = element_text(size = 16))  +
  facet_wrap(vars(treatment)) +
  labs(y = "Proportion of species", x = "Vegetation index (G-R)") 
  # geom_text(
  #   data    = dat_text,
  #   mapping = aes(x = c(-4,-4), y = c(1,1), label = label)
  # )
 disp
 
dat_text <- data.frame(
  label = c("A", "B"),
  treatment   = c('With neighbors','Without neighbors'))
  
 
disp_full <- disp + geom_text(
  data    = dat_text,
  mapping = aes(x = -4, y = .97, label = label),
  hjust   = -0.1,
  vjust   = -0.1,
  size = 7)


pdf(
  'Figures/dispersal.pdf',
  width = 8,
  height = 5
)
disp_full
dev.off()

png(
  'Figures/dispersal.png',
  width = 8,
  height = 5,
  units = "in",
  res = 600
)
disp_full
dev.off()
