
##################################
# SUPPLEMENTAL MATERIAL
##################################

# DESCRIPTION:
# Supplementary material to see how effectively we reduced the cover in cleared plots in comparison to the intact paired plot. Also to see how total % cover changed across the productivity gradient.

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
# download read in csv version of data

cover1 <- read_csv("Data/Plot photo cover estimates 2023.csv")

# rename some things
names(cover1)[1] <- "plot"
names(cover1)[2] <- 'image'
names(cover1)[3] <- "treatment"
cover1 <- cover1[,-8]

cover1$plot <- as.factor(cover1$plot)
cover1$treatment <- as.factor(cover1$treatment)
cover1$image <- as.factor(cover1$image)

# summary stats:

cover <- cover1 %>%
  pivot_longer(., cols = 4:7, names_to = "category", values_to = 'cover')
cover$category <- as.factor(cover$category)

# magnitudes of difference in cover
cover %>%
  group_by(treatment, category) %>%
  dplyr::summarize(min = min(cover), max = max(cover))

# magnitudes of difference in productivity
green_all %>%
  group_by(treatment) %>%
  dplyr::summarize(min = min(green_index), max = max(green_index))

# total cover includes litter this year, last year, and live biomass
total_cover <- cover1
total_cover$total <- NA
total_cover$total <- rowSums(total_cover[5:7])

# Link total_cover with another key to get grid using image and plot.
key_for_cover <- green_all %>%
  dplyr::select(image, grid, site, plot, green_index_scaled, block) %>%
  distinct()
key_for_cover$image <- as.factor(key_for_cover$image) 

# link
cover_green_all <- left_join(total_cover, key_for_cover, by = c("plot", "image"),
                             relationship = "many-to-many")
cover_green_all # doesn't have plot 341, 343 or 345 

# to use beta family need 0 < y < 1
cover_green_all$prop <- NA
cover_green_all$prop <- cover_green_all$total/100

cover_green_all$prop[cover_green_all$prop == 0] <- 0.01
cover_green_all$prop[cover_green_all$prop >= 1] <- 0.99

cover_green_all <- na.omit(cover_green_all) # can't fix missing plots n was 892 -> 850

########################
# MODELLING

# relationship between total cover by treatment across productivity

mod <- glmmTMB(prop ~ green_index_scaled*treatment + (1|site) + (1|site:grid), 
               data = cover_green_all,
               family = beta_family()) # having block created singular fit
 
summary(mod)

vis1 <- ggpredict(mod, terms = c('green_index_scaled[all]', 'treatment')); plot(vis1)

supp_cover <- ggplot() +
  geom_jitter(cover_green_all, mapping = aes(x = green_index_scaled, y = prop, color = treatment), alpha = 0.3 ) +
  geom_line(vis1, mapping = aes(x = x, y = predicted, color = group),
            linewidth = 1) +
  geom_ribbon(vis1, mapping = aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
              show.legend = FALSE,
              alpha = 0.3) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'top',
        text = element_text(size = 16)) +
  labs(x = "Vegetation index (G-R)", y = "Total proportion cover") +
  scale_color_manual(values = c("sienna1","darkgreen"), labels = c("Without neighbors", "With neighbors")) +
  scale_fill_manual(values = c("sienna1","darkgreen"), labels = c("Without neighbors", "With neighbors")) +
  ylim(0, 1)

pdf(
  'Figures/supp_cover.pdf',
  width = 5,
  height = 5
)
supp_cover
dev.off()

png(
  'Figures/supp_cover.png',
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
supp_cover
dev.off()

###############

# pairwise reduction stats (new column % reduced = total biotic - total abiotic cover)
new_cover_diff <- cover_green_all %>%
  select(-Bare, -Live, -`Old litter`, -Dead, -prop) %>%
  pivot_wider(., names_from = treatment, values_from = total) %>%
  mutate(reduced = B - A) %>%
  mutate(n = n())

# summarize average reduction in cover
new_cover_diff %>%
  dplyr::summarize(min = min(reduced), max = max(reduced), average = mean(reduced), sd = sd(reduced))

# find SE
new_cover_diff %>%
  mutate(se = sd(reduced)/sqrt(n)) %>%
  select(se) %>%
  distinct() # 0.885

# relationship between our reduction of biomass across productivity
mod2 <- lmer(reduced ~ green_index_scaled + (1|site) + (1|site:grid), data = new_cover_diff) # block irrelevant

summary(mod2)

vis <- ggpredict(mod2, terms = c('green_index_scaled[all]')); plot(vis)

supp_reduced <- ggplot() +
  geom_jitter(new_cover_diff, mapping = aes(x = green_index_scaled, y = reduced),
              color = 'midnightblue', alpha = 0.5) +
  geom_line(vis, mapping = aes(x = x, y = predicted),
            color = 'midnightblue') +
  geom_ribbon(vis, mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = 'midnightblue',
              alpha = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  labs(x = "Vegetation index (G-R)", y = "Reduced % cover for blocks") +
  geom_hline(yintercept = 0, linetype = 'dashed')
supp_reduced

pdf(
  'Figures/supp_reduced.pdf',
  width = 5,
  height = 5
)
supp_reduced
dev.off()

png(
  'Figures/supp_reduced.png',
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
supp_reduced
dev.off()



