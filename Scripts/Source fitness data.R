
###########################################
## Source code for clean data files
###########################################

# if running source script independently, load these packages:

# library(tidyverse)
# library(stringr)
# library(car)
# library(lme4)
# library(glmmTMB)
# library(visreg)
# library(patchwork)
# library(ggeffects)

############ Clean persistence dataset

seeds_unclean <- read_csv("Data/Transplant_master data_Oct2020.csv")

seeds_unclean$status <- as.factor(as.character(seeds_unclean$status))
levels(seeds_unclean$status) # four levels

seeds_unclean$status <- as.character(as.factor(seeds_unclean$status)) # make status = NA into factor
seeds_unclean$status[is.na(seeds_unclean$status)] <- c("normal")
seeds_unclean <- filter(seeds_unclean, !status%in%c('missing','herbivory')) # remove uninformative transplants
seeds_unclean <- na.omit(seeds_unclean)


# summarize by scale to get means, n, and cis grouped by treatments
seeds_clean <- ungroup(seeds_unclean)

############ Clean greenness data (from image J Fiji)

greenness_plot_data <- read_csv("Data/greenness_plot_data.csv") # updeated 2023
green <- dplyr::select(greenness_plot_data, grid, site, plot, blue_mean, red_mean, green_mean) 
rm(greenness_plot_data)

############ Clean abundance & occurrence data 

# join abundance data
# abundance categories = 0,1,2,3
# key = 0, 1-10, 11-100, 101-1000

photos1 <- read_csv("Data/Abundance & occurrence 2019_final.csv")
photos1 <- select(photos1, -notes)
photos <- gather(photos1, species, ab_cat, 6:9) # abundance categories
photos <- na.omit(photos)
# join tags to photos using:
tag_labels <- read_csv("Data/Tag labels.csv") 

ab <- left_join(tag_labels, photos, by = c("species", "treatment", "grid", "site", "plot")) 

green_ab <- left_join(green, ab, c("site","grid","plot")) 

# join with seeds_clean
green_seeds <- left_join(green_ab,seeds_clean, by =c("site","grid","treatment","species","tag"))
green_seeds <- filter(na.omit(green_seeds))

############ Assign same abundance in treatment A as paired treatment B

plots <- levels(factor(tag_labels$plot))
grids <- as.numeric(levels(factor(tag_labels$grid)))
spp <- levels(factor(tag_labels$species))

green_seeds$occ <- green_seeds$ab_cat

green_seeds %>% filter(species == "plaere") %>% na.omit(green_seeds) -> green_seeds1

for(j in grids) {
  plots <- levels(factor(green_seeds1$plot[green_seeds1$grid == j]))
  for(i in plots) {
    if(length(green_seeds1$ab_cat[green_seeds1$treatment %in% c("A","B") & green_seeds1$plot == i & green_seeds1$grid == j]) == 2) { 
      # first filter for observations I have where A and B exist
      if(green_seeds1$ab_cat[green_seeds1$treatment == "A" & green_seeds1$plot == i & green_seeds1$grid == j] <= green_seeds1$ab_cat[green_seeds1$treatment == "B" & green_seeds1$plot == i & green_seeds1$grid == j]){
        
        green_seeds1$occ[green_seeds1$treatment == "A" & green_seeds1$plot == i & green_seeds1$grid == j] <- green_seeds1$occ[green_seeds1$treatment == "B" & green_seeds1$plot == i & green_seeds1$grid == j] 
        # if A occurrence is < B occurrence, assign A to occurrence B's.
        
      }
      if(green_seeds1$ab_cat[green_seeds1$treatment == "B" & green_seeds1$plot == i & green_seeds1$grid == j] <= green_seeds1$ab_cat[green_seeds1$treatment == "A" & green_seeds1$plot == i & green_seeds1$grid == j]){
        
        green_seeds1$occ[green_seeds1$treatment == "B" & green_seeds1$plot == i & green_seeds1$grid == j] <- green_seeds1$occ[green_seeds1$treatment == "A" & green_seeds1$plot == i & green_seeds1$grid == j] 
        # if B occurrence is < A occurrence, assign B occurrence to A's.
      }
    }
  }
}

# run loop again for each species then rbind
green_seeds %>% filter(species == "miccal") %>% na.omit(ab_cat) -> green_seeds2

for(j in grids) {
  plots <- levels(factor(green_seeds2$plot[green_seeds2$grid == j]))
  for(i in plots) {
    if(length(green_seeds2$ab_cat[green_seeds2$treatment %in% c("A","B") & green_seeds2$plot == i & green_seeds2$grid == j]) == 2) {
      if(green_seeds2$ab_cat[green_seeds2$treatment == "A" & green_seeds2$plot == i & green_seeds2$grid == j] <= green_seeds2$ab_cat[green_seeds2$treatment == "B" & green_seeds2$plot == i & green_seeds2$grid == j]){
        
        green_seeds2$occ[green_seeds2$treatment == "A" & green_seeds2$plot == i & green_seeds2$grid == j] <- green_seeds2$occ[green_seeds2$treatment == "B" & green_seeds2$plot == i & green_seeds2$grid == j]
      }
      if(green_seeds2$ab_cat[green_seeds2$treatment == "B" & green_seeds2$plot == i & green_seeds2$grid == j] <= green_seeds2$ab_cat[green_seeds2$treatment == "A" & green_seeds2$plot == i & green_seeds2$grid == j]){
        
        green_seeds2$occ[green_seeds2$treatment == "B" & green_seeds2$plot == i & green_seeds2$grid == j] <- green_seeds2$occ[green_seeds2$treatment == "A" & green_seeds2$plot == i & green_seeds2$grid == j] 
        
      }
    }
  }
}


green_seeds %>% filter(species == "brohor") %>% na.omit(ab_cat) -> green_seeds3

for(j in grids) {
  plots <- levels(factor(green_seeds3$plot[green_seeds3$grid == j]))
  for(i in plots) {
    if(length(green_seeds3$ab_cat[green_seeds3$treatment %in% c("A","B") & green_seeds3$plot == i & green_seeds3$grid == j]) == 2) {
      # if(green_seeds3$ab_cat[green_seeds3$treatment == "A" & green_seeds3$plot == i & green_seeds3$grid == j] <= green_seeds3$ab_cat[green_seeds3$treatment == "B" & green_seeds3$plot == i & green_seeds3$grid == j]){
      #
      #   green_seeds3$occ[green_seeds3$treatment == "A" & green_seeds3$plot == i & green_seeds3$grid == j] <- green_seeds3$occ[green_seeds3$treatment == "B" & green_seeds3$plot == i & green_seeds3$grid == j]
      # }
      if(green_seeds3$ab_cat[green_seeds3$treatment == "B" & green_seeds3$plot == i & green_seeds3$grid == j] <= green_seeds3$ab_cat[green_seeds3$treatment == "A" & green_seeds3$plot == i & green_seeds3$grid == j]){

        green_seeds3$occ[green_seeds3$treatment == "B" & green_seeds3$plot == i & green_seeds3$grid == j] <- green_seeds3$occ[green_seeds3$treatment == "A" & green_seeds3$plot == i & green_seeds3$grid == j]

       }
    }
  }
}

green_seeds %>% filter(species == "vulmic") %>% na.omit(ab_cat) -> green_seeds4

for(j in grids) {
  plots <- levels(factor(green_seeds4$plot[green_seeds4$grid == j]))
  for(i in plots) {
    if(length(green_seeds4$ab_cat[green_seeds4$treatment %in% c("A","B") & green_seeds4$plot == i & green_seeds4$grid == j]) == 2) {
      if(green_seeds4$ab_cat[green_seeds4$treatment == "A" & green_seeds4$plot == i & green_seeds4$grid == j] <= green_seeds4$ab_cat[green_seeds4$treatment == "B" & green_seeds4$plot == i & green_seeds4$grid == j]){
        
        green_seeds4$occ[green_seeds4$treatment == "A" & green_seeds4$plot == i & green_seeds4$grid == j] <- green_seeds4$occ[green_seeds4$treatment == "B" & green_seeds4$plot == i & green_seeds4$grid == j]
      }
      if(green_seeds4$ab_cat[green_seeds4$treatment == "B" & green_seeds4$plot == i & green_seeds4$grid == j] <= green_seeds4$ab_cat[green_seeds4$treatment == "A" & green_seeds4$plot == i & green_seeds4$grid == j]){
        
        green_seeds4$occ[green_seeds4$treatment == "B" & green_seeds4$plot == i & green_seeds4$grid == j] <- green_seeds4$occ[green_seeds4$treatment == "A" & green_seeds4$plot == i & green_seeds4$grid == j] 
        
      }
    }
  }
}


green_all <- rbind(green_seeds1, green_seeds2, green_seeds3, green_seeds4)

rm(photos1,green_seeds1, green_seeds2, green_seeds3, green_seeds4, green_seeds, green)

############ Calculate occurrence and persistence thresholds 

green_all$persistence <- NULL
green_all$persistence[green_all$seed >= 2] <- 1
green_all$persistence[green_all$seed < 2 ] <- 0
green_all$occurrence <- NULL
green_all$occurrence[green_all$occ > 0] <- 1
green_all$occurrence[green_all$occ == 0] <- 0
# add an informative name for occ (column where abundance categories are intact and shared within block)
names(green_all)[14] <- "ab_cat_shared"

############  Calculate greenness index (Red - Green) 

green_all$green_index <- green_all$green_mean - green_all$red_mean
# REVERSE ORDER 
green_all$green_index <- green_all$green_index*(-1)
# this makes productive areas higher numbers :) 

# find number of species with viable data in each plot for future analyses: 'replicates'
rep_dat <- green_all %>% 
  select(plot, grid, species, treatment) %>%
  group_by(grid, plot, treatment) %>%
  mutate(replicates = n()) %>%
  select(-species) %>%
  distinct()
max(rep_dat$replicates) # max is 4 -> good! 

# turn necessary variables to factors
rep_dat$grid <- as.factor(rep_dat$grid)
rep_dat$plot <- as.factor(rep_dat$plot)
rep_dat$treatment <- as.factor(rep_dat$treatment)

############ Categorize productivity

green_all$green_cat <- NA
green_all <- green_all[order(green_all$green_index),] # performs correct
# 2931 observations, split data frame into 3 categories: 1-977-1954,1955-2931
x1 <- filter(green_all[c(1:977),])
x2 <- filter(green_all[c(978:1954),])
x3 <- filter(green_all[c(1955:2931),])
x1$green_cat <- "Harsh"
x2$green_cat <- "Mid"
x3$green_cat <- "Prod"

green_all <- rbind(x1,x2,x3)
unique(green_all$green_cat)

############ Data for persistence & occurrence figure

temp_p <- green_all %>% # temporary data frame to manipulate persistence (yes or no)
  dplyr::select(plot, species, treatment, persistence, green_index, green_cat, grid, site) %>%
  dplyr::distinct()
temp_p <- pivot_wider(temp_p, names_from = "species", values_from = "persistence")
temp_p[is.na(temp_p)] <- 0
temp_p$type <- "P"
temp_p$species_no <- rowSums(temp_p[c(7:10)])

temp_o <- green_all %>% # temporary data frame to manipulate occurrence (yes or no)
  dplyr::select(plot, species, treatment, occurrence, green_index, green_cat, grid, site) %>%
  dplyr::distinct()
temp_o <- pivot_wider(temp_o, names_from = "species", values_from = "occurrence")
temp_o[is.na(temp_o)] <- 0
temp_o$type <- "O"
temp_o$species_no <- rowSums(temp_o[c(7:10)])
fig_dat <- rbind(temp_p, temp_o) # now has persistence and occurrence as binary
rm(temp_o, temp_p)

############ Data for dispersal limitation vs sinks figure

# categories ME, SS_n, SS_y, DL
green_all$contingency <- with(green_all, paste0(occurrence, persistence))

green_all$contingency[green_all$contingency == "11"] <- c('SS_y') # species sorting
green_all$contingency[green_all$contingency == "00"] <- c('SS_n') # species sorting
green_all$contingency[green_all$contingency == "10"] <- c('ME') # mass effects aka sinks
green_all$contingency[green_all$contingency == "01"] <- c("DL") # dispersal limitation

############ Dataframe transformations for statistical analyses

# make characters factors

fig_dat$plot <- as.factor(fig_dat$plot)
fig_dat$grid <- as.factor(fig_dat$grid)
fig_dat$site <- as.factor(fig_dat$site)
fig_dat$treatment <- as.factor(fig_dat$treatment)
fig_dat$type <- as.factor(fig_dat$type)

green_all$grid <- as.factor(green_all$grid)
green_all$site <- as.factor(green_all$site)
green_all$plot <- as.factor(green_all$plot)
green_all$treatment <- as.factor(green_all$treatment)
green_all$contingency <- as.factor(green_all$contingency)

# add a replicates column on fig_dat and green_all for no_species with viable data in each plot
fig_dat <- left_join(fig_dat, rep_dat, by = c('plot','grid','treatment'))
green_all <- left_join(green_all, rep_dat, by = c('plot','grid','treatment'))

# Create scaled variables for stats analyses:

#1. scale x for whole dataset:
green_all$green_index_scaled <- NA
green_all$green_index_scaled <- scale(green_all$green_index)
# check that values are the same per plot:
check <- green_all %>%
  select(plot, grid, site, treatment, green_index_scaled) %>%
  distinct() %>%
  pivot_wider(., names_from = 'treatment', values_from = 'green_index_scaled')
rm(check) # good. 
# make green_index_scaled numeric
green_all$green_index_scaled <- as.numeric(green_all$green_index_scaled)

fig_dat$green_index_scaled <- NA
fig_dat$green_index_scaled <- scale(fig_dat$green_index)
# check that values are the same per plot:
check <- fig_dat %>%
  select(plot, grid, site, treatment, green_index_scaled) %>%
  distinct() %>%
  pivot_wider(., names_from = 'treatment', values_from = 'green_index_scaled')
rm(check) # good. 
# make green_index_scaled numeric
fig_dat$green_index_scaled <- as.numeric(fig_dat$green_index_scaled)

#2. log seed for exponential models with all species first
# when log transforming zeros, worth thinking about other statistical consequences of log(y + c) https://aosmith.rbind.io/2018/09/19/the-log-0-problem/
# in reality, the version of exponential using log(y) ~ x, family = Poisson, will not be mathematically different than y ~ x, poison. So don't stress over constant as model likely won't be chosen
green_all$seed_half <- NA
green_all$seed_half <- round(green_all$seed/2) # rounds seeds up to integer, not a problem bc this would only lead to overestimating persistence, which is a conservative approach
green_all$log_seed <- NA
green_all$log_seed <- log(green_all$seed_half + 1) 

# create block to name group of paired plots:

green_all$block <- NA
green_all$block <- paste0(green_all$grid, green_all$plot)

green_all$block <- as.factor(green_all$block)

fig_dat$block <- NA
fig_dat$block <- paste0(fig_dat$grid, fig_dat$plot)

fig_dat$block <- as.factor(fig_dat$block)

# checks that all loaded:
# for posting reproducable code start with these files
# View(green_all)
# View(fig_dat)
# View(fig_cat)


