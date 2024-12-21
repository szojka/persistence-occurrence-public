################################################################################
# CODE ASSOCIATED WITH MANUSCRIPT: "Sink populations decouple species occupancy and persistence across a productivity gradient"

## AUTHORS: Emilie F. Craig, Megan Szojka, Rachel M. Germain, Lauren G. Shoemaker 
################################################################################

# OVERVIEW:

This project addressed how persistence and occurrence of four species change across a productivity gradient in Californian serpentine grasslands. These four species were analyzed independently in the population level results (Fig 2, S6) and were aggregated for community levels results (Fig 3, 4). We tested the effects of environmental filtering and biotic interactions by clearing vegetated materials out of a plot, while leaving an adjacent plot intact creating a paired design (n paired plots = 450). We additionally tested the effects of dispersal by comparing species' natural occupancy patterns to transplant persistence within both cleared and intact plots. The naming of scripts below indicate either which organismal level (i.e., spp or community) or biological process (e.g., dispersal) is being analyzed within the script.

# HOW TO RUN PROJECT:

All scripts that start with 'Stats & figures...' depend on 'Source fitness data.R'. All the 'Stats & figures...' scripts are independent, save 'Stats & figures - spp occurrence patterns.R' and 'Stats & figures - spp seed production.R', as the first script is necessary to create Figure 2 which is coded in the latter script. Specific descriptions of each script's function are found below.

# DESCRIPTION OF RELEVENT SCRIPTS:

## Model_simulations.R

Running is not necessary. This script was used to simulate our nested random effect structure in order to make sure the random effects are specified in the most statistically correct fashion. 

## Source fitness data.R

Reading in and cleaning all the data necessary for figures. This needs to be sourced before every other script--this source code is included at the top of each script.

## Stats & figures - community models.R

Used to create Figures 3, Tables S3, S5, S6, and associated occurrence & persistence models. 

## Stats & figures - cover estimates.R

Used to create supplementary materials for Figures S3, S4, specifically to model the percent cover estimates between cleared vs uncleared paired plots across the productivity gradient.

## Stats & figures - dispersal effects.R

Used to create Figure 4, Table S7, and associated multinomial model.

## Stats & figures - spp occurrence patterns.R

Used to create Figure 2 and associated binomial occurrence models for each species.

## Stats & figures - spp seed production.R

Used to create Figures 2, S6, Tables S2, S4 and associated species specific seed production, population growth, and zero-inflation models.

## Stats & figures - pseudo sinks method 1.R

Used to create Figure S7, looking at how abundance categories for each species predict seed production, as evidence for or against density dependence.

## Stats & figures - pseudo sinks method 2.R

Used to create Figure S8, looking at how abundance categories for each species predict the proportion of sinks, potentially showing evidence for or against pseudo-sinks.

# DESCRIPTION OF DATA FILES:

## Transplant_master data_Oct2020.csv

Seed production data for each transplant. Used in script "Source fitness data.R".

## Tag labels.csv

Connects each transplant tag to the plot, grid, and site it is contained in. Used in script "Source fitness data.R".

## Abundance & occurrence 2019_final.csv

Observed abundance categories for each species in each plot. Used in script "Source fitness data.R".

## greenness_plot_data.csv

Connects plot identity to amount of green, blue, and red color indexes. Used to create vegetaion index (Green - Red). Used in script "Source fitness data.R".

## Plot photo cover estimates 2023.csv

Supplemental data that estimates percent cover of all vegetation for each plot. Used in script "Supp stats & figures - cover estimates.R".


# SUPPORTING INFORMATION AND CONTACT:

1. code DOI: 10.5281/zenodo.8381381

2. contact author Megan Szojka mszojka@uwyo.edu with any concerns or questions.

3. for use of these data in other projects, the authors request that you please let us know, as there are other ongoing projects with these data.

##################################################################


