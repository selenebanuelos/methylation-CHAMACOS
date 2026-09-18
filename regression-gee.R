# Author: Selene Banuelos
# Date: 9/10/2026
# Description: Fit mixed effects model using generalized estimating equations 
# to estimate unadjusted population average DNAm age deviation between asthma
# trajectory classes

# DNAm age ~ chrono age + asthma traj class + chrono age*asthma traj class

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(gee)
# function to get post-estimate associations from GEE (by Alan Hubbard)
#source('code/gee_post_estimate.R')

# import data ------------------------------------------------------------------
# participant information
demo <- read.dta13('data-raw/de_la_Rosa_07.dta',
                   nonint.factors = TRUE,
                   generate.factors = TRUE)

# estimated epigenetic age from various clocks
dnam_age <- read.dta13("data-raw/de_la_Rosa_epigenetic_07.dta", 
                       nonint.factors=TRUE, 
                       generate.factors=TRUE)

# asthma classification from latent class variable analysis
asthma <- read.csv('data-raw/class-labels-k3.csv')

# data wrangling ---------------------------------------------------------------
# vector of non-DNAm age variable names
non_DNAm <- c('newid', 'Chrono_Age', 'Timepoint', 'Array', 'Method')

# make DNAm age data longer for easier data wrangling
dnam_long <- dnam_age %>%
  # pivot all DNAm age columns into longer format
  pivot_longer(cols = -any_of(non_DNAm), 
               names_to = 'clock',
               values_to = 'DNAm_age')

# check chrono age of participants ---------------------------------------------
# looks like some participants might have 2 chrono ages for same timepoint?
ages_two <- dnam_age %>%
  select(newid, Chrono_Age, Timepoint, Array) %>%
  # keep only distinct rows
  distinct(.) %>%
  # identify participants that have more than 1 row per timepoint
  group_by(newid, Timepoint) %>%
  summarise(age_count = n_distinct(Chrono_Age)) %>%
  filter(age_count >1)

# get dataframe of participants that have 2 chrono ages per timepoint
df <- ages_two %>%
  select(-age_count) %>%
  # join original DNAm age data to list of participants with >1 chrono age
  left_join(., dnam_age, by = c('newid', 'Timepoint'))
  

# check that all participants are included in all methods ----------------------
# interested in methods: clock foundation, methscore CpG, methscore PC
check <- dnam_long %>%
  filter(Method %in% c('Clock Foundation', 'Methscore CpG', 'Methscore PC')) %>%
  pivot_wider(names_from = c(clock, Method),
              names_glue = '{clock}_{Method}',
              values_from = DNAm_age)

# ------------------------------------------------------------------------------
# use GrimAge and IEAA calculated using Clock Foundation method
clock_found <- filter(dnam_long, Method == 'Clock Foundation',
                      clock %in% c('GrimAge', 'IEAA')) %>%
  # specify that these are CpG-based clocks
  mutate(type = 'CpG')

methscore_cpg <- filter(dnam_long, Method == 'Methscore CpG',
                        clock %in% c('Horvath', 
                                     'SkinBlood', 
                                     'Hannum', 
                                     'PhenoAge')) %>%
  # specify that these are CpG-based clocks
  mutate(type = 'CpG')

methscore_pc <- filter(dnam_long, Method == 'Methscore PC',
                       clock %in% c('Horvath', 
                                    'SkinBlood', 
                                    'Hannum', 
                                    'PhenoAge',
                                    'GrimAge')) %>%
  # specify that these are PC-based clocks
  mutate(type = 'PC')

# merge selected clocks
select_dnam <- rbind(clock_found, methscore_cpg, methscore_pc) %>%
  # remove DNAm age calculation method & Array
  #select(-c(Method)) %>%
  # make data wider
  pivot_wider(names_from = c(clock, type),
              names_glue = '{clock}_{type}',
              values_from = DNAm_age)

# START HERE -------------------------------------------------------------------
# merge epigenetic age and asthma classifications
methscore_cpg <- dnam_age %>%
  # keep DNAm age generated with method of interest
  filter(Method == 'Methscore CpG') %>%
  # change variable type for merging
  mutate(newid = as.integer(newid)) %>%
  # merge asthma class to estimated epi ages
  left_join(., asthma, by = 'newid') %>%
  # keep only timepoints of interest (9Y-18Y)
  filter(Timepoint %in% (c('Age 9', 'Age 12', 'Age 14', 'Age 18'))) %>%
  # remove anyone missing asthma trajectory classification
  filter(!is.na(LCA)) 
# need to double check that everyone only has 1 measure per time

# factor variables
methscore_cpg$LCA <- factor(methscore_cpg$LCA,
                           levels = c('never/infrequent', # reference: N/I
                                      'late onset', 
                                      'persistent'))

# merge epigenetic age and asthma classifications
methscore_pc <- dnam_age %>%
  # keep DNAm age generated with method of interest
  filter(Method == 'Methscore PC') %>%
  # change variable type for merging
  mutate(newid = as.integer(newid)) %>%
  # merge asthma class to estimated epi ages
  left_join(., asthma, by = 'newid') %>%
  # keep only timepoints of interest (9Y-18Y)
  filter(Timepoint %in% (c('Age 9', 'Age 12', 'Age 14', 'Age 18'))) %>%
  # remove anyone missing asthma trajectory classification
  filter(!is.na(LCA)) 
# need to double check that everyone only has 1 measure per time

# factor variables
methscore_pc$LCA <- factor(methscore_pc$LCA,
                            levels = c('never/infrequent', # reference: N/I
                                       'late onset', 
                                       'persistent'))

# regression using GEE ---------------------------------------------------------

# should i be using exchangeable + chronological age: less efficient than AR?

# or autoregressive + timepoint (rounded age): will result in potentially larger
# epigenetic age deviation than if i just used accurate chronological age. Can
# only use timepoint here because AR structure assumes equally spaced time points

horvath_cpg <- gee(Horvath ~ Chrono_Age + LCA,
               id = newid,
               data = methscore_cpg,
               family = gaussian,
               corstr = 'exchangeable')

skinblood_cpg <- gee(SkinBlood ~ Chrono_Age + LCA,
                   id = newid,
                   data = methscore_cpg,
                   family = gaussian,
                   corstr = 'exchangeable')

hannum_cpg <- gee(Hannum ~ Chrono_Age + LCA,
                   id = newid,
                   data = methscore_cpg,
                   family = gaussian,
                   corstr = 'exchangeable')

pheno_cpg <- gee(PhenoAge ~ Chrono_Age + LCA,
                   id = newid,
                   data = methscore_cpg,
                   family = gaussian,
                   corstr = 'exchangeable')

horvath_pc <- gee(Horvath ~ Chrono_Age + LCA,
                   id = newid,
                   data = methscore_pc,
                   family = gaussian,
                   corstr = 'exchangeable')

skinblood_pc <- gee(SkinBlood ~ Chrono_Age + LCA,
                     id = newid,
                     data = methscore_pc,
                     family = gaussian,
                     corstr = 'exchangeable')

hannum_pc <- gee(Hannum ~ Chrono_Age + LCA,
                  id = newid,
                  data = methscore_pc,
                  family = gaussian,
                  corstr = 'exchangeable')

pheno_pc <- gee(PhenoAge ~ Chrono_Age + LCA,
                 id = newid,
                 data = methscore_pc,
                 family = gaussian,
                 corstr = 'exchangeable')

# 95% confidence intervals -----------------------------------------------------
# function that creates 95% CI using robust SE
robust_ci <- function(model, # gee object
                      label # character string describing model
                      ){
  
  # extract coefficients (1) and robust SE (4)
  coef_data <- summary(model)$coefficients[, c(1,4)]
  
  # calculate bounds
  lower_bound <- coef_data[,1] - 1.96 * coef_data[,2]
  upper_bound <- coef_data[,1] + 1.96 * coef_data[,2]
  
  # format 95% CI with estimate and bounds
  ci_table <- data.frame(
    Estimate = round(coef_data[,1], digits = 4),
    Lower_95_CI = round(lower_bound, digits = 4),
    Upper_95_CI = round(upper_bound, digits = 4),
    clock = label
  )
  
  print(paste('95% CIs for', label))
  return(ci_table)
  
}

results <- rbind(robust_ci(horvath_cpg, 'horvath cpg'),
      robust_ci(skinblood_cpg, 'skinblood cpg'),
      robust_ci(hannum_cpg, 'hannum cpg'),
      robust_ci(pheno_cpg, 'phenoage cpg'),
      robust_ci(horvath_pc, 'horvath pc'),
      robust_ci(skinblood_pc, 'skinblood pc'),
      robust_ci(hannum_pc, 'hannum pc'),
      robust_ci(pheno_pc, 'phenoage pc')
      )

# output -----------------------------------------------------------------------
write.csv(results, 'data-processed/gee-results.csv', row.names = TRUE)