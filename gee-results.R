# Author: Selene Banuelos
# Date: 9/10/2026
# Description: Fit mixed effects model using generalized estimating equations 
# to estimate unadjusted population average DNAm age deviation between asthma
# trajectory classes

# DNAm age ~ chrono age + asthma traj class

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(gee)
library(purrr)
library(tibble)
# function to get post-estimate associations from GEE (by Alan Hubbard)
source('code/gee_post_estimate.R')

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
# timepoints of interest
timepoints <- c('Age 9', 'Age 12', 'Age 14', 'Age 18')

# vector of non-DNAm age variable names
non_DNAm <- c('newid', 'Chrono_Age', 'Timepoint', 'Array', 'Method')

# figure out how to select DNAm age estimate to use based on following hierarchy
# horvath CpG: 450, EPICv2, EPICv1
# all other CpG clocks & all PC clocks: 450, EPICv1, EPICv2
hierarchy_horvath <- dnam_age %>%
  # make data longer for easier data manipulation
  pivot_longer(cols = -any_of(non_DNAm),
               names_to = 'clock',
               values_to = 'DNAm_age') %>%
  # specify if PC-based or CpG-based estimates
  mutate(type = case_when(
    Method == 'Morgan Levine PC' | Method == 'Methscore PC' ~ 'PC',
    .default = 'CpG')) %>%
  # specify heirarchy for everything other than CpG-based horvath estimates
  filter(type == 'CpG' & clock == 'Horvath') %>%
  # create variable that indicates the use priority of each clock/array combo
  mutate(priority = case_when(
    # for all clocks, prioritize using data from 450K array first
    Array == '450K' ~ 1,
    Array == 'EPICv2' ~ 2,
    Array == 'EPICv1' ~ 3))

# select estimates to keep using hierarchy given above for horvath CpG-based
hierarchy_else <- dnam_age %>%
  # make data longer for easier data manipulation
  pivot_longer(cols = -any_of(non_DNAm),
               names_to = 'clock',
               values_to = 'DNAm_age') %>%
  # specify if PC-based or CpG-based estimates
  mutate(type = case_when(
    Method == 'Morgan Levine PC' | Method == 'Methscore PC' ~ 'PC',
    .default = 'CpG')) %>%
  # specify heirarchy for everything other than CpG-based horvath estimates
  filter(!(type == 'CpG' & clock == 'Horvath')) %>%
  # create variable that indicates the use priority of each clock/array combo
  mutate(priority = case_when(
    # for all clocks, prioritize using data from 450K array first
    Array == '450K' ~ 1,
    Array == 'EPICv1' ~ 2,
    Array == 'EPICv2' ~ 3))

# select estimates to keep using hierarchy given above for all estimates other
# than horvath CpG-based estimates
hierarchy_all <- rbind(hierarchy_else, hierarchy_horvath) %>%
  # keep highest priority DNAm age available for each participant
  group_by(newid, Timepoint, Method, clock) %>%
  # keep highest priority estimate for each group (priority == 1 = highest)
  slice_min(priority) %>%
  ungroup()
# HIGHEST PRIORITY IS NOT THE ONLY ONE KEPT?????

View(filter(hierarchy_all, newid == '1004'))

# check that everyone has 1 estimate per timepoint/method/clock
hierarchy_all %>%
  # identify participants that have more than 1 measure per timepoint
  group_by(newid, Timepoint, Method, clock) %>%
  summarise(measure_count = n_distinct(Chrono_Age)) %>%
  View()
# everyone has only 1 estimate per timepoint/method/clock

# check that all estimates for participants were retained
# original data 
original <- dnam_age %>%
  # make data longer for easier data manipulation
  pivot_longer(cols = -any_of(non_DNAm),
               names_to = 'clock',
               values_to = 'DNAm_age') %>%
  select(newid, Timepoint, Method, clock) %>%
  distinct() %>%
  # arrange all columns identically for comparison
  arrange(across(everything(), desc))

# prioritized data
prioritized <- hierarchy_all %>%
  select(newid, Timepoint, Method, clock) %>%
  # arrange all columns identically for comparison
  arrange(across(everything(), desc))
  
# show that dataframes are the same
arsenal::comparedf(original, prioritized) %>% summary(.)
# dataframes are exactly the same, all estimates were retained for participants

# prepare dataframe format for use with gee ------------------------------------
# factor variables
asthma$LCA <- factor(asthma$LCA,
                            levels = c('never/infrequent', # reference: N/I
                                       'late onset', 
                                       'persistent'))

# merge asthma trajectories with CpG-based DNAm age estimated with Methscore
dnam_asthma <- hierarchy_all %>% # n=599
  select(-priority) %>%
  # make data wider
  pivot_wider(names_from = clock,
              values_from = DNAm_age) %>%
  filter(Timepoint %in% timepoints) %>% # n=453
  # change variable type for merging
  mutate(newid = as.integer(newid)) %>%
  # merge asthma class to estimated epi ages
  left_join(., asthma, by = 'newid') %>%
  # remove participants missing main outcome: asthma trajectory class
  filter(!is.na(LCA)) # n=452

# fit marginal model with GEE --------------------------------------------------

# should i be using exchangeable + chronological age: less efficient than AR?

# or autoregressive + timepoint (rounded age): will result in potentially larger
# epigenetic age deviation than if i just used accurate chronological age. Can
# only use timepoint here because AR structure assumes equally spaced time points

# function that fits marginal model using GEE
fit_model <- function(clock, # name of clock
                      method, # method used to estimate DNAm age 
                      df # dataframe
                      ){

  # subset data to isolate given clock/method/type combination
  subset_df <- filter(df, Method == method)
  
  # create formula
  form <- reformulate(c('Chrono_Age', 'LCA'), response = clock)
  
  # fit marginal model using GEE with exchangeable working correlation
  gee(formula = form,
               id = newid,
               data = subset_df,
               family = gaussian,
               corstr = 'exchangeable')
  
}

# function that creates 95% CI using robust SE
robust_ci <- function(model # gee object
                      ){
  
  # extract clock name (LHS) from terms
  clock_name <- as.character(model$terms[[2]])
  
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
    clock = clock_name
  )
  
  return(ci_table)
  
}


methscore_cpg <- lapply(c('Horvath', 'SkinBlood', 'Hannum', 'PhenoAge'),
                        fit_model, 
                        method = 'Methscore CpG', 
                        df = dnam_asthma
                        ) %>%
  map(robust_ci) %>%
  list_rbind() %>%
  mutate(type = 'CpG')
  

methscore_pc <- lapply(c('Horvath', 'SkinBlood', 'Hannum', 'PhenoAge', 'GrimAge'),
                       fit_model, 
                       method = 'Methscore PC', 
                       df = dnam_asthma
                       ) %>%
  map(robust_ci) %>%
  list_rbind() %>%
  mutate(type = 'PC')

clockfound_cpg <- lapply(c('GrimAge', 'IEAA'),
                         fit_model, 
                         method = 'Clock Foundation', 
                         df = dnam_asthma
                         ) %>%
  map(robust_ci) %>%
  list_rbind() %>%
  mutate(type = 'CpG')

results <- do.call(rbind, list(methscore_cpg, methscore_pc, clockfound_cpg))

# output -----------------------------------------------------------------------
write.csv(results, 'data-processed/gee-results.csv', row.names = TRUE)