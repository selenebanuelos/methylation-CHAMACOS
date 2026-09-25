# Author: Selene Banuelos
# Date: 9/24/2026
# Description: Fit linear model to compare unadjusted average DNAm age between
# asthma trajectory classes at 18Y visit. 

# DNAm age ~ chrono age + asthma traj class

# setup
library(readstata13)
library(dplyr)

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

# factor variables
asthma$LCA <- factor(asthma$LCA,
                     levels = c('never/infrequent', # reference: N/I
                                'late onset', 
                                'persistent'))

# merge asthma trajectories with CpG-based DNAm age estimated with Methscore
dnam_asthma <- dnam_age %>% # n=599
  # cross-sectional analysis using 18Y visit data
  filter(Timepoint == 'Age 18') %>% # n=378
  # change variable type for merging
  mutate(newid = as.integer(newid)) %>%
  # merge asthma class to estimated epi ages
  left_join(., asthma, by = 'newid') %>%
  # remove participants missing main outcome: asthma trajectory class
  filter(!is.na(LCA)) # n=378

# fit linear model -------------------------------------------------------------
# function that fits linear model
fit_model <- function(clock, # name of clock
                      method, # method used to estimate DNAm age 
                      df # dataframe
){
  
  # subset data to isolate given clock/method/type combination
  subset_df <- filter(df, Method == method)
  
  # create formula
  form <- reformulate(c('Chrono_Age', 'LCA'), response = clock)
  
  # fit marginal model using GEE with exchangeable working correlation
  lm(formula = form, data = subset_df)
  
}

methscore_cpg <- lapply(c('Horvath', 'SkinBlood', 'Hannum', 'PhenoAge'),
                        fit_model, 
                        method = 'Methscore CpG', 
                        df = dnam_asthma
                        ) 

methscore_pc <- lapply(c('Horvath', 'SkinBlood', 'Hannum', 'PhenoAge', 'GrimAge'),
                       fit_model, 
                       method = 'Methscore PC', 
                       df = dnam_asthma
                       )

clockfound_cpg <- lapply(c('GrimAge', 'IEAA'),
                         fit_model, 
                         method = 'Clock Foundation', 
                         df = dnam_asthma
                         )

  map(robust_ci) %>%
  list_rbind() %>%
  mutate(type = 'CpG')
  


# use sandwich estimator to get robust SE