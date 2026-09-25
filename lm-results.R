# Author: Selene Banuelos
# Date: 9/24/2026
# Description: Fit linear model to compare unadjusted average DNAm age between
# asthma trajectory classes at 18Y visit. 

# DNAm age ~ chrono age + asthma traj class

# setup
library(readstata13)
library(dplyr)
library(sandwich)

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

# fit linear model & calculate 95% CI ------------------------------------------
# function that fits linear model
fit_model <- function(clock, # name of clock
                      method, # method used to estimate DNAm age 
                      df # dataframe
){
  
  # subset data to isolate given clock/method/type combination
  subset_df <- filter(df, Method == method)
  
  # create formula ('clock' ~ Chrono_Age + LCA)
  form <- reformulate(c('Chrono_Age', 'LCA'), response = clock)
  
  # fit marginal model using GEE with exchangeable working correlation
  lm(formula = form, data = subset_df)
  
}

# Huber-White robust 'sandwich' estimator, per FDA recommendations
# example, using CpG-based Horvath clock estimated using methscore
# robust CI is typically wider than classical

# function that creates 95% CI using robust SE (sandwich estimator)
robust_ci <- function(model){
  
  # extract clock name (LHS) from terms
  clock_name <- as.character(model$terms[[2]])
  
  # estimate robust covariance matrix with sandwich estimator
  hat_cov <- sandwich::vcovHC(model, 
                    type = 'HC3' # recommended in Long & Ervin (2000) for lm
                    )

  # extract estimated variances for all terms
  hat_var <- diag(hat_cov)
  
  # calculate SEs for all terms
  hat_se <- sqrt(hat_var)
  
  # extract coefficients
  hat_beta <- model$coefficients
  
  # calculate lower and upper bounds of 95% CIs
  lower <- hat_beta - 2*hat_se
  upper <- hat_beta + 2*hat_se
  
  # create string 'beta(lower, upper)' for easy copy & paste if needed
  beta_ci <- paste0(round(hat_beta, digits = 2),
                    ' (', 
                    round(lower, digits = 2), 
                    ', ', 
                    round(upper, digits = 2), 
                    ')')
  
  # build 95% CIs
  ci_table <- data.frame(
    term = names(hat_beta),
    beta_95_ci = beta_ci, 
    clock = clock_name,
    beta = hat_beta,
    lower_95_ci = lower,
    upper_95_ci = upper,
    row.names = NULL
  )
}

# get results for CpG-based DNAm estimated using Methscore
methscore_cpg <- lapply(c('Horvath', 'SkinBlood', 'Hannum', 'PhenoAge'),
                        fit_model, 
                        method = 'Methscore CpG', 
                        df = dnam_asthma
                        ) %>%
  map(robust_ci) %>%
  list_rbind() %>%
  mutate(type = 'CpG')

# get results for PC-based DNAm estimated using Methscore
methscore_pc <- lapply(c('Horvath', 'SkinBlood', 'Hannum', 'PhenoAge', 'GrimAge'),
                       fit_model, 
                       method = 'Methscore PC', 
                       df = dnam_asthma
                       ) %>%
  map(robust_ci) %>%
  list_rbind() %>%
  mutate(type = 'PC')

# get results for CpG-based DNAm estimated using Clock Foundation
clockfound_cpg <- lapply(c('GrimAge', 'IEAA'),
                         fit_model, 
                         method = 'Clock Foundation', 
                         df = dnam_asthma
                         ) %>%
  map(robust_ci) %>%
  list_rbind() %>%
  mutate(type = 'CpG')

# combine all results together
results <- do.call(rbind, list(methscore_cpg, methscore_pc, clockfound_cpg))

# output -----------------------------------------------------------------------
write.csv(results, 'data-processed/lm-results.csv', row.names = FALSE)