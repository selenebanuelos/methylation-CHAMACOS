# Author: Selene Banuelos
# Date: 6/30/2026
# Description: Compare correlation and median absolute error between 
# chronological age and DNAm age estimated with several methods

# setup
library(readstata13) # work with STATA DTA files
library(dplyr)
library(janitor)
library(tidyr)
library(Metrics) # calculate median absolute error
library(gtsummary)

# import data ------------------------------------------------------------------
# estimated epigenetic age using various methods
dnam_age <- read.dta13("data-raw/de_la_Rosa_epigenetic_06.dta", 
                       nonint.factors=TRUE, 
                       generate.factors=TRUE)

# data wrangling ---------------------------------------------------------------
# create vector of methods we want to compare
methods <- c('methylCIPHER', 
             'Morgan Levine PC', 
             'DunedinPACE',
             'Methscore CpG', 
             'Methscore PC')

# create vector of timepoints of interest
timepoints <- c('Age 9', 'Age 12', 'Age 14', 'Age 18')

# reformat data frame
reformat <- dnam_age %>%
  # keep estimates from methods of interest
  filter(Method %in% methods,
         # keep estimates from tiempoints of interest
         Timepoint %in% timepoints) %>%
  # remove empty columns
  remove_empty('cols') %>%
  # make data longer to calculate correlations and MAE
  pivot_longer(cols = c('Horvath',
                        'SkinBlood',
                        'Hannum',
                        'PhenoAge',
                        'GrimAge',
                        'DNAmTL',
                        'DunedinPACE'),
               names_to = 'clock',
               values_to = 'dnam_age'
               )

# correlation and MAE ----------------------------------------------------------
# wrapper function to use stats::cor.test() within mutate
calculate_corr <- function(x, # numeric vector
                           y # numeric vector
                           ){
  
  # Pearson's product-moment correlation
  result <- cor.test(x, y, method = c('pearson'))
  
  # get correlation coefficient, r
  r <- result$estimate %>% round(digits = 2)
  
  # get 95% CI for r
  ci <- result$conf.int %>% round(digits = 2)
  
  # return string containing: "r,95% CI lower bound,95% CI upper bound"
  paste(c(r, ci), collapse = ',')
  
}

# calculate these within each timepoint, method, and clock combination
corr_mae <- reformat %>%
  # calculate correlation coefficient within each combination
  group_by(Timepoint, Method, clock) %>%
  # calculate cross-sectional Pearson correlation coefficient
  mutate(pearson = calculate_corr(.$Chrono_Age, .$dnam_age)) %>%
  # separate Pearson's r and 95% CI, split by ',' delimiter
  separate_wider_delim(cols = pearson, 
                       delim = ',',
                       names = c('pearson_r', 'ci_low', 'ci_high')) %>%
  # convert r and 95% CI into numeric
  mutate(across(.cols = pearson_r:ci_high, 
                .fns = as.numeric)) %>%
  # calculate cross-sectional median absolute error
  mutate(mae = round(mdae(actual = Chrono_Age, predicted = dnam_age),
                     digits = 2)
         )%>%
  ungroup(.) %>%
  # create column with complete 95% CI
  mutate(ci_95 = paste0('(', ci_low, ', ', ci_high, ')')) %>%
  # remove separate 95% CI bounds columns
  select(-c(ci_low, ci_high)) %>%
  # make data wider to move clock names into column names
  pivot_wider(names_from = clock,
              names_sep = '_',
              values_from = c('pearson_r', 'mae', 'ci_95')) %>%
  # make data longer to move stats from names into one column
  pivot_longer(cols = contains('pearson_r_'),
               names_to = 'statistic'
               )

# create tables ----------------------------------------------------------------
corr_mae %>%
  tbl_summary(include = c(pearson_r, mae))