# Author: Selene Banuelos
# Date: 7/15/2026
# Description: Try to separate out CHAM1 from CHAM2 participants

# Notes: I know that CHAM2 recruitment started at time point 9Y. Following this, 
# I can assume that there is a group of participant IDs is only present between 
# time points 9Y-18Y (CHAM2) and another group of IDs that are present across
# all times (CHAM1)

# setup
library(readstata13)
library(dplyr)

# import data ------------------------------------------------------------------
# survey data
asthma <- read.dta13('data-raw/de_la_Rosa_06.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
clean <- asthma %>%
  # focus on asthma variables 
  select(newid, contains('asth_')) %>%
  # remove rows that have NA in all columns other than newid column
  filter(!if_all(-newid, is.na))

# separate time points into two groups
# 5Y-7Y
five_7 <- clean %>%
  # keep vars only corresponding to time points of interest
  select(newid, contains(c('_5Y', '_7Y'))) %>%
  # remove rows that have NA in all columns other than newid column
  filter(!if_all(-newid, is.na))

# 9Y-18Y
nine_18 <- clean %>%
  # keep vars only corresponding to time points of interest
  select(newid, contains(c('_9Y', '_10Y', '_12Y', '_14Y', '_16Y', '_18Y'))) %>%
  # remove rows that have NA in all columns other than newid column
  filter(!if_all(-newid, is.na))

# which IDs are in 9Y-18Y but missing from 5Y/7Y
# I assume these are CHAM2 participants
cham2 <- setdiff(nine_18$newid, five_7$newid)

# the remaining participants with asthma data collected at at least 1 time point
# are assumed to be CHAM1
cham1 <- clean %>%
  filter(!newid %in% cham2)
# don't know how I feel about this... think this through more tomorrow

  