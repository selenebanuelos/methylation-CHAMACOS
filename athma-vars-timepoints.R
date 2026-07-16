# Author: Selene Banuelos
# Date: 7/15/2026
# Description: Compare asthma-related variables at different time points

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(gtsummary)
library(ggplot2)

# import data ------------------------------------------------------------------
# asthma variables
asthma <- read.dta13('data-raw/de_la_Rosa_06.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
# reformat data to use with gtsummary::tbl_summary()
long <- asthma %>%
  # make table with asthma related variables only
  select(newid, contains('asth_')) %>%
  # make data longer
  pivot_longer(cols = contains('asth_'),
               names_to = c('variable', 'age'),
               names_pattern = 'asth_(.+)_(\\d+Y)$')

# spread data wider
wide <- long %>%
  pivot_wider(id_cols = c(newid, age),
              names_from = variable,
              values_from = value)

# factor time point variable to control order displayed in table
wide$age <- factor(wide$age,
                   levels = c('5Y','7Y','9Y','10Y','12Y','14Y','16Y','18Y')
                   )

# data visualization -----------------------------------------------------------
# table
wide %>%
  tbl_summary(include = c(sym, symed, med, diag_ever), 
              by = age,
              missing_text = 'Missing',
              # show missing as n (%)
              missing_stat = '{N_miss} ({p_miss}%)')
# looks like there's different rates of missingness across vars at 9Y, 12Y, 14Y,
# 16Y, 18Y

# use plots below to visualize this info
# spaghetti plots
# plot 'yes' responses over time for each of 3 vars
long %>%
  ggplot(aes(x = age, y = )


# plot 'missing' responses over time for each of 3 vars
