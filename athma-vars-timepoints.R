# Author: Selene Banuelos
# Date: 7/15/2026
# Description: Compare asthma-related variables at different time points

# setup
library(readstata13)
library(dplyr)
library(tidyr)
library(gtsummary)
library(gt)
library(ggplot2)

# import data ------------------------------------------------------------------
# asthma variables
asthma <- read.dta13('data-raw/de_la_Rosa_06.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# unconfirmed cohort membership
cohorts <- read.csv('data-processed/cohorts-unconfirmed.csv',
                   colClasses = c(newid = 'character'))

# data wrangling ---------------------------------------------------------------
# add in cohort information to asthma dataset
asthma_cohort <- left_join(asthma, cohorts, by = 'newid')

# reformat data to use with gtsummary::tbl_summary()
long <- asthma_cohort %>%
  # make table with asthma related variables only
  select(newid, cohort, contains('asth_')) %>%
  # make data longer
  pivot_longer(cols = contains('asth_'),
               names_to = c('variable', 'age'),
               names_pattern = 'asth_(.+)_(\\d+Y)$')

# spread data wider
wide <- long %>%
  pivot_wider(id_cols = c(newid, cohort, age),
              names_from = variable,
              values_from = value)

# factor time point variable to control order displayed in table
long$age <- factor(long$age,
                            levels = c('5Y','7Y','9Y','10Y','12Y','14Y','16Y','18Y')
)
wide$age <- factor(wide$age,
                   levels = c('5Y','7Y','9Y','10Y','12Y','14Y','16Y','18Y')
)

# table ------------------------------------------------------------------------
# non-stratified table
wide %>%
  filter(cohort != 'no data') %>%
  tbl_summary(include = c(sym, symed, med, diag_ever), 
              by = age,
              missing_text = 'Missing',
              # show missing as n (%)
              missing_stat = '{N_miss} ({p_miss}%)')


# table stratified by (unconfirmed) CHAMACOS cohort
wide %>%
  filter(cohort != 'no data') %>%
  tbl_strata(
    strata = cohort,
    .tbl_fun = 
      ~ .x %>%
      tbl_summary(include = c(sym, symed, med, diag_ever), 
                  by = age,
                  missing_text = 'Missing',
                  # show missing as n (%)
                  missing_stat = '{N_miss} ({p_miss}%)'),
    .header = '**{strata}**'
  )
# looks like there's different rates of missingness across vars at 9Y, 12Y, 14Y,
# 16Y, 18Y

# table without data on missingness
wide %>%
  tbl_summary(include = c(sym, med, diag_ever), 
              by = age,
              label = list(sym = 'asth_sym_*y',
                           med = 'asth_med_*y',
                           diag_ever = 'asth_diag_ever_*y'),
              # show yes and no n(%)
              type = all_dichotomous() ~ 'categorical',
              missing = 'no') %>%
  # change header to 'Variable'
  modify_header(label = '**Variable**') %>%
  modify_spanning_header(all_stat_cols() ~ "**Overall**")

# table without data on missingness, stratified by (unconfirmed) CHAMACOS cohort
wide %>%
  filter(cohort != 'no data') %>%
  tbl_strata(
    strata = cohort,
    .tbl_fun = 
      ~ .x %>%
      tbl_summary(include = c(sym, med, diag_ever), 
                  by = age,
                  label = list(sym = 'asth_sym_*y',
                               med = 'asth_med_*y',
                               diag_ever = 'asth_diag_ever_*y'),
                  # show yes and no n(%)
                  type = all_dichotomous() ~ 'categorical',
                  missing = 'no'),
    .header = '**{strata}**'
  ) %>%
  # change header to 'Variable'
  modify_header(label = '**Variable**')

# spaghetti plots --------------------------------------------------------------
# plot 'yes' responses over time for each of 3 vars, stratified by cohort
long %>%
  # remove participants with no data across 5Y-18Y
  filter(cohort != 'no data') %>%
  group_by(cohort, variable, age) %>%
  summarise(yes_count = sum(value == 'Yes', na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = age, y = yes_count, group = variable)) +
  geom_line(aes(color = variable)) +
  facet_wrap(vars(cohort)) +
  theme_minimal()
# there are no 'Yes' responses to med, sym, symed at 9Y in CHAM2, which is weird

# plot 'no' responses over time for each of 3 vars, stratified by cohort
long %>%
  # remove participants with no data across 5Y-18Y
  filter(cohort != 'no data') %>%
  group_by(cohort, variable, age) %>%
  summarise(no_count = sum(value == 'No', na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = age, y = no_count, group = variable)) +
  geom_line(aes(color = variable)) +
  facet_wrap(vars(cohort)) +
  theme_minimal()
# everyone responded 'No' to med, sym, symed at 9Y in CHAM2

# plot 'missing' responses over time for each  # remove participants with no data across 5Y-18Y
long %>%
  # remove participants with no data across 5Y-18Y
  filter(cohort != 'no data') %>%
  group_by(cohort, variable, age) %>% 
  summarise(missing_count = sum(is.na(value))) %>%
  ungroup() %>%
  ggplot(aes(x = age, y = missing_count, group = variable)) +
  geom_line(aes(color = variable)) +
  ylim(0, 65) +
  facet_wrap(vars(cohort)) +
  theme_minimal()

# output -----------------------------------------------------------------------
# save tables
#gtsave(no_missing, 'figures/asthma-vars-overall.png')