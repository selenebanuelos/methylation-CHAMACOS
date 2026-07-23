# Author: Selene Banuelos
# Date: 7/20/2026
# Description: Latent class analysis to estimate asthma trajectories

# setup
library(readstata13)
library(dplyr)
library(poLCA) # previously used for LCA
library(purrr)
library(stringr)
library(ggplot2)
options(scipen = 999)

# import data ------------------------------------------------------------------
# asthma-related data
asthma <- read.dta13('data-raw/de_la_Rosa_07.dta',
                     nonint.factors = TRUE,
                     generate.factors = TRUE)

# data wrangling ---------------------------------------------------------------
# do some reformatting to prepare for LCA
current_asthma <- asthma %>%
  # keep asthma related variables only
  dplyr::select(newid, cham, contains('asth_')) %>%
  # remove participants that are missing all data between 9Y-18Y
  filter(!if_all(contains(c('9y', '10Y', '12Y', '14Y', '18Y')), is.na)) %>%
  # make data longer for some data manipulation
  pivot_longer(cols = contains('asth_'),
               names_to = c('.value', 'age'),
               names_pattern = '(asth_.*)_([0-9]+Y)$') %>%
  # remove data from ages 5 and 7
  filter(age != '5Y' & age != '7Y') %>%
  # remove data from participants with only 1 visit (can't create trajectory?)
  group_by(newid) %>%
  filter(n_distinct(age) > 1) %>%
  ungroup() %>%
  # create a 'current asthma' variable for each timepoint
  # current asthma defined as 2/3 of the following: current asthma symptoms, 
  # current asthma medication, ever asthma diagnosis
  # asth_symed = current asthma symptoms OR current asthma medication use
  mutate(
    current_asthma = case_when(
      asth_symed == 'Yes' & asth_diag_ever == 'Yes' ~ 1, # yes
      asth_symed == 'No' | asth_diag_ever == 'No' ~ 0, # no
      is.na(asth_symed) | is.na(asth_diag_ever) ~ NA),
    # convert to factor to use in LCA
    current_asthma = factor(current_asthma,
                            levels = c(1, 0),
                            labels = c('Yes', 'No')
                            )) %>%
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
              names_from = age,
              values_from = c(asth_sym, 
                              asth_symed, 
                              asth_med, 
                              asth_diag_ever, 
                              current_asthma),
              names_glue = '{.value}_{age}')

# conditional independence assumption: latent class membership explains all of
# the shared variance among the observed indicators

# correlation of candidate indicators ------------------------------------------
# correlation between current asthma at different ages
asthma %>%
  # keep asthma related variables only
  dplyr::select(newid, cham, contains('asth_')) %>%
  # remove participants that are missing all data between 9Y-18Y
  filter(!if_all(contains(c('9y', '10Y', '12Y', '14Y', '18Y')), is.na)) %>%
  # make data longer for some data manipulation
  pivot_longer(cols = contains('asth_'),
               names_to = c('.value', 'age'),
               names_pattern = '(asth_.*)_([0-9]+Y)$') %>%
  # remove data from ages 5 and 7
  filter(age != '5Y' & age != '7Y') %>%
  # remove data from participants with only 1 visit (can't create trajectory?)
  group_by(newid) %>%
  filter(n_distinct(age) > 1) %>%
  ungroup() %>%
  # create a 'current asthma' variable for each timepoint
  # current asthma defined as 2/3 of the following: current asthma symptoms, 
  # current asthma medication, ever asthma diagnosis
  # asth_symed = current asthma symptoms OR current asthma medication use
  mutate(
    current_asthma = case_when(
      asth_symed == 'Yes' & asth_diag_ever == 'Yes' ~ 1, # yes
      asth_symed == 'No' | asth_diag_ever == 'No' ~ 0, # no
      is.na(asth_symed) | is.na(asth_diag_ever) ~ NA)) %>%
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
              names_from = age,
              values_from = c(asth_sym, 
                              asth_symed, 
                              asth_med, 
                              asth_diag_ever, 
                              current_asthma),
              names_glue = '{.value}_{age}') %>%
  dplyr::select(contains('current_asthma')) %>%
  # calculate tetrachoric correlations
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho
# all correlations are >0.5

# correlation between current asthma sym, current med use, ever asthma diagnosis
asthma %>%
  # keep asthma related variables only
  dplyr::select(newid, cham, contains('asth_')) %>%
  # remove participants that are missing all data between 9Y-18Y
  filter(!if_all(contains(c('9y', '10Y', '12Y', '14Y', '18Y')), is.na)) %>%
  # make data longer for some data manipulation
  pivot_longer(cols = contains('asth_'),
               names_to = c('.value', 'age'),
               names_pattern = '(asth_.*)_([0-9]+Y)$') %>%
  # remove data from ages 5 and 7
  filter(age != '5Y' & age != '7Y') %>%
  # remove data from participants with only 1 visit (can't create trajectory?)
  group_by(newid) %>%
  filter(n_distinct(age) > 1) %>%
  ungroup() %>%
  # make data wider to use with poLCA::poLCA()
  pivot_wider(id_cols = c('newid', 'cham'),
            names_from = age,
            values_from = c(asth_sym, 
                            asth_symed, 
                            asth_med, 
                            asth_diag_ever),
            names_glue = '{.value}_{age}') %>%
  dplyr::select(contains(c('asth_sym_', 'asth_med_', 'asth_diag_'))) %>%
  # change all 'Yes' to 1, 'No to 0
  mutate(across(everything(), ~ case_when(. == 'Yes' ~ 1,
                                          . == 'No' ~ 0,
                                          is.na(.) ~ NA))) %>%
  # calculate tetrachoric correlations
  sirt::tetrachoric2(.) %>%
  # view correlation matrix
  .$rho

# check data missingness among analytic sample ---------------------------------


# Rosie's code -----------------------------------------------------------------
##Model function
function_9_18Y <- as.formula(
  paste("cbind(", 
        paste(c("current_asthma_9Y",
                "current_asthma_10Y",
                "current_asthma_12Y",
                "current_asthma_14Y",
                "current_asthma_16Y",
                "current_asthma_18Y"), 
              collapse = ","),
        ") ~ 1")
  )
# so, we're not using covariates to stimate latent class membership?
# is ever asthma variable correct?

#Run models for 2 to 5 classes
set.seed(1234)
lca_9_18Y <- lapply(2:5, function(k) {
  message("Currently estimating model with ", k, " classes...")
  poLCA(function_9_18Y ,
        data = current_asthma,
        nclass = k,
        nrep = 100,
        maxiter = 5000,
        na.rm = FALSE,
        verbose = FALSE)
})

#Create the Comparison Table of Models
lca_2_timepoint <- data.frame(
  Classes = 2:5,
  Log_Likelihood = sapply(lca_9_18Y, function(m) round(m$llik, 2)),
  BIC = sapply(lca_9_18Y, function(m) round(m$bic, 2)),
  AIC = sapply(lca_9_18Y, function(m) round(m$aic, 2)),
  Smallest_Class_Pct = sapply(lca_9_18Y, function(m) {
    round(min(m$P) * 100, 1)
  })
)

#Output results
print(lca_2_timepoint)

set.seed(1234)
two_class <- poLCA(function_9_18Y,
                   data = current_asthma,
                   nclass = 2,
                   nrep = 100,
                   maxiter = 5000,
                   na.rm = FALSE,
                   verbose = FALSE)

set.seed(1234)
three_class <- poLCA(function_9_18Y,
                    data = current_asthma,
                    nclass = 3,
                    nrep = 100,
                    maxiter = 5000,
                    na.rm = FALSE,
                    verbose = FALSE)

set.seed(1234)
four_class <- poLCA(function_9_18Y,
                     data = current_asthma,
                     nclass = 4,
                     nrep = 100,
                     maxiter = 5000,
                     na.rm = FALSE,
                     verbose = FALSE)

set.seed(1234)
five_class <- poLCA(function_9_18Y,
                    data = current_asthma,
                    nclass = 5,
                    nrep = 100,
                    maxiter = 5000,
                    na.rm = FALSE,
                    verbose = FALSE)

# check class sizes - don't want small classes
# check entropy - higher entropy indicates better class separation
table(two_class$predclass)
poLCA.entropy(two_class)

table(three_class$predclass)
poLCA.entropy(three_class)

table(four_class$predclass)
poLCA.entropy(four_class)

table(five_class$predclass)
poLCA.entropy(five_class)

# Visualize trajectories -------------------------------------------------------
get_probs <- function(age_list) {

  as.data.frame(age_list) %>%
    mutate(class = str_extract(row.names(.), '[0-9]+'))
}

map_df(three_class$probs, get_probs, .id = 'age') %>%
  mutate(age = str_extract(age, '[0-9]+'),
         age = factor(age,levels = c('9', '10', '12', '14', '16', '18'))) %>%
  ggplot(aes(x = age, 
             y = `Pr(1)`,
             color = class)) +
  geom_line(aes(group = class))

map_df(four_class$probs, get_probs, .id = 'age')  %>%
  mutate(age = str_extract(age, '[0-9]+'),
         age = factor(age,levels = c('9', '10', '12', '14', '16', '18'))) %>%
  ggplot(aes(x = age, 
             y = `Pr(1)`,
             color = class)) +
  geom_line(aes(group = class))

map_df(five_class$probs, get_probs, .id = 'age')  %>%
  mutate(age = str_extract(age, '[0-9]+'),
         age = factor(age,levels = c('9', '10', '12', '14', '16', '18'))) %>%
  ggplot(aes(x = age, 
             y = `Pr(1)`,
             color = class)) +
  geom_line(aes(group = class))

