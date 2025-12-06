library(tidyverse)

# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                     col_names = F, 
                     n_max = 2, 
                     col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

'interpretation: reads the first two rows of biomarker-raw.csv (excludes first 
2 columns), row 1 = full protein names, row 2 = abbreviations, transposes t()
so each protein is a row with name and abbreviation' 


# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# this is winsorization: any value with |x| > 3 gets capped

# read in data
biomarker_clean <- read_csv('data/biomarker-raw.csv', 
         skip = 2,
         col_select = -2L,
         col_names = c('group', 
                       'empty',
                       pull(var_names, abbreviation),
                       'ados'),
         na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())

'interpretation: skips the first 2 header rows and reads the actual numeric
data, drops column 2 (col_select = -2L), sets column names: group, empty, all
protein abbreviations, and ados

filters out rows where group is missing, for every protein column: log10()
transformation, scale() (mean 0, SD 1), trim values with |z| > 3, reorders
columns to put group and ados first, and saves the processed tibble'

# export as r binary
save(list = 'biomarker_clean', 
     file = 'data/biomarker-clean.RData')


