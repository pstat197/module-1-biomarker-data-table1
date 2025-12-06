library(tidyverse)
# read in the raw data 
raw_data <- read_csv('data/biomarker-raw.csv', 
                     col_names = F, 
                     skip = 2, # names + abbreviations
                     col_select = -(1:2), # drops first 2 columns
                     na = c('-', '')) # tells to treat these values as missing

# output data makes each row a subject, each column a protein
## and values are raw concentration levels 

set.seed(123) # make sample reproducible

view(raw_data)

# pick 6 protein columns randomly
sample_proteins <- sample(names(raw_data), 6)
sample_proteins

raw_data %>% 
  select(all_of(sample_proteins)) %>% 
  pivot_longer( # reshapes data to long format
    everything(),
    names_to = 'protein',
    values_to = 'value'
  ) %>% 
  ggplot(aes(x=value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ protein, scales = 'free') + # creates separate mini-plots
  ## for each protein (~protein means one facet per unique values of protein)
  ## scales = free means x-axis varies by facet
  labs(title = 'raw protein dist')

## these plots from the sample proteins
## are very right-skewed which means we may want to consider
## a log transformation

log_data <- raw_data %>% 
  select(all_of(sample_proteins)) %>% 
  mutate(across(everything(), log10)) %>% 
  pivot_longer(
    everything(), 
    names_to = 'protein',
    values_to = 'value'
  ) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ protein, scales = 'free') + 
  labs(title = 'log transformed protein dist')

log_data
# output data is a lot more symmetric/normal and less 
## extreme, proving why log transformation is useful here