library(tidyverse)

# Question 01: Why log?

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

# Question 02: Why trim?

# read in numeric protein data
raw_numeric <- read_csv(
  "data/biomarker-raw.csv",
  skip = 1,# skip the name row
  na = c("-", "")
)

names(raw_numeric)[1:5]
tail(names(raw_numeric))

biomarker_raw <- raw_numeric %>% 
  rename(
    group = '...1',
    ados = '...1320'
  ) %>% 
  select(-Target) %>% # drops empty target column
  filter(!is.na(group))

dim(biomarker_raw)
head(biomarker_raw[, 1:6])
# want z-scores with values w mean 0, sd 1

biomarker_z <- biomarker_raw %>% 
  mutate(across(
    .cols = -c(group, ados), 
    .fns = ~ scale(log10(.x))[,1] # log10 and standardize
  ))

# log transform reduces how skewed raw data is 
# scale converts each protein to z-scores
# basically preprocessing script before trim()

# flag outliers
outlier_counts <- biomarker_z %>% 
  mutate(across(
    .cols = -c(group, ados),
    .fns = ~ abs(.x) > 3 # true if val is outlier
  )) %>% 
  mutate(
    n_outliers = rowSums(across(-c(group, ados)))
  ) %>% select(group, ados, n_outliers)

# rowwise counting because we want to know which subjects
# have extreme values and whether that differs between groups

# if most subjects have just a few outliers but a few have
## over 100, those subjects may be important and change the models

summary(outlier_counts$n_outliers)

outlier_counts %>% 
  ggplot(aes(x = n_outliers)) +
  geom_histogram(bins = 30) +
  labs(
    title = 'number of outlier protein values per subject (|z| > 3)',
    x = 'number of outlying protein values',
    y = 'number of subjects'
  )

# compare asd vs td
group_outliers <- outlier_counts %>%
  group_by(group) %>%
  summarise(
    mean_outliers   = mean(n_outliers),
    median_outliers = median(n_outliers),
    max_outliers    = max(n_outliers),
    .groups = "drop"
  )

group_outliers

top_outlier_subjects <- outlier_counts %>% 
  arrange(desc(n_outliers)) %>% 
  slice(1:10)

top_outlier_subjects

# removing the trimming step and examining the outliers defined as |z| > 3 after 
# log transformation and standardization, i found that most subjects had 
# relatively few protein outliers. however, a small subset of subjects showed 
# extremely large numbers of outliers across many different proteins. 
# specifically, they are outliers at the subject level and not just individual
# protein-levels. 

# the comparison shows that ASD subjects have 13.2 mean outliers and 126 
# maximum outliers. the TD subjects have 17.6 mean outliers and 
# maximum 157.

# both groups include subjects with very high outlier counts, but TD subjects 
# show a slightly higher mean and maximum number of outliers in the 
# untrimmed dataset.

# overall, after removing trimming, it seems that although most subjects had
# few outliers, a small number had over 120-150 extreme values across proteins.
# these were mostly TD participants and were clear subject-level outlines because
# such cases will disproportionately affect variable selection and classification
# so trimming is essential to stabilize