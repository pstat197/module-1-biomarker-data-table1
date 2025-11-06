# Katie's Script
library(tidyverse)
library(ggplot2)

set.seed(101422)

# QUESTION 1
# loads preprocessed data
load('data/biomarker-clean.RData')
#/Users/katiepyo/Desktop/PSTAT/PSTAT197/module-1-biomarker-data-table1"

# check what objects we loaded
# ls()

# view our dataset
view(biomarker_clean)

# randomly take a sample of proteins
# get col names
n_cols_to_select <- 5 # Example: select 3 columns
selected_df_dplyr <- biomarker_clean %>% 
  select(-group, -ados) %>%
  select(where(is.numeric)) %>%
  select(sample(seq_len(ncol(.)), size = n_cols_to_select))

# turn dataframe into : protein, value df. 
long_df <- selected_df_dplyr %>%
  pivot_longer(cols = everything(), names_to = "protein", values_to = "value")


ggplot(long_df, aes(x = value)) +
  geom_histogram(bins = 30, color = "black", na.rm = TRUE) +
  facet_wrap(~ protein, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Binned frequency of protein values", x = "Value range", y = "Count")


ggplot(long_df, aes(x = value)) +
  geom_histogram(bins = 30, color = "black", na.rm = TRUE) +
  facet_wrap(~ protein, scales = "free_x") +
  scale_y_log10() +  # log scale for counts
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Frequency of protein values (log scale counts)", x = "Value", y = "Count (log scale)")

# log transforming seems to make the distribution more rounded and normal, which means we can assume normal distribution behavior. 
# QUESTION 2


# Only numeric columns for counting
protein_cols <- biomarker_clean %>%
  select(-group, -ados) %>%
  select(where(is.numeric)) %>%
  colnames()
view(biomarker_clean)
# Add row ID and count 3 or -3 in each row
outliers <- biomarker_clean %>%
  mutate(id = row_number(),
         n = rowSums(across(all_of(protein_cols), 
                                     ~ . >= 3 | . <= -3))) %>%
  select(id, group, n)

view(outliers)


# QUESTION 3


## run exactly data from other script
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')

# partition into training and test set
set.seed(123)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.7)

biomarker.train <- training(biomarker_split)
biomarker.test <- testing(biomarker_split)

setwd("/Users/katiepyo/Desktop/PSTAT/PSTAT197/module-1-biomarker-data-table1")

load('data/biomarker-clean.RData')
setwd("/Users/katiepyo/Desktop/PSTAT/PSTAT197/module-1-biomarker-data-table1")
## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker.train %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker.train %>%
  select(-c(group, ados))

response <- biomarker.train %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)


biomarker.train.transform <- biomarker.train %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

biomarker.test.transform <- biomarker.test %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)



# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = biomarker.train.transform, 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

biomarker.test.transform %>%
  add_predictions(fit, type = 'response') %>%
  mutate(estimate = factor(pred > 0.5),
         truth = factor(class)) %>% 
  class_metrics(estimate = estimate,
                truth = truth, pred,
                event_level = 'second')


# excluding somehting
# A tibble: 4 × 3
#.metric     .estimator .estimate
#<chr>       <chr>          <dbl>
#  1 sensitivity binary         0.76 
#2 specificity binary         0.682
#3 accuracy    binary         0.723
#4 roc_auc     binary         0.793




