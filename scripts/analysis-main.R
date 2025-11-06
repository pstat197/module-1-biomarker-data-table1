
# ------------------ ## QUESTION 1 ## --------------- 

#load raw data
load('/data/biomarker-clean.RData')
rawbiodata <- read.csv('/data/biomarker-raw.csv',
                       header = TRUE, 
                       skip = 1)

# change headings
rawbiodata <- subset(rawbiodata, select = -Target)
view(rawbiodata)
rawsampledata <- rawbiodata %>% sample_n(size = 50)
view(rawsampledata)

# label the different axes
bio_long <- rawsampledata %>%
  select(X, CEBPB, CHIP, NSE, PIAS4) %>%  
  pivot_longer(cols = -X,
               names_to = "Protein",
               values_to = "Value")
view(bio_long)

# without log transformation
ggplot(bio_long, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "gray", color = "white", na.rm = TRUE) +
  facet_wrap(~Protein, scales = "free") +
  theme_minimal() +
  labs(title = "Distributions of protein levels before log transformation",
       x = "# of Protiens")

# with log transformation
ggplot(bio_long, aes(x = log(Value))) +
  geom_histogram(bins = 30, fill = "gray", color = "white", na.rm = TRUE) +
  facet_wrap(~Protein, scales = "free") +
  theme_minimal() +
  labs(title = "Distributions of protein levels before log transformation",
       x = "# of Protiens")


# ------------------ ## QUESTION 2 ## ---------------

# get names
var_names <- read_csv('/data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# clean & normalize data
biomarker_clean <- read_csv('/data/biomarker-raw.csv', 
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
                ~ as.numeric(scale(log10(.x))[, 1]))) %>%
  # reorder columns
  select(group, ados, everything())

# protein column names
protein_cols <- pull(var_names, abbreviation)
protein_cols

# add sample ID (represents customers)
bio_z <- biomarker_clean %>%
  mutate(SampleID = row_number())

# tabulate the outliers for each subject
subject_outliers <- bio_z %>%
  mutate(n_outliers = rowSums(across(all_of(protein_cols),
                                     ~ is.finite(.x) & abs(.x) > 3))) %>%
  select(SampleID, group, n_outliers)

view(subject_outliers)

# boxplot for ASD vs TD
ggplot(subject_outliers, aes(x = group, y = n_outliers, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of Outlier Counts (ASD vs TD)",
       x = "Group",
       y = "Number of Outlying Proteins (|z| > 3)") +
  scale_fill_manual(values = c("ASD" = "gray", "TD" = "gray")) +
  theme(legend.position = "none")


# ------------------ ## QUESTION 3A ## ---------------

load('/data/biomarker-clean.RData')
biomarker_clean

set.seed(123)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.7)
biomarker.train <- training(biomarker_split)
biomarker.test <- testing(biomarker_split)

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

proteins_s1

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

proteins_s2

## LOGISTIC REGRESSION
#######################
set.seed(123)
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

fit

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

finalpanel1 <- names(coef(fit)[-1])



# ------------------ ## QUESTION 3B ## ---------------

load('/data/biomarker-clean.RData')
biomarker_clean

set.seed(123)
n <- nrow(biomarker_clean)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.7)
biomarker.train <- training(biomarker_split)
biomarker.test <- testing(biomarker_split)

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
  slice_min(p.adj, n = 15) %>%
  pull(protein)

proteins_s1

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
  slice_max(MeanDecreaseGini, n = 15) %>%
  pull(protein)

proteins_s2

## LOGISTIC REGRESSION
#######################
set.seed(123)
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

fit

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

finalpanel2 <- names(coef(fit)[-1])


# ------------------ ## QUESTION 3C ## ---------------

load('/data/biomarker-clean.RData')
biomarker_clean

set.seed(123)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.7)
biomarker.train <- training(biomarker_split)
biomarker.test <- testing(biomarker_split)

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

proteins_s1

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

proteins_s2

## LOGISTIC REGRESSION
#######################
set.seed(123)
# select subset of interest
proteins_sstar <- union(proteins_s1, proteins_s2)

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

fit

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

summary(fit)

finalpanel3 <- names(coef(fit)[-1])
finalpanel3


# ----------------------- ## QUESTION 4 ## ------------------ 

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')
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

ttests_out <- biomarker_clean %>%
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
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

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
proteins_sstar <- intersection(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')
# View the summary of fit, see if there are any proteins with high p values
summary(fit)

# keep only proteins with low p values
fit2 <- glm(class ~ -1 + DERM + RELT + IgD, 
           data = training(biomarker_split), 
           family = 'binomial')
summary(fit2)

# compare fit vs. fit2 performance
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)


testing(biomarker_split) %>%
  add_predictions(fit2, type = 'response') %>%
  mutate(estimate = factor(pred > 0.5),
         truth = factor(class)) %>% 
  class_metrics(estimate = estimate,
                truth = truth, pred,
                event_level = 'second')

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(estimate = factor(pred > 0.5),
         truth = factor(class)) %>% 
  class_metrics(estimate = estimate,
                truth = truth, pred,
                event_level = 'second')










