library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)


# Question 3 - large n top predictive proteins + fuzzy interaction

# get names
var_names <- read_csv('C:/Users/prasa/Desktop/ds/pstat197/197a/module_1/module-1-biomarker-data-table1/data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean <- read_csv('C:/Users/prasa/Desktop/ds/pstat197/197a/module_1/module-1-biomarker-data-table1/data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  mutate(across(.cols = -c(group, ados), 
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  select(group, ados, everything())

set.seed(123456)
biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.7)

# Create training and testing sets
train_data <- training(biomarker_split)
test_data <- testing(biomarker_split)

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

ttests_out <- train_data %>%
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

# select top proteins (n = 20)
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 20) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- train_data %>%
  select(-c(group, ados))

response <- train_data %>% pull(group) %>% factor()

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
  slice_max(MeanDecreaseGini, n = 20) %>%
  pull(protein)

## FUZZY INTERSECTION
####################

# Keep proteins appearing in both selection methods
all_proteins <- list(proteins_s1, proteins_s2)
protein_counts <- table(unlist(all_proteins))
fuzzy_proteins <- names(protein_counts[protein_counts >= 2])

## LOGISTIC REGRESSION
#######################

# Prepare training data for logistic regression
train_fuzzy <- train_data %>%
  select(group, any_of(fuzzy_proteins)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# Prepare test data with same proteins
test_fuzzy <- test_data %>%
  select(group, any_of(fuzzy_proteins)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = train_fuzzy, 
           family = 'binomial')

# evaluate on the held-out test set (30%)
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

# Create predictions for test set
test_predictions <- test_fuzzy %>%
  mutate(
    pred_prob = predict(fit, newdata = ., type = "response"),
    pred_class = factor(pred_prob > 0.5, levels = c(FALSE, TRUE), labels = c("TD", "ASD")),
    truth = factor(class, levels = c(FALSE, TRUE), labels = c("TD", "ASD"))
  )

# Calculate metrics
test_results <- test_predictions %>%
  class_metrics(
    truth = truth,
    estimate = pred_class,
    pred_prob,
    event_level = 'second'
  )

# Print test set results
print("Model performance on test set (fuzzy intersection):")
print(test_results)