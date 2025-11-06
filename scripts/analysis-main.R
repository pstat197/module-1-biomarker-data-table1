
# ------------------ ## QUESTION 1 ## --------------- 

#load raw data
load('/Users/kaeya/PSTAT 197A/module-1-biomarker-data-table1/data/biomarker-clean.RData')
rawbiodata <- read.csv('/Users/kaeya/PSTAT 197A/module-1-biomarker-data-table1/data/biomarker-raw.csv',
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
var_names <- read_csv('/Users/kaeya/PSTAT 197A/module-1-biomarker-data-table1/data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# clean & normalize data
biomarker_clean <- read_csv('/Users/kaeya/PSTAT 197A/module-1-biomarker-data-table1/data/biomarker-raw.csv', 
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

# ------------------ ## QUESTION 4 ## --------------- 

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










