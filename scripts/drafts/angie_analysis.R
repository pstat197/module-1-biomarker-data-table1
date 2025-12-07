# Question 03

library(tidyverse)
library(randomForest)
library(glmnet)
library(pROC)

df <- biomarker_z %>%
  select(-ados) %>% # don’t use ADOS for prediction
  mutate(group = factor(group))  

set.seed(123)
train_idx <- sample(seq_len(nrow(df)), size = 0.7 * nrow(df))

train <- df[train_idx, ]
test  <- df[-train_idx, ]

corr_scores <- train %>% 
  summarize(across(
    -group,
    ~ abs(cor(as.numeric(group == 'ASD'), .x, use='pairwise.complete.obs'))
  )) %>% 
  pivot_longer(everything(), names_to = 'protein', values_to = 'score') %>% 
  arrange(desc(score))

ttest_fn <- function(x) {
  t.test(x ~ train$group)$statistic %>% abs()
}

ttest_scores <- train %>%
  summarize(across(
    -group,
    ttest_fn
  )) %>%
  pivot_longer(everything(), names_to = "protein", values_to = "score") %>%
  arrange(desc(score))

x_train <- train %>% select(-group)
y_train <- train$group

rf_model <- randomForest(
  x = x_train, 
  y = y_train,
  importance = TRUE
)

rf_imp <- importance(rf_model)

rf_scores <- data.frame(
  protein = rownames(rf_imp),
  score = rf_imp[, 'MeanDecreaseGini']
  ) %>% 
  arrange(desc(score))

head(rf_scores)

logit_scores <- train %>%
  summarize(across(
    -group,
    ~ abs(coef(glm(train$group ~ .x, family = "binomial"))[2])
  )) %>%
  pivot_longer(everything(), names_to = "protein", values_to = "score") %>%
  arrange(desc(score))

top_n <- 20
top_corr  <- corr_scores  %>% slice_head(n = top_n)
top_ttest <- ttest_scores %>% slice_head(n = top_n)
top_rf    <- rf_scores    %>% slice_head(n = top_n)
top_logit <- logit_scores %>% slice_head(n = top_n)

all_methods <- bind_rows(
  mutate(top_corr,  method = "corr"),
  mutate(top_ttest, method = "ttest"),
  mutate(top_rf,    method = "rf"),
  mutate(top_logit, method = "logit")
)

fuzzy_panel <- all_methods %>%
  count(protein) %>%
  filter(n >= 2) %>% # appears in at least 2 methods
  pull(protein)

fuzzy_panel

fit <- glm(group ~ ., data = train %>% select(group, all_of(fuzzy_panel)),
           family = "binomial")

test_pred <- predict(fit, newdata = test, type = "response")

auc_fuzzy <- roc(test$group, test_pred)$auc
auc_fuzzy

# compare hard vs fuzzy intersection and different top-N choices

hard_panel <- all_methods %>%
  count(protein) %>%
  filter(n == 4) %>% # appears in all 4 selection methods
  pull(protein)

hard_panel
length(hard_panel)

top_n10 <- 10

top_corr10  <- corr_scores  %>% slice_head(n = top_n10)
top_ttest10 <- ttest_scores %>% slice_head(n = top_n10)
top_rf10    <- rf_scores    %>% slice_head(n = top_n10)
top_logit10 <- logit_scores %>% slice_head(n = top_n10)

all_methods10 <- bind_rows(
  mutate(top_corr10,  method = "corr"),
  mutate(top_ttest10, method = "ttest"),
  mutate(top_rf10,    method = "rf"),
  mutate(top_logit10, method = "logit")
)

fuzzy_panel10 <- all_methods10 %>%
  count(protein) %>%
  filter(n >= 2) %>%
  pull(protein)

fuzzy_panel10
length(fuzzy_panel10)

library(pROC)

compute_auc <- function(panel) {
  if (length(panel) == 0) return(NA_real_)
  
  fit <- glm(
    group ~ .,
    data = train %>% select(group, all_of(panel)),
    family = "binomial"
  )
  
  test_pred <- predict(fit, newdata = test, type = "response")
  
  roc(test$group, test_pred)$auc
}

auc_hard        <- compute_auc(hard_panel)
auc_fuzzy20     <- compute_auc(fuzzy_panel)
auc_fuzzy10     <- compute_auc(fuzzy_panel10)

results_models <- tibble::tibble(
  model        = c("Hard intersection (top 20)", 
                   "Fuzzy intersection (top 20)", 
                   "Fuzzy intersection (top 10)"),
  n_proteins   = c(length(hard_panel),
                   length(fuzzy_panel),
                   length(fuzzy_panel10)),
  auc          = c(auc_hard,
                   auc_fuzzy20,
                   auc_fuzzy10)
)

results_models

# Question 04: simpler model

# simpler (LASSO) method
library(glmnet)

# design matrix and response for training
x_train <- model.matrix(group ~ ., data = train)[, -1]
y_train <- train$group

set.seed(123)
cv_fit <- cv.glmnet(
  x_train,
  y_train,
  family = "binomial",
  alpha = 1 # LASSO
)

lasso_coef <- coef(cv_fit, s = "lambda.min")

selected_idx <- which(lasso_coef != 0)[-1]  # drop intercept
lasso_panel <- rownames(lasso_coef)[selected_idx]

lasso_panel
length(lasso_panel)

x_test <- model.matrix(group ~ ., data = test)[, -1]

lasso_pred <- predict(cv_fit, newx = x_test, s = "lambda.min", type = "response")

auc_lasso <- roc(test$group, as.numeric(lasso_pred))$auc
auc_lasso

results_models <- results_models %>%
  add_row(
    model      = "LASSO panel",
    n_proteins = length(lasso_panel),
    auc        = auc_lasso
  )

results_models

# write-up

'i repeated the feature-selection using a 70/30 train/test split, performing
the analysis on the training set only and evaluating final models on the test
set. compared with the in-class analysis, the test AUC values were lower
and more conservative.

when increasing the number of top proteins from 10 to 20, i found that the
fuzzy intersection panel (proteins selected by at least 2/4 methods)
went from 10 to 23 proteins. however, the larger panel actually produced worse
test performance (with an AUC of 0.651) compared to the smaller fuzzy panel
(which had an AUC of 0.742). the hard intersection panel (proteins 
selected by all four methods) based on the top 20 per method contained only 7 
proteins and had an AUC of 0.736.adding more proteins tended to seem to 
increase noise rather than help with prediction, reducing the accuracy overall. 

comparing the intersection strategies, it seems like the hard intersection 
output a small panel but had good performance overall. on the other hand,
the fuzzy intersection including the top 10 per method had a slightly higher 
AUC than even the hard intersection panel. the top-10 fuzzy approach seems to 
overall be the better choice compared with the other intersections. overall,
these results seem to highlight the lack of stability when it comes to
biomarker selection in higher dimensions but lower sample size situations. 
small changes to intersection criteria, etc. can greatly affect the final 
panel and its predictive performance.


in terms of the simpler panel, i used LASSO logistic regression. the resulting
model selected 47 proteins but still achieved the strongest AUC of 0.7582. This
means that it did better than even the fuzzy-10 panel. Although larger, LASSO 
panels provide improved classification accuracy and presents a good alternative
to intersection-based approaches. Overall, since LASSO looks at multiple 
proteins at once, it can decide which combination works best. the earlier 
intersection methods looked at each protein independently, so it couldnt 
properly capture these combined effects.
