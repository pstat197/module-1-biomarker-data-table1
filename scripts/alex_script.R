#QUESTION 3

library(tidyverse)
library(randomForest)

# load data
load("team assignment_two/module-1-biomarker-data-table1/data/biomarker-clean.RData")

#Train/test split
set.seed(123)
n <- nrow(biomarker_clean)
train_index <- sample(1:n, size = 0.8 * n)  
train_data <- biomarker_clean[train_index, ]
test_data  <- biomarker_clean[-train_index, ]

protein_names <- setdiff(names(train_data), c("group", "ados"))

t_results <- map_dfr(protein_names, function(p) {
  t <- t.test(train_data[[p]] ~ train_data$group)
  tibble(protein = p, p_value = t$p.value)
}) %>%
  arrange(p_value)

top_ttest <- t_results %>%
  slice_head(n = 20) %>%    
  pull(protein)

#Random forest on train
# predictors and response
x_train <- train_data %>% select(-c(group, ados))
y_train <- factor(train_data$group)

set.seed(123)
rf_model <- randomForest(x = x_train, y = y_train, ntree = 500, importance = TRUE)

# top 20 proteins
rf_importance <- as.data.frame(rf_model$importance)
rf_importance$protein <- rownames(rf_importance)

top_rf <- rf_importance %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = 20) %>%
  pull(protein)

shared_proteins <- intersect(top_ttest, top_rf)

cat("Number of shared proteins:", length(shared_proteins), "\n")
print(shared_proteins)


# Make binary class variable (1 = ASD, 0 = TD)
train_sub <- train_data %>%
  select(group, all_of(shared_proteins)) %>%
  mutate(group = ifelse(group == "ASD", 1, 0))

test_sub <- test_data %>%
  select(group, all_of(shared_proteins)) %>%
  mutate(group = ifelse(group == "ASD", 1, 0))

# fit model
fit <- glm(group ~ ., data = train_sub, family = "binomial")

#evaluate on test set
pred_probs <- predict(fit, newdata = test_sub, type = "response")
pred_class <- ifelse(pred_probs > 0.5, 1, 0)

# confusion matrix
confusion <- table(True = test_sub$group, Predicted = pred_class)
print(confusion)

#accuracy
accuracy <- mean(pred_class == test_sub$group)
cat("Test accuracy:", round(accuracy, 3), "\n")

cat("\nSummary:\n")
cat("- Selection done only on training data (no data leakage)\n")
cat("- Top 20 proteins chosen by each method\n")
cat("- Logistic model built on intersection of selected proteins\n")
cat("- Final accuracy evaluated on untouched test set\n")