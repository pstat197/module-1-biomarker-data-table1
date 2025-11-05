#QUESTION 1
raw_data <- read.csv("data/biomarker-raw.csv", header = TRUE, skip = 1)

raw_data <- raw_data %>% select(-Target)

bio_long <- raw_data %>%
  pivot_longer(cols = CHIP:PIAS4,
               names_to = "protein",
               values_to = "level")

ggplot(bio_long, aes(x = level)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", na.rm = TRUE) +
  facet_wrap(~protein, scales = "free") +
  theme_minimal() +
  labs(title = "Protein level distributions (raw scale)",
       x = "Protein measurement", y = "Count")

ggplot(bio_long, aes(x = log(level))) +
  geom_histogram(bins = 30, fill = "gray40", color = "white", na.rm = TRUE) +
  facet_wrap(~protein, scales = "free") +
  theme_minimal() +
  labs(title = "Protein level distributions (log scale)",
       x = "Log-transformed measurement", y = "Count")


# QUESTION 3 
library(tidyverse)
library(randomForest)

load("data/biomarker-clean.RData")

set.seed(123)
n <- nrow(biomarker_clean)
train_index <- sample(1:n, size = 0.7 * n)
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

x_train <- train_data %>% select(-c(group, ados))
y_train <- factor(train_data$group)

set.seed(101422)
rf_model <- randomForest(x = x_train, y = y_train,
                         ntree = 500, importance = TRUE)

rf_importance <- as.data.frame(rf_model$importance)
rf_importance$protein <- rownames(rf_importance)

top_rf <- rf_importance %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = 20) %>%
  pull(protein)

fuzzy_shared <- unique(unlist(
  lapply(top_ttest, function(name1) {
    matches <- agrep(name1, top_rf, max.distance = 0.1, ignore.case = TRUE)
    if (length(matches) > 0) return(name1)
    else return(NULL)
  })
))

cat("Number of fuzzy-matched shared proteins:", length(fuzzy_shared), "\n")
print(fuzzy_shared)

shared_proteins <- fuzzy_shared

train_sub <- train_data %>%
  select(group, all_of(shared_proteins)) %>%
  mutate(group = ifelse(group == "ASD", 1, 0))

test_sub <- test_data %>%
  select(group, all_of(shared_proteins)) %>%
  mutate(group = ifelse(group == "ASD", 1, 0))

fit <- glm(group ~ ., data = train_sub, family = "binomial")

pred_probs <- predict(fit, newdata = test_sub, type = "response")
pred_class <- ifelse(pred_probs > 0.5, 1, 0)

confusion <- table(True = test_sub$group, Predicted = pred_class)
print(confusion)

accuracy <- mean(pred_class == test_sub$group)
cat("Test accuracy:", round(accuracy, 3), "\n")
