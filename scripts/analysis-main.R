
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











