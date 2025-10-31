
# ------------------ ## QUESTION 1 ## --------------- 

#load raw data
load('/Users/kaeya/PSTAT 197A/module-1-biomarker-data-table1/data/biomarker-clean.RData')
rawbiodata <- read.csv('/Users/kaeya/PSTAT 197A/module-1-biomarker-data-table1/data/biomarker-raw.csv',
                       header = TRUE, 
                       skip = 1)

# change headings
rawbiodata <- subset(rawbiodata, select = -Target)
view(rawbiodata)

# label the different axes
bio_long <- rawbiodata %>%
  pivot_longer(cols = CHIP:PIAS4, names_to = "Protein", values_to = "Value")
bio_long
view(bio_long)

# without log transformation
ggplot(bio_long, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", na.rm = TRUE) +
  facet_wrap(~Protein, scales = "free") +
  theme_minimal() +
  labs(title = "Distributions of protein levels before log transformation",
       x = "# of Protiens")

# with log transformation
ggplot(bio_long, aes(x = log(Value))) +
  geom_histogram(bins = 30, fill = "gray", color = "white", na.rm = TRUE) +
  facet_wrap(~Protein, scales = "free") +
  theme_minimal() +
  labs(title = "Distributions of protein levels before log transformation")


# ------------------ ## QUESTION 2 ## ---------------



