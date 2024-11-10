#--CHIP ANALYSIS----------------------------------------------------------------
#--CHIP, brain lesions and cognitive performance in AF--------------------------
# Author: Pascal B. Meyre
# Date: 07/27/24
# Location: Reinach, Baselland, Switzerland

# Use combined dataset of SWISS-and BEAT-AF pts.
# Use CHIP dataset of SWISS-and BEAT-AF pts.
#-------------------------------------------------------------------------------

################################################################################
# Merge CHIP and SWISS/BEAT-AF datasets
################################################################################

# Import SWISS/BEAT-AF dataset
library(data.table)
swiss.beat <- fread("/Users/pascalmeyre/Desktop/Research/1_Projects_Analysis/19_CHIP_brain_lesions_cognition/analysis/datasets/chip.20241003.csv")

# Make a dataset only including CHIP variable and pat_id
library(dplyr)

chip <- swiss.beat %>%
  select(pat.id, chip.status.carrier, nb.chip.driver.genes, chip.genes)

# Merge with the datset of Swiss-and BEAT-AF that includes MRI brain lesions
# Import datafile
swiss.beat.old <- fread("/Users/pascalmeyre/Desktop/Research/1_Projects_Analysis/19_CHIP_brain_lesions_cognition/analysis/datasets/swiss.beat.20240315.csv")
dat <- merge(swiss.beat.old, chip, by = "pat.id")

# Merge with dataset including alcohol consumption
CHIP.outcomes <- fread("/Users/pascalmeyre/Desktop/Research/1_Projects_Analysis/19_CHIP_brain_lesions_cognition/analysis/datasets/chip.20241003.csv")

alcohol <- CHIP.outcomes %>%
  select(pat.id, total.alcohol.consumption)

dat.final <- merge(dat, alcohol, by = "pat.id")

# Make a dataset with only Swiss-AF patients
dat_swiss <- subset(dat.final, study == "Swiss")

# Make a new with the dataset from 20241003
dat_swiss_filtered <- merge(dat_swiss_filtered, swiss.beat, by = "pat.id")

# Exlcude patients who did not undergo MRI screening
dat_swiss_filtered <- dat_swiss[!(
  is.na(dat_swiss$mb.nr) &
    is.na(dat_swiss$snci.count) &
    is.na(dat_swiss$snci.vol) &
    is.na(dat_swiss$lncci.count) &
    is.na(dat_swiss$lncci.vol) &
    is.na(dat_swiss$mod.faz) &
    is.na(dat_swiss$wml.count) &
    is.na(dat_swiss$wml.vol)
), ]

colnames(dat_swiss_filtered)[133] <- "CHIP_carrier"
colnames(dat_swiss_filtered)[134] <- "nr_CHIP_genes"

###
colnames(dat)[133] <- "CHIP_carrier"
colnames(dat)[134] <- "nr_CHIP_genes"

# Frequency of each CHIP mutation
chip_filtered <- subset(dat_swiss_filtered, CHIP_carrier == 1)
genes_list <- strsplit(as.character(chip_filtered$chip.genes), ", ")
genes_unlist <- unlist(genes_list)
genes_freq <- table(genes_unlist)
genes_freq_df <- as.data.frame(genes_freq)
colnames(genes_freq_df) <- c("Gene", "Frequency")
genes_freq_df <- genes_freq_df[order(-genes_freq_df$Frequency), ]
genes_freq_df$Gene <- factor(genes_freq_df$Gene, levels = genes_freq_df$Gene)

library(ggplot2)
ggplot(genes_freq_df, aes(x = reorder(Gene, -Frequency), y = Frequency, fill = Gene)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Frequency), vjust = -0.5, size = 3.5) +  # Add numbers on top of bars
  labs(x = "CHIP Genes", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")) +
  scale_fill_manual(values = scales::hue_pal()(length(genes_freq_df$Gene)))

# Frequency of patients with multiple CHIP mutations
ggplot(chip_filtered, aes(x = factor(nr_CHIP_genes), y = (..count..)/sum(..count..), fill = factor(nr_CHIP_genes))) +
  geom_bar() +
  geom_text(aes(label = scales::percent((..count..)/sum(..count..), accuracy = 0.1)),
            stat = "count", 
            vjust = -0.5, 
            size = 3) +  # Display percentage labels above the bars
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  labs(x = "Number of CHIP-driver Genes", y = "Percentage") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black")
  ) +
  scale_fill_manual(values = scales::hue_pal()(length(unique(chip_filtered$nr_CHIP_genes))))

# Number of CHIP mutations per individual across quartiles of age
dat_swiss_filtered$age_quartiles <- cut(dat_swiss_filtered$age.bl, breaks = quantile(dat_swiss_filtered$age.bl, probs = seq(0, 1, 0.25), na.rm = TRUE),
                         include.lowest = TRUE,
                         labels = c("44.3-67.7", "67.7-72.8", "72.8-77.9", "78.0-92.2"))

library(dplyr)
quartile_ranges <- dat_swiss_filtered %>%
  group_by(age_quartiles) %>%
  summarise(Min_Age = min(age.bl, na.rm = TRUE), 
            Max_Age = max(age.bl, na.rm = TRUE))
print(quartile_ranges)

# Convert nr_CHIP_genes to factor for better control of the order in the plot
dat_swiss_filtered$nr_CHIP_genes <- factor(dat_swiss_filtered$nr_CHIP_genes, levels = c(0, 1, 2, 3, 4),
                            labels = c("0", "1", "2", "3", "≥4"))

dat_swiss_filtered$nr_CHIP_genes <- factor(dat_swiss_filtered$nr_CHIP_genes, levels = c("≥4", "3", "2", "1", "0"))

ggplot(dat_swiss_filtered, aes(x = age_quartiles, fill = nr_CHIP_genes)) +
  geom_bar(position = "fill", aes(y = (..count..)/sum(..count..))) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Age quartiles (years)", y = "% of patients", fill = "No. of CHIP genes") +
  theme_bw() +
  theme(legend.position = "right") + 
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black"))

# Prevalence of DNMT3A, TET2 and other genes across age quartiles
dat_swiss_filtered$DNMT3A <- ifelse(grepl("DNMT3A", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$TET2 <- ifelse(grepl("TET2", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$Other_genes <- ifelse(dat_swiss_filtered$chip.status.carrier == 1 & !grepl("DNMT3A|TET2", dat_swiss_filtered$chip.genes), 1, 0)

library(ggplot2)
library(dplyr)
library(tidyr)

# Convert DNMT3A, TET2, and other_genes to factors for proper labeling
dat_swiss_filtered$DNMT3A <- factor(dat_swiss_filtered$DNMT3A, levels = c(0, 1), labels = c("No DNMT3A", "DNMT3A"))
dat_swiss_filtered$TET2 <- factor(dat_swiss_filtered$TET2, levels = c(0, 1), labels = c("No TET2", "TET2"))
dat_swiss_filtered$Other_genes <- factor(dat_swiss_filtered$Other_genes, levels = c(0, 1), labels = c("No Other Genes", "Other Genes"))

# Create a summary table for stacked bar chart
stacked_data <- dat_swiss_filtered %>%
  group_by(age_quartiles) %>%
  summarise(DNMT3A = mean(DNMT3A == "DNMT3A") * 100,    # Proportion within each quartile
            TET2 = mean(TET2 == "TET2") * 100,
            Other_genes = mean(Other_genes == "Other Genes") * 100) %>%
  pivot_longer(cols = c(DNMT3A, TET2, Other_genes),
               names_to = "Gene_Type", 
               values_to = "Percentage")

stacked_data <- stacked_data %>%
  group_by(age_quartiles) %>%
  arrange(age_quartiles, desc(Gene_Type)) %>%
  mutate(cumulative_percentage = cumsum(Percentage) - (Percentage / 2))

stacked_data$Gene_Type <- factor(stacked_data$Gene_Type, levels = c("DNMT3A", "TET2", "Other_genes"))

# Plot the stacked bar chart
ggplot(stacked_data, aes(x = age_quartiles, y = Percentage, fill = Gene_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cumulative_percentage, label = paste0(round(Percentage, 1), "%")), 
            color = "black", size = 3.5) +  # Adjust size as needed
  labs(x = "Age quartiles (years)", y = "Prevalence", fill = NULL) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_bw() +
  theme(legend.position = "right") + 
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black"))

# Generate the sex variable (male=1, female=2)
dat_swiss_filtered$pat.sex <- ifelse(dat_swiss_filtered$pat.sex =="Male", 1, 2)

library(patchwork)
library(forestmodel)

# Logistic regression for CHIP carrier across age quartiles
plot_chip <- forest_model(glm(CHIP_carrier ~ age_quartiles + pat.sex, data = dat_swiss_filtered, family = binomial)) +
  labs(title = "CHIP_carrier")

plot_dnmt3a <- forest_model(glm(DNMT3A ~ age_quartiles + pat.sex, data = dat_swiss_filtered, family = binomial)) +
  labs(title = "DNMT3A")

plot_tet2 <- forest_model(glm(TET2 ~ age_quartiles + pat.sex, data = dat_swiss_filtered, family = binomial)) +
  labs(title = "TET2")

plot_other <- forest_model(glm(Other_genes ~ age_quartiles + pat.sex, data = dat_swiss_filtered, family = binomial)) +
  labs(title = "Other_genes")

combined_plot1 <- plot_chip / plot_dnmt3a
combined_plot1

combined_plot2 <- plot_tet2 / plot_other
combined_plot2

# Preparation of covariates
# Generate the sex variable (male=1, female=2)
dat_swiss_filtered$pat.sex <- ifelse(dat_swiss_filtered$pat.sex =="Male", 1, 2)

# Smoking status
dat_swiss_filtered$rauchen <- ifelse(dat_swiss_filtered$rauchen == "nie geraucht", 0,
                                     ifelse(dat_swiss_filtered$rauchen == "fr\xfcher", 1,
                                            ifelse(dat_swiss_filtered$rauchen == "Ja, aktiv", 2, NA)))

# Current smoker
dat_swiss_filtered$current.smoker <- ifelse(dat_swiss_filtered$rauchen == 2, 1, 0)
num_missing_current.smoker <- sum(is.na(dat_swiss_filtered$current.smoker))
print(num_missing_current.smoker)
dat_swiss_filtered$current.smoker <- ifelse(is.na(dat_swiss_filtered$current.smoker), 0, dat_swiss_filtered$current.smoker)

# BMI
dat_swiss_filtered$bmi <- dat_swiss_filtered$gewicht/(dat_swiss_filtered$groesse/100)^2
dat_swiss_filtered$bmi <- as.numeric(dat_swiss_filtered$bmi)

# Hypertension
num_missing_hypertonie <- sum(is.na(dat_swiss_filtered$prev.hypertonie))
print(num_missing_hypertonie)
dat_swiss_filtered$prev.hypertonie <- ifelse(is.na(dat_swiss_filtered$prev.hypertonie), 0, dat_swiss_filtered$prev.hypertonie)

# Diabetes
num_missing_diabetes <- sum(is.na(dat_swiss_filtered$prev.diabetes))
print(num_missing_diabetes)
dat_swiss_filtered$prev.diabetes <- ifelse(is.na(dat_swiss_filtered$prev.diabetes), 0, dat_swiss_filtered$prev.diabetes)

# Prior stroke/TIA
num_missing_prior.stroke <- sum(is.na(dat_swiss_filtered$prev.stroke.tia))
print(num_missing_prior.stroke)
dat_swiss_filtered$prev.stroke.tia <- ifelse(is.na(dat_swiss_filtered$prev.stroke.tia), 0, dat_swiss_filtered$prev.stroke.tia)

# Oral anticoagulation
dat_swiss_filtered$med.oak.yn <- ifelse(dat_swiss_filtered$med.oak.yn =="No", 0, 1)

num_missing_med.oak.yn <- sum(is.na(dat_swiss_filtered$med.oak.yn))
print(num_missing_med.oak.yn)
dat_swiss_filtered$med.oak.yn <- ifelse(is.na(dat_swiss_filtered$med.oak.yn), 0, dat_swiss_filtered$med.oak.yn)

# Chronic kidney disease
num_missing_prev.niereninsuff <- sum(is.na(dat_swiss_filtered$prev.niereninsuff))
print(num_missing_prev.niereninsuff)
dat_swiss_filtered$prev.niereninsuff <- ifelse(is.na(dat_swiss_filtered$prev.niereninsuff), 0, dat_swiss_filtered$prev.niereninsuff)


baseline_characteristics <- merge(dat_swiss_filtered, swiss.beat, by = "pat.id")

library(tableone)
vars <- c("age.bl", "pat.sex", "prev.heart.failure", "prev.schlaf.apnoe", "prev.hypertonie", 
          "prev.diabetes", "prev.pavk", "prev.niereninsuff", "prev.mi", 
          "prev.sys.embolism", "prev.stroke.tia", "chads", "bmi", "rr.sys.liegend", "rr.dia.liegend", 
          "heart.rate", "trigly", "lpa", "alat", "creatinine.gfr", "ldl", "hdl", "chol", "crphs")

# Specify the categorical variables
catVars <- c("pat.sex", "prev.heart.failure", "prev.schlaf.apnoe", "prev.hypertonie", 
             "prev.diabetes", "prev.pavk", "prev.niereninsuff", "prev.mi", 
             "prev.sys.embolism", "prev.stroke.tia")

# Create the Table 1
table1 <- CreateTableOne(vars = vars, 
                         data = baseline_characteristics, 
                         strata = "chip.status.carrier.y", 
                         factorVars = catVars, 
                         test = TRUE)

# Print the table with p-values
print(table1, showAllLevels = TRUE, pDigits = 3)

################################################################################
# Association of CHIP with biomarkers
################################################################################

# Import SWISS/BEAT-AF biomaker dataset
library(data.table)
swiss.beat <- fread("/Users/pascalmeyre/Desktop/Research/1_Projects_Analysis/18_Biomarkers_MACE_bleeding/analysis/datasets/swiss.beat.biomarker.20240202.csv")

CHIP_biomarkers <- merge(swiss.beat, dat_swiss_filtered, by = "pat.id")

# Box plots for IL6 according to CHIP mutations
library(ggplot2)
library(ggpubr)

# Create a box plot with jittered points for hsCRP according to CHIP_carrier
ggplot(CHIP_biomarkers, aes(x = as.factor(chip.status.carrier), y = crphs, fill = as.factor(chip.status.carrier))) +  # Map fill and color to CHIP_carrier
  geom_boxplot(outlier.shape = NA) +  # Box plot without outliers
  geom_jitter(width = 0.2, alpha = 0.3, color = "darkgray") +  # Add jittered points (dots)
  labs(x = "CHIP carrier", y = "hsCRP", fill = "CHIP carrier") +  # Label the axes and legend
  theme_bw() +  # Use a black-and-white theme
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text.x = element_text(size = 9, color = "black"),  # Adjust x-axis text
    axis.text.y = element_text(size = 9, color = "black"),  # Adjust y-axis text
    text = element_text(size = 14)  # Increase overall text size
  ) +
  stat_compare_means(method = "wilcox.test") +  # Set custom fill colors
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 20) +
  scale_x_discrete(labels = c("0" = "No", "1" = "Yes"))

################################################################################
# Create a box plot with jittered points for IL6 according to CHIP_carrier
ggplot(CHIP_biomarkers, aes(x = as.factor(CHIP_carrier), y = il6, fill = as.factor(CHIP_carrier))) +  # Map fill and color to CHIP_carrier
  geom_boxplot(outlier.shape = NA) +  # Box plot without outliers
  geom_jitter(width = 0.2, alpha = 0.3, color = "darkgray") +  # Add jittered points (dots)
  labs(x = "CHIP carrier", y = "IL6", fill = "CHIP carrier") +  # Label the axes and legend
  theme_bw() +  # Use a black-and-white theme
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text.x = element_text(size = 9, color = "black"),  # Adjust x-axis text
    axis.text.y = element_text(size = 9, color = "black"),  # Adjust y-axis text
    text = element_text(size = 14)  # Increase overall text size
  ) +
  stat_compare_means(method = "wilcox.test") +  # Set custom fill colors
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +  # Set custom point colors
  ylim(0, 20) +
  scale_x_discrete(labels = c("0" = "No", "1" = "Yes"))

ggplot(CHIP_biomarkers, aes(x = as.factor(DNMT3A), y = il6, fill = as.factor(DNMT3A))) +  # Map fill and color to CHIP_carrier
  geom_boxplot(outlier.shape = NA) +  # Box plot without outliers
  geom_jitter(width = 0.2, alpha = 0.3, color = "darkgray") +  # Add jittered points (dots)
  labs(x = "DNMT3A carrier", y = "IL6", fill = "DNMT3A carrier") +  # Label the axes and legend
  theme_bw() +  # Use a black-and-white theme
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text.x = element_text(size = 9, color = "black"),  # Adjust x-axis text
    axis.text.y = element_text(size = 9, color = "black"),  # Adjust y-axis text
    text = element_text(size = 14)  # Increase overall text size
  ) +
  stat_compare_means(method = "wilcox.test") +  # Set custom fill colors
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 20) +
  scale_x_discrete(labels = c("No DNMT3A" = "No", "DNMT3A" = "Yes"))

################################################################################
# Association between CHIP status and brain lesions
################################################################################

# Make a dataset with only Swiss-AF patients
dat_swiss <- subset(dat, study == "Swiss")

# Exlcude patients who did not undergo MRI screening
dat_swiss_filtered <- dat_swiss[!(
  is.na(dat_swiss$mb.nr) &
    is.na(dat_swiss$snci.count) &
    is.na(dat_swiss$snci.vol) &
    is.na(dat_swiss$lncci.count) &
    is.na(dat_swiss$lncci.vol) &
    is.na(dat_swiss$mod.faz) &
    is.na(dat_swiss$wml.count) &
    is.na(dat_swiss$wml.vol)
), ]

dat_swiss_filtered$chip.status.carrier <- as.factor(dat_swiss_filtered$chip.status.carrier)

na_count_snci_count <- sum(is.na(dat_swiss_filtered$snci.count))
# Exclude patients with snci.count==NA
dat_swiss_filtered <- dat_swiss_filtered[!is.na(dat_swiss_filtered$snci.count), ]

dat_swiss_filtered$DNMT3A <- ifelse(grepl("DNMT3A", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$TET2 <- ifelse(grepl("TET2", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$Other_genes <- ifelse(dat_swiss_filtered$chip.status.carrier == 1 & !grepl("DNMT3A|TET2", dat_swiss_filtered$chip.genes), 1, 0)

################################################################################
# Association between CHIP status and SNCI lesion
# Generate SNCI presence variable (1=yes, 0=no)
dat_swiss_filtered$snci_status <- ifelse(dat_swiss_filtered$snci.count > 0, 1, 0)

# Create Gene Category: No CHIP carrier, CHIP carrier, DNMT3A, TET2, Other_genes
dat_swiss_filtered <- dat_swiss_filtered %>%
  mutate(Gene_Category = case_when(
    CHIP_carrier == 0 ~ "No CHIP carrier",
    DNMT3A == 1 ~ "DNMT3A",
    TET2 == 1 ~ "TET2",
    Other_genes == 1 ~ "Other_genes",
    CHIP_carrier == 1 ~ "CHIP carrier"
  ))

# Create a summary table for the stacked bar chart
stacked_data_snci <- dat_swiss_filtered %>%
  group_by(Gene_Category, snci_status) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Gene_Category) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Compute cumulative percentage for label positioning
stacked_data_snci <- stacked_data_snci %>%
  group_by(Gene_Category) %>%
  arrange(Gene_Category, desc(snci_status)) %>%
  mutate(cumulative_percentage = cumsum(Percentage) - (Percentage / 2))

stacked_data_snci$Gene_Category <- factor(stacked_data_snci$Gene_Category,
                                          levels = c("No CHIP carrier", "DNMT3A", "TET2", "Other_genes"))

# Perform a trend test
# Create a contingency table of Gene_Category and SNCI status
contingency_table <- xtabs(Count ~ Gene_Category + snci_status, data = stacked_data_snci)
print(contingency_table)

chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)

p_value <- chi_square_result$p.value

# Plot the stacked bar chart
library(ggplot2)
library(ggsignif)

# Your plot with the significance bracket
ggplot(stacked_data_snci, aes(x = Gene_Category, y = Percentage, fill = factor(snci_status))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cumulative_percentage, label = paste0(round(Percentage, 1), "%")), 
            color = "black", size = 3.5, vjust = 1.5) + 
  labs(x = "", y = "% of patients", fill = "SNCI Presence") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_discrete(labels = c("No", "Yes")) +  # Use default colors and change labels
  theme_classic() +
  theme(legend.position = "top", 
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")) +
  
  # Adding the significance bracket
  geom_signif(comparisons = list(c("No CHIP carrier", "Other_genes")),
              map_signif_level = FALSE,
              y_position = 100,  # Place the bracket lower to avoid overlap
              tip_length = 0, textsize = 0) +
  
  # Adding p-value annotation above the bars
  annotate("text", x = 2.5, y = 105, label = "P=0.01718", size = 3.5, vjust = -0.5)

################################################################################
# Logistic regression models according to CHIP carrier status

dat_swiss_filtered$age_category <- ifelse(dat_swiss_filtered$age.bl >= 70, 1, 0)

# Univariable logistic regression
logistic_model <- glm(snci_status ~ chip.status.carrier, data = dat_swiss_filtered, family = binomial)
summary(logistic_model)

library(forestmodel)
print(forest_model(glm(snci_status ~ chip.status.carrier, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(snci_status ~ chip.status.carrier + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(snci_status ~ chip.status.carrier + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariale adjusted
print(forest_model(glm(snci_status ~ chip.status.carrier + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to DNMT3A carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(snci_status ~ DNMT3A, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(snci_status ~ DNMT3A + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(snci_status ~ DNMT3A + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(snci_status ~ DNMT3A + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to TET2 carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(snci_status ~ TET2, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(snci_status ~ TET2 + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(snci_status ~ TET2 + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(snci_status ~ TET2 + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to other genes carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(snci_status ~ Other_genes, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(snci_status ~ Other_genes + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(snci_status ~ Other_genes + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
dat_swiss_filtered$age_squared <- dat_swiss_filtered$age.bl^2
dat_swiss_filtered$age_standardized <- dat_swiss_filtered$age.bl / 10
dat_swiss_filtered$age_category <- ifelse(dat_swiss_filtered$age.bl >= 70, 1, 0)
print(forest_model(glm(snci_status ~ Other_genes + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Data for the forest plot
labels <- c("No CHIP carrier (Reference)", 
            "CHIP carrier", 
            "DNMT3A", 
            "TET2", 
            "Other genes")

n_values <- c("1222/1562", 
              "340/1562", 
              "167/1562", 
              "96/1562", 
              "99/1562")

or_with_ci <- c("", 
                "1.24 (0.93 - 1.65)", 
                "0.98 (0.65 - 1.45)", 
                "1.09 (0.66 - 1.74)", 
                "1.63 (1.02 - 2.54)")

p_values <- c(NA, "0.14", "0.93", "0.74", "0.03")

odds_ratios <- c(NA, 1.24, 0.98, 1.09, 1.63)
lower_ci <- c(NA, 0.93, 0.65, 0.66, 1.02)
upper_ci <- c(NA, 1.65, 1.45, 1.74, 2.54)

# Add column headers
header <- c("", "N", "OR (95% CI)", "P value")

# Create a matrix for the labels and data
label_matrix <- rbind(header, cbind(labels, n_values, or_with_ci, p_values))

library(forestplot)

# Create the forest plot
forestplot(
  labeltext = label_matrix,  # Use the matrix for labels
  mean = c(NA, odds_ratios),
  lower = c(NA, lower_ci),
  upper = c(NA, upper_ci),
  is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Mark headers as summary
  xlog = TRUE,  # No log scale for OR
  zero = 1,  # Reference line for OR = 1
  xlab = "Odds Ratio (95% CI)",
  boxsize = 0.2,  # Size of boxes for OR
  ci.vertices = TRUE,  # Show CI vertices
  ci.vertices.height = 0.05,  # Height of CI vertices
  txt_gp = forestplot::fpTxtGp(
    label = gpar(fontsize = 12),          # Font size for labels
    ticks = gpar(fontsize = 16),          # Font size for x-axis tick labels
    xlab = gpar(fontsize = 14)  # Font size and bold for x-axis label
  )
)

# Association between number of CHIP mutations and SNI lesion
print(forest_model(glm(snci_status ~ nb.chip.driver.genes + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Data for the forest plot
labels <- "Nr. of CHIP mutations"

n_values <- "1562"

or_with_ci <- "1.15 (0.92 - 1.42)"

p_values <- "0.20"

odds_ratios <- 1.15
lower_ci <- 0.92
upper_ci <- 1.42

# Add column headers
header <- c("", "N", "OR (95% CI)", "P value")

# Create a matrix for the labels and data
label_matrix <- rbind(header, cbind(labels, n_values, or_with_ci, p_values))

library(forestplot)

# Create the forest plot
forestplot(
  labeltext = label_matrix,  # Use the matrix for labels
  mean = c(NA, odds_ratios),
  lower = c(NA, lower_ci),
  upper = c(NA, upper_ci),
  is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Mark headers as summary
  xlog = TRUE,  # No log scale for OR
  zero = 1,  # Reference line for OR = 1
  xlab = "Odds Ratio (95% CI)",
  boxsize = 0.2,  # Size of boxes for OR
  ci.vertices = TRUE,  # Show CI vertices
  ci.vertices.height = 0.05,  # Height of CI vertices
  txt_gp = forestplot::fpTxtGp(
    label = gpar(fontsize = 12),          # Font size for labels
    ticks = gpar(fontsize = 16),          # Font size for x-axis tick labels
    xlab = gpar(fontsize = 14)            # Font size and bold for x-axis label
  ),
  xticks = seq(0.5, 2, by = 0.5)  # Define x-axis ticks from 0.5 to 2
)

################################################################################
# Association between CHIP status and LNCCI lesion

dat_swiss_filtered$CHIP_carrier <- as.factor(dat_swiss_filtered$chip.status.carrier)

na_count_lncci_count <- sum(is.na(dat_swiss_filtered$lncci.count))
# Exclude patients with lncci.count==NA
dat_swiss_filtered <- dat_swiss_filtered[!is.na(dat_swiss_filtered$lncci.count), ]

dat_swiss_filtered$DNMT3A <- ifelse(grepl("DNMT3A", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$TET2 <- ifelse(grepl("TET2", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$Other_genes <- ifelse(dat_swiss_filtered$CHIP_carrier == 1 & !grepl("DNMT3A|TET2", dat_swiss_filtered$chip.genes), 1, 0)

# Generate LNCCI presence variable (1=yes, 0=no)
dat_swiss_filtered$lncci_status <- ifelse(dat_swiss_filtered$lncci.count > 0, 1, 0)

# Create Gene Category: No CHIP carrier, CHIP carrier, DNMT3A, TET2, Other_genes
dat_swiss_filtered <- dat_swiss_filtered %>%
  mutate(Gene_Category = case_when(
    CHIP_carrier == 0 ~ "No CHIP carrier",
    DNMT3A == 1 ~ "DNMT3A",
    TET2 == 1 ~ "TET2",
    Other_genes == 1 ~ "Other_genes",
    CHIP_carrier == 1 ~ "CHIP carrier"
  ))

# Create a summary table for the stacked bar chart
stacked_data_lncci <- dat_swiss_filtered %>%
  group_by(Gene_Category, lncci_status) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Gene_Category) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Compute cumulative percentage for label positioning
stacked_data_lncci <- stacked_data_lncci %>%
  group_by(Gene_Category) %>%
  arrange(Gene_Category, desc(lncci_status)) %>%
  mutate(cumulative_percentage = cumsum(Percentage) - (Percentage / 2))

stacked_data_lncci$Gene_Category <- factor(stacked_data_lncci$Gene_Category,
                                           levels = c("No CHIP carrier", "DNMT3A", "TET2", "Other_genes"))

# Perform a trend test
# Create a contingency table of Gene_Category and SNCI status
contingency_table <- xtabs(Count ~ Gene_Category + lncci_status, data = stacked_data_lncci)
print(contingency_table)

chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)

p_value <- chi_square_result$p.value

# Plot the stacked bar chart
library(ggplot2)
library(ggsignif)

# Your plot with the significance bracket
ggplot(stacked_data_lncci, aes(x = Gene_Category, y = Percentage, fill = factor(lncci_status))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cumulative_percentage, label = paste0(round(Percentage, 1), "%")),
            color = "black", size = 3.5, vjust = 1.5) + 
  labs(x = "", y = "% of patients", fill = "LNCCI Presence") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_discrete(labels = c("No", "Yes")) +  # Use default colors and change labels
  theme_classic() +
  theme(legend.position = "top", 
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")) +
  
  # Adding the significance bracket
  geom_signif(comparisons = list(c("No CHIP carrier", "Other_genes")),
              map_signif_level = FALSE,
              y_position = 100,  # Place the bracket lower to avoid overlap
              tip_length = 0, textsize = 0) +
  
  # Adding p-value annotation above the bars
  annotate("text", x = 2.5, y = 105, label = "P=0.9013", size = 3.5, vjust = -0.5)

################################################################################
# Logistic regression models according to CHIP carrier status
# Univariable logistic regression
logistic_model <- glm(lncci_status ~ chip.status.carrier, data = dat_swiss_filtered, family = binomial)
summary(logistic_model)

dat_swiss_filtered$age_squared <- dat_swiss_filtered$age.bl^2
dat_swiss_filtered$age_standardized <- dat_swiss_filtered$age.bl / 10
dat_swiss_filtered$age_category <- ifelse(dat_swiss_filtered$age.bl >= 70, 1, 0)

library(forestmodel)
print(forest_model(glm(lncci_status ~ chip.status.carrier, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(lncci_status ~ chip.status.carrier + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(lncci_status ~ chip.status.carrier + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(lncci_status ~ chip.status.carrier + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to DNMT3A carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(lncci_status ~ DNMT3A, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(lncci_status ~ DNMT3A + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(lncci_status ~ DNMT3A + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(lncci_status ~ DNMT3A + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to TET2 carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(lncci_status ~ TET2, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(lncci_status ~ TET2 + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(lncci_status ~ TET2 + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(lncci_status ~ TET2 + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to other genes carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(lncci_status ~ Other_genes, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(lncci_status ~ Other_genes + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(lncci_status ~ Other_genes + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(lncci_status ~ Other_genes + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Data for the forest plot
labels <- c("No CHIP carrier (Reference)", 
            "CHIP carrier", 
            "DNMT3A", 
            "TET2", 
            "Other genes")

n_values <- c("1222/1562", 
              "340/1562", 
              "167/1562", 
              "96/1562", 
              "99/1562")

or_with_ci <- c("", 
                "0.96 (0.72 - 1.29)", 
                "0.91 (0.61 - 1.34)", 
                "0.82 (0.48 - 1.34)", 
                "1.10 (0.67 - 1.74)")

p_values <- c(NA, "0.81", "0.64", "0.44", "0.71")

odds_ratios <- c(NA, 0.96, 0.91, 0.82, 1.10)
lower_ci <- c(NA, 0.72, 0.61, 0.48, 0.67)
upper_ci <- c(NA, 1.29, 1.34, 1.34, 1.74)

# Add column headers
header <- c("", "N", "OR (95% CI)", "P value")

# Create a matrix for the labels and data
label_matrix <- rbind(header, cbind(labels, n_values, or_with_ci, p_values))

library(forestplot)

# Create the forest plot
forestplot(
  labeltext = label_matrix,  # Use the matrix for labels
  mean = c(NA, odds_ratios),
  lower = c(NA, lower_ci),
  upper = c(NA, upper_ci),
  is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Mark headers as summary
  xlog = TRUE,  # No log scale for OR
  zero = 1,  # Reference line for OR = 1
  xlab = "Odds Ratio (95% CI)",
  boxsize = 0.2,  # Size of boxes for OR
  ci.vertices = TRUE,  # Show CI vertices
  ci.vertices.height = 0.05,  # Height of CI vertices
  txt_gp = forestplot::fpTxtGp(
    label = gpar(fontsize = 12),          # Font size for labels
    ticks = gpar(fontsize = 16),          # Font size for x-axis tick labels
    xlab = gpar(fontsize = 14)  # Font size and bold for x-axis label
  )
)

# Association between number of CHIP mutations and SNI lesion
print(forest_model(glm(lncci_status ~ nb.chip.driver.genes + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Data for the forest plot
labels <- "Nr. of CHIP mutations"

n_values <- "1562"

or_with_ci <- "0.99 (0.79 - 1.23)"

p_values <- "0.95"

odds_ratios <- 0.99
lower_ci <- 0.79
upper_ci <- 1.23

# Add column headers
header <- c("", "N", "OR (95% CI)", "P value")

# Create a matrix for the labels and data
label_matrix <- rbind(header, cbind(labels, n_values, or_with_ci, p_values))

library(forestplot)

# Create the forest plot
forestplot(
  labeltext = label_matrix,  # Use the matrix for labels
  mean = c(NA, odds_ratios),
  lower = c(NA, lower_ci),
  upper = c(NA, upper_ci),
  is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Mark headers as summary
  xlog = TRUE,  # No log scale for OR
  zero = 1,  # Reference line for OR = 1
  xlab = "Odds Ratio (95% CI)",
  boxsize = 0.2,  # Size of boxes for OR
  ci.vertices = TRUE,  # Show CI vertices
  ci.vertices.height = 0.05,  # Height of CI vertices
  txt_gp = forestplot::fpTxtGp(
    label = gpar(fontsize = 12),          # Font size for labels
    ticks = gpar(fontsize = 16),          # Font size for x-axis tick labels
    xlab = gpar(fontsize = 14)            # Font size and bold for x-axis label
  ),
  xticks = seq(0.5, 2, by = 0.5)  # Define x-axis ticks from 0.5 to 2
)

################################################################################
# Association between CHIP status and microbleeds

dat_swiss_filtered$chip.status.carrier <- as.factor(dat_swiss_filtered$chip.status.carrier)

na_count_mb_count <- sum(is.na(dat_swiss_filtered$mb.nr))
# Exclude patients with snci.count==NA
dat_swiss_filtered <- dat_swiss_filtered[!is.na(dat_swiss_filtered$mb.nr), ]

# Generate Mb presence variable (1=yes, 0=no)
dat_swiss_filtered$mb_status <- ifelse(dat_swiss_filtered$mb.nr > 0, 1, 0)

dat_swiss_filtered$DNMT3A <- ifelse(grepl("DNMT3A", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$TET2 <- ifelse(grepl("TET2", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$Other_genes <- ifelse(dat_swiss_filtered$CHIP_carrier == 1 & !grepl("DNMT3A|TET2", dat_swiss_filtered$chip.genes), 1, 0)

# Create Gene Category: No CHIP carrier, CHIP carrier, DNMT3A, TET2, Other_genes
dat_swiss_filtered <- dat_swiss_filtered %>%
  mutate(Gene_Category = case_when(
    CHIP_carrier == 0 ~ "No CHIP carrier",
    DNMT3A == 1 ~ "DNMT3A",
    TET2 == 1 ~ "TET2",
    Other_genes == 1 ~ "Other_genes",
    CHIP_carrier == 1 ~ "CHIP carrier"
  ))

# Create a summary table for the stacked bar chart
stacked_data_mb <- dat_swiss_filtered %>%
  group_by(Gene_Category, mb_status) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Gene_Category) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Compute cumulative percentage for label positioning
stacked_data_mb <- stacked_data_mb %>%
  group_by(Gene_Category) %>%
  arrange(Gene_Category, desc(mb_status)) %>%
  mutate(cumulative_percentage = cumsum(Percentage) - (Percentage / 2))

stacked_data_mb$Gene_Category <- factor(stacked_data_mb$Gene_Category,
                                        levels = c("No CHIP carrier", "DNMT3A", "TET2", "Other_genes"))

# Perform a trend test
# Create a contingency table of Gene_Category and SNCI status
contingency_table <- xtabs(Count ~ Gene_Category + mb_status, data = stacked_data_mb)
print(contingency_table)

chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)

p_value <- chi_square_result$p.value

# Plot the stacked bar chart
library(ggplot2)
library(ggsignif)

# Your plot with the significance bracket
ggplot(stacked_data_mb, aes(x = Gene_Category, y = Percentage, fill = factor(mb_status))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cumulative_percentage, label = paste0(round(Percentage, 1), "%")), 
            color = "black", size = 3.5, vjust = 1.5) + 
  labs(x = "", y = "% of patients", fill = "Mb Presence") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_discrete(labels = c("No", "Yes")) +  # Use default colors and change labels
  theme_classic() +
  theme(legend.position = "top", 
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")) +
  
  # Adding the significance bracket
  geom_signif(comparisons = list(c("No CHIP carrier", "Other_genes")),
              map_signif_level = FALSE,
              y_position = 100,  # Place the bracket lower to avoid overlap
              tip_length = 0, textsize = 0) +
  
  # Adding p-value annotation above the bars
  annotate("text", x = 2.5, y = 105, label = "P=0.001545", size = 3.5, vjust = -0.5)

################################################################################
# Logistic regression models according to CHIP carrier status
# Univariable logistic regression
logistic_model <- glm(mb_status ~ chip.status.carrier, data = dat_swiss_filtered, family = binomial)
summary(logistic_model)

dat_swiss_filtered$age_squared <- dat_swiss_filtered$age.bl^2
dat_swiss_filtered$age_standardized <- dat_swiss_filtered$age.bl / 10
dat_swiss_filtered$age_category <- ifelse(dat_swiss_filtered$age.bl >= 70, 1, 0)

library(forestmodel)
print(forest_model(glm(mb_status ~ chip.status.carrier, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(mb_status ~ chip.status.carrier + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(mb_status ~ chip.status.carrier + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(mb_status ~ chip.status.carrier + age_category + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to DNMT3A carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(mb_status ~ DNMT3A, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(mb_status ~ DNMT3A + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(mb_status ~ DNMT3A + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking, bmi, hypertension and diabetes
print(forest_model(glm(mb_status ~ DNMT3A + age_category + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to TET2 carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(mb_status ~ TET2, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(mb_status ~ TET2 + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(mb_status ~ TET2 + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking, bmi, hypertension and diabetes
print(forest_model(glm(mb_status ~ TET2 + age_category + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to other genes carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(mb_status ~ Other_genes, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(mb_status ~ Other_genes + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(mb_status ~ Other_genes + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking, bmi, hypertension and diabetes
print(forest_model(glm(mb_status ~ Other_genes + age_category + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Data for the forest plot
labels <- c("No CHIP carrier (Reference)", 
            "CHIP carrier", 
            "DNMT3A", 
            "TET2", 
            "Other genes")

n_values <- c("1179/1516", 
              "333/1516", 
              "164/1516", 
              "91/1516", 
              "98/1516")

or_with_ci <- c("", 
                "1.45 (1.09 - 1.93)", 
                "1.03 (0.69 - 1.50)", 
                "1.50 (0.92 - 2.38)", 
                "1.70 (1.07 - 2.65)")

p_values <- c(NA, "0.010", "0.90", "0.09", "0.021")

odds_ratios <- c(NA, 1.45, 1.03, 1.50, 1.70)
lower_ci <- c(NA, 1.09, 0.69, 0.92, 1.07)
upper_ci <- c(NA, 1.93, 1.50, 2.38, 2.65)

# Add column headers
header <- c("", "N", "OR (95% CI)", "P value")

# Create a matrix for the labels and data
label_matrix <- rbind(header, cbind(labels, n_values, or_with_ci, p_values))

library(forestplot)

# Create the forest plot
forestplot(
  labeltext = label_matrix,  # Use the matrix for labels
  mean = c(NA, odds_ratios),
  lower = c(NA, lower_ci),
  upper = c(NA, upper_ci),
  is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Mark headers as summary
  xlog = TRUE,  # No log scale for OR
  zero = 1,  # Reference line for OR = 1
  xlab = "Odds Ratio (95% CI)",
  boxsize = 0.2,  # Size of boxes for OR
  ci.vertices = TRUE,  # Show CI vertices
  ci.vertices.height = 0.05,  # Height of CI vertices
  txt_gp = forestplot::fpTxtGp(
    label = gpar(fontsize = 12),          # Font size for labels
    ticks = gpar(fontsize = 16),          # Font size for x-axis tick labels
    xlab = gpar(fontsize = 14)  # Font size and bold for x-axis label
  )
)

# Association between number of CHIP mutations and SNI lesion
print(forest_model(glm(mb_status ~ nb.chip.driver.genes +  age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Data for the forest plot
labels <- "Nr. of CHIP mutations"

n_values <- "1516"

or_with_ci <- "1.36 (1.10 - 1.69)"

p_values <- "0.005"

odds_ratios <- 1.36
lower_ci <- 1.10
upper_ci <- 1.69

# Add column headers
header <- c("", "N", "OR (95% CI)", "P value")

# Create a matrix for the labels and data
label_matrix <- rbind(header, cbind(labels, n_values, or_with_ci, p_values))

library(forestplot)

# Create the forest plot
forestplot(
  labeltext = label_matrix,  # Use the matrix for labels
  mean = c(NA, odds_ratios),
  lower = c(NA, lower_ci),
  upper = c(NA, upper_ci),
  is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Mark headers as summary
  xlog = FALSE,  # No log scale for OR
  zero = 1,  # Reference line for OR = 1
  xlab = "Odds Ratio (95% CI)",
  boxsize = 0.2,  # Size of boxes for OR
  ci.vertices = TRUE,  # Show CI vertices
  ci.vertices.height = 0.05,  # Height of CI vertices
  txt_gp = forestplot::fpTxtGp(
    label = gpar(fontsize = 12),          # Font size for labels
    ticks = gpar(fontsize = 16),          # Font size for x-axis tick labels
    xlab = gpar(fontsize = 14)            # Font size and bold for x-axis label
  ),
  xticks = seq(0.5, 2, by = 0.5)  # Define x-axis ticks from 0.5 to 2
)

################################################################################
# Association between CHIP status and white matter lesions

dat_swiss_filtered$CHIP_carrier <- as.factor(dat_swiss_filtered$chip.status.carrier)

na_count_mb_count <- sum(is.na(dat_swiss_filtered$mod.faz))
# Exclude patients with snci.count==NA
dat_swiss_filtered <- dat_swiss_filtered[!is.na(dat_swiss_filtered$mod.faz), ]

dat_swiss_filtered$DNMT3A <- ifelse(grepl("DNMT3A", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$TET2 <- ifelse(grepl("TET2", dat_swiss_filtered$chip.genes), 1, 0)
dat_swiss_filtered$Other_genes <- ifelse(dat_swiss_filtered$CHIP_carrier == 1 & !grepl("DNMT3A|TET2", dat_swiss_filtered$chip.genes), 1, 0)

# Generate WML presence variable according to the mod. Fazekas >2 (1=yes, 0=no)
dat_swiss_filtered$wml_status <- ifelse(dat_swiss_filtered$mod.faz == "Yes", 1, 0)

# Create Gene Category: No CHIP carrier, CHIP carrier, DNMT3A, TET2, Other_genes
dat_swiss_filtered <- dat_swiss_filtered %>%
  mutate(Gene_Category = case_when(
    CHIP_carrier == 0 ~ "No CHIP carrier",
    DNMT3A == 1 ~ "DNMT3A",
    TET2 == 1 ~ "TET2",
    Other_genes == 1 ~ "Other_genes",
    CHIP_carrier == 1 ~ "CHIP carrier"
  ))

# Create a summary table for the stacked bar chart
stacked_data_wml <- dat_swiss_filtered %>%
  group_by(Gene_Category, wml_status) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Gene_Category) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Compute cumulative percentage for label positioning
stacked_data_wml <- stacked_data_wml %>%
  group_by(Gene_Category) %>%
  arrange(Gene_Category, desc(wml_status)) %>%
  mutate(cumulative_percentage = cumsum(Percentage) - (Percentage / 2))

stacked_data_wml$Gene_Category <- factor(stacked_data_wml$Gene_Category,
                                         levels = c("No CHIP carrier", "DNMT3A", "TET2", "Other_genes"))

# Perform a trend test
# Create a contingency table of Gene_Category and SNCI status
contingency_table <- xtabs(Count ~ Gene_Category + wml_status, data = stacked_data_wml)
print(contingency_table)

chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)

# Plot the stacked bar chart
library(ggplot2)
library(ggsignif)

# Your plot with the significance bracket
ggplot(stacked_data_wml, aes(x = Gene_Category, y = Percentage, fill = factor(wml_status))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cumulative_percentage, label = paste0(round(Percentage, 1), "%")), 
            color = "black", size = 3.5, vjust = 1.5) + 
  labs(x = "", y = "% of patients", fill = "WML Presence (mod. Fazekas >2)") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_discrete(labels = c("No", "Yes")) +
  theme_classic() +
  theme(legend.position = "top", 
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")) +
  
  # Adding the significance bracket
  geom_signif(comparisons = list(c("No CHIP carrier", "Other_genes")),
              map_signif_level = FALSE,
              y_position = 100,  # Place the bracket lower to avoid overlap
              tip_length = 0, textsize = 0) +
  
  # Adding p-value annotation above the bars
  annotate("text", x = 2.5, y = 105, label = "P=1.857e-06", size = 3.5, vjust = -0.5)

################################################################################
# Logistic regression models according to CHIP carrier status
# Univariable logistic regression
logistic_model <- glm(wml_status ~ chip.status.carrier, data = dat_swiss_filtered, family = binomial)
summary(logistic_model)

dat_swiss_filtered$age_squared <- dat_swiss_filtered$age.bl^2
dat_swiss_filtered$age_standardized <- dat_swiss_filtered$age.bl / 10
dat_swiss_filtered$age_category <- ifelse(dat_swiss_filtered$age.bl >= 70, 1, 0)

library(forestmodel)
print(forest_model(glm(wml_status ~ chip.status.carrier, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(wml_status ~ chip.status.carrier + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(wml_status ~ chip.status.carrier + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(wml_status ~ chip.status.carrier + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to DNMT3A carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(wml_status ~ DNMT3A, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(wml_status ~ DNMT3A + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(wml_status ~ DNMT3A + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(wml_status ~ DNMT3A + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to TET2 carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(wml_status ~ TET2, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(wml_status ~ TET2 + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(wml_status ~ TET2 + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(wml_status ~ TET2 + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Logistic regression models according to other genes carrier status
# Univariable logistic regression
library(forestmodel)
print(forest_model(glm(wml_status ~ Other_genes, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age and sex
print(forest_model(glm(wml_status ~ Other_genes + age.bl + pat.sex, data = dat_swiss_filtered, family = binomial)))

# Adjusted for age, sex, current smoking and bmi
print(forest_model(glm(wml_status ~ Other_genes + age.bl + pat.sex + current.smoker + 
                         bmi, data = dat_swiss_filtered, family = binomial)))

# Multivariable adjusted
print(forest_model(glm(wml_status ~ Other_genes + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Data for the forest plot
labels <- c("No CHIP carrier (Reference)", 
            "CHIP carrier", 
            "DNMT3A", 
            "TET2", 
            "Other genes")

n_values <- c("1225/1566", 
              "341/1566", 
              "168/1566", 
              "97/1566", 
              "99/1566")

or_with_ci <- c("", 
                "1.56 (1.20 - 2.04)", 
                "1.23 (0.87 - 1.75)", 
                "1.33 (0.84 - 2.13)", 
                "1.70 (1.07 - 2.74)")

p_values <- c(NA, "0.001", "0.25", "0.23", "0.03")

odds_ratios <- c(NA, 1.56, 1.23, 1.33, 1.70)
lower_ci <- c(NA, 1.20, 0.87, 0.84, 1.07)
upper_ci <- c(NA, 2.04, 1.75, 2.13, 2.74)

# Add column headers
header <- c("", "N", "OR (95% CI)", "P value")

# Create a matrix for the labels and data
label_matrix <- rbind(header, cbind(labels, n_values, or_with_ci, p_values))

library(forestplot)

# Create the forest plot
forestplot(
  labeltext = label_matrix,  # Use the matrix for labels
  mean = c(NA, odds_ratios),
  lower = c(NA, lower_ci),
  upper = c(NA, upper_ci),
  is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Mark headers as summary
  xlog = TRUE,  # No log scale for OR
  zero = 1,  # Reference line for OR = 1
  xlab = "Odds Ratio (95% CI)",
  boxsize = 0.2,  # Size of boxes for OR
  ci.vertices = TRUE,  # Show CI vertices
  ci.vertices.height = 0.05,  # Height of CI vertices
  txt_gp = forestplot::fpTxtGp(
    label = gpar(fontsize = 12),          # Font size for labels
    ticks = gpar(fontsize = 16),          # Font size for x-axis tick labels
    xlab = gpar(fontsize = 14)  # Font size and bold for x-axis label
  )
)

# Association between number of CHIP mutations and SNI lesion
print(forest_model(glm(wml_status ~ nb.chip.driver.genes + age_category + pat.sex + bmi +
                         current.smoker + prev.hypertonie + med.oak.yn + 
                         krk.malignom, data = dat_swiss_filtered, family = binomial)))

# Data for the forest plot
labels <- "Nr. of CHIP mutations"

n_values <- "1566"

or_with_ci <- "1.41 (1.15 - 1.75)"

p_values <- "0.001"

odds_ratios <- 1.41
lower_ci <- 1.15
upper_ci <- 1.75

# Add column headers
header <- c("", "N", "OR (95% CI)", "P value")

# Create a matrix for the labels and data
label_matrix <- rbind(header, cbind(labels, n_values, or_with_ci, p_values))

library(forestplot)

# Create the forest plot
forestplot(
  labeltext = label_matrix,  # Use the matrix for labels
  mean = c(NA, odds_ratios),
  lower = c(NA, lower_ci),
  upper = c(NA, upper_ci),
  is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),  # Mark headers as summary
  xlog = TRUE,  # No log scale for OR
  zero = 1,  # Reference line for OR = 1
  xlab = "Odds Ratio (95% CI)",
  boxsize = 0.2,  # Size of boxes for OR
  ci.vertices = TRUE,  # Show CI vertices
  ci.vertices.height = 0.05,  # Height of CI vertices
  txt_gp = forestplot::fpTxtGp(
    label = gpar(fontsize = 12),          # Font size for labels
    ticks = gpar(fontsize = 16),          # Font size for x-axis tick labels
    xlab = gpar(fontsize = 14)            # Font size and bold for x-axis label
  ),
  xticks = seq(0.5, 2, by = 0.5)  # Define x-axis ticks from 0.5 to 2
)

################################################################################
# Association between CHIP status and history of stroke

# Create Gene Category: No CHIP carrier, CHIP carrier, DNMT3A, TET2, Other_genes
dat_swiss_filtered <- dat_swiss_filtered %>%
  mutate(Gene_Category = case_when(
    CHIP_carrier == 0 ~ "No CHIP carrier",
    DNMT3A == 1 ~ "DNMT3A",
    TET2 == 1 ~ "TET2",
    Other_genes == 1 ~ "Other_genes",
    CHIP_carrier == 1 ~ "CHIP carrier"
  ))

# Create a summary table for the stacked bar chart
stacked_data_stroke <- dat_swiss_filtered %>%
  group_by(Gene_Category, prev.stroke.tia) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Gene_Category) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Compute cumulative percentage for label positioning
stacked_data_stroke <- stacked_data_stroke %>%
  group_by(Gene_Category) %>%
  arrange(Gene_Category, desc(prev.stroke.tia)) %>%
  mutate(cumulative_percentage = cumsum(Percentage) - (Percentage / 2))

stacked_data_stroke$Gene_Category <- factor(stacked_data_stroke$Gene_Category,
                                          levels = c("No CHIP carrier", "DNMT3A", "TET2", "Other_genes"))

# Perform a trend test
# Create a contingency table of Gene_Category and SNCI status
contingency_table <- xtabs(Count ~ Gene_Category + prev.stroke.tia, data = stacked_data_stroke)
print(contingency_table)

chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)

p_value <- chi_square_result$p.value

# Plot the stacked bar chart
library(ggplot2)
library(ggsignif)

# Your plot with the significance bracket
ggplot(stacked_data_stroke, aes(x = Gene_Category, y = Percentage, fill = factor(prev.stroke.tia))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cumulative_percentage, label = paste0(round(Percentage, 1), "%")), 
            color = "black", size = 3.5, vjust = 1.5) + 
  labs(x = "", y = "% of patients", fill = "Stroke/TIA Prevalent") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_discrete(labels = c("No", "Yes")) +  # Use default colors and change labels
  theme_classic() +
  theme(legend.position = "top", 
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")) +
  
  # Adding the significance bracket
  geom_signif(comparisons = list(c("No CHIP carrier", "Other_genes")),
              map_signif_level = FALSE,
              y_position = 100,  # Place the bracket lower to avoid overlap
              tip_length = 0, textsize = 0) +
  
  # Adding p-value annotation above the bars
  annotate("text", x = 2.5, y = 105, label = "P=0.007368", size = 3.5, vjust = -0.5)

################################################################################
# Association between CHIP status and new ischemic stroke

# Create Gene Category: No CHIP carrier, CHIP carrier, DNMT3A, TET2, Other_genes
dat_swiss_filtered <- dat_swiss_filtered %>%
  mutate(Gene_Category = case_when(
    CHIP_carrier == 0 ~ "No CHIP carrier",
    DNMT3A == 1 ~ "DNMT3A",
    TET2 == 1 ~ "TET2",
    Other_genes == 1 ~ "Other_genes",
    CHIP_carrier == 1 ~ "CHIP carrier"
  ))

# Create a summary table for the stacked bar chart
stacked_data_newstroke <- dat_swiss_filtered %>%
  group_by(Gene_Category, stroke.tia) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Gene_Category) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Compute cumulative percentage for label positioning
stacked_data_newstroke <- stacked_data_newstroke %>%
  group_by(Gene_Category) %>%
  arrange(Gene_Category, desc(stroke.tia)) %>%
  mutate(cumulative_percentage = cumsum(Percentage) - (Percentage / 2))

stacked_data_newstroke$Gene_Category <- factor(stacked_data_newstroke$Gene_Category,
                                            levels = c("No CHIP carrier", "DNMT3A", "TET2", "Other_genes"))

# Perform a trend test
# Create a contingency table of Gene_Category and SNCI status
contingency_table <- xtabs(Count ~ Gene_Category + stroke.tia, data = stacked_data_newstroke)
print(contingency_table)

chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)

p_value <- chi_square_result$p.value

# Plot the stacked bar chart
library(ggplot2)
library(ggsignif)

# Your plot with the significance bracket
ggplot(stacked_data_newstroke, aes(x = Gene_Category, y = Percentage, fill = factor(stroke.tia))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = cumulative_percentage, label = paste0(round(Percentage, 1), "%")), 
            color = "black", size = 3.5, vjust = 1.5) + 
  labs(x = "", y = "% of patients", fill = "New stroke/TIA") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_discrete(labels = c("No", "Yes")) +  # Use default colors and change labels
  theme_classic() +
  theme(legend.position = "top", 
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black")) +
  
  # Adding the significance bracket
  geom_signif(comparisons = list(c("No CHIP carrier", "Other_genes")),
              map_signif_level = FALSE,
              y_position = 100,  # Place the bracket lower to avoid overlap
              tip_length = 0, textsize = 0) +
  
  # Adding p-value annotation above the bars
  annotate("text", x = 2.5, y = 105, label = "P=0.8528", size = 3.5, vjust = -0.5)

################################# END ##########################################
