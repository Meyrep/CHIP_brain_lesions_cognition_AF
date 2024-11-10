#--CHIP ANALYSIS----------------------------------------------------------------
#--CHIP, brain lesions and cognitive performance in AF--------------------------
# Author: Pascal B. Meyre
# Date: 10/21/24
# Location: Reinach, Baselland, Switzerland

# Use combined dataset of SWISS-and BEAT-AF pts.
# Use CHIP dataset of SWISS-and BEAT-AF pts.
# Use the brain lesions follow-up file
#-------------------------------------------------------------------------------

################################################################################
# Merge SWISS/BEAT-AF and cognitive performance datasets
################################################################################

# Import SWISS/BEAT-AF dataset
library(data.table)
swiss.beat <- fread("/Users/pascalmeyre/Desktop/Research/1_Projects_Analysis/19_CHIP_brain_lesions_cognition/analysis/datasets/chip.20241003.csv")

chip <- swiss.beat %>%
  filter(study == "Swiss-AF") %>%
  select(pat.id, study, chip.status.carrier, nb.chip.driver.genes, chip.genes)

# Write the filtered data to a file (e.g., CSV)
write.csv(filtered_data, file = "filtered_swiss_af_patients.csv", row.names = FALSE)

# Make a dataset only including CHIP variable and pat_id
library(dplyr)

chip <- swiss.beat %>%
  select(pat.id, chip.status.carrier, nb.chip.driver.genes, chip.genes)

# Import cognition dataset
library(data.table)
cognition <- fread("/Users/pascalmeyre/Desktop/Research/1_Projects_Analysis/19_CHIP_brain_lesions_cognition/analysis/datasets/saf.cog.20241104.csv")

# Make a dataset only including cognitive variables and pat_id
library(dplyr)

cognitive <- cognition %>%
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

