# The purpose of this script is to generate TPM values for the entire delta14
# sequencing library.

# Clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)

# Initialize a new function to create directories
newDir <- function(x) {
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

################################################################################

# Read in the file containing mmu gene lengths
lengths <- read.csv("Data_Files/Mmu_Gene_Lengths/Mus_musculus_exon_lengths_in_BP.csv") %>%
  column_to_rownames(var = "X")

# Convert the gene lengths from bp to kbp
lengths$exon.lengths <- lengths$exon.lengths / 1000
lengths$Gene <- NULL

# Read in the full, non-normalized counts matrix
counts <- read.csv("Data_Files/Master_D14_Counts.csv") %>%
  column_to_rownames(var = "X")

# Save the formatted gene names in a vector for later
genesFormatted <- rownames(counts)

# Convert the formatted gene names to Ensembl
rownames(counts) <- gsub("^(.*?)\\s-.*", "\\1", rownames(counts))

# Function to generate TPM values
countToTPM <- function(count_df, length_df) {
  
  # Only take the genes in the length df that we have in our expression mat
  genesKeep <- rownames(count_df)
  
  # Subset the lengths
  length_df <- length_df[rownames(length_df) %in% genesKeep,]
  
  # Scale each count by its respective gene length
  rate <- count_df/length_df
  
  return(rate)
  # Initialize an empty vector
  tpm <- c()
  
  # Iterate through the counts and scale by library size
  for (i in 1:nrow(count_df)) {
    tpm[i] <- rate[i]/sum(rate) * 1e6
  }
  
  return(tpm)
}

# Function to generate log2(tpm) values
log2TPMs <- function(TPMs) {
  
  # Add a small constant
  TPMs <- TPMs + 0.01
  
  # Log transform the tpm values
  logTPMs <- as.data.frame(t(apply(TPMs, 1, log2)))
  
  return(logTPMs)
}

################################################################################

# Convert raw, non-normalized counts to TPM values
tpm_vals <- countToTPM(counts, lengths)
rownames(tpm_vals) <- genesFormatted

# Save tpm_values as a csv
write.csv(tpm_vals, file = "Data_Files/TPM_Values/Non_Transformed_TPM_Values_Delta14_Project.csv")

# Log2 transform the tpm values
logTPMs <- log2TPMs(tpm_vals)
rownames(logTPMs) <- genesFormatted

# Write the log2 transformed TPM values as a csv (typically would use these for heatmap generation)
write.csv(logTPMs, file = "Data_Files/TPM_Values/Log2_Transformed_TPM_Values_Delta14_project.csv")


