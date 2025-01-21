#!/usr/bin/env Rscript

library(readxl)
library(readr)
library(dplyr)

# Load the data
gwas_sheet <- read_tsv("${gwas_sheet}") # GWAS sheet with sample, phenotype, and covariates
phenos <- read_excel("${phenos}") # Clinical data with sample ID and measurements

print(gwas_sheet)
print(phenos)

# Iterate through each row of the gwas_sheet
apply(gwas_sheet, 1, function(row) {
  # Extract the phenotype and covariates
  pheno <- as.character(row["phenotype"]) # Extract phenotype name
  covariates <- unlist(strsplit(row["covariates"], ",")) # Split covariates into separate columns

  # Check if the phenotype exists in phenos
  if (!(pheno %in% colnames(phenos))) {
    stop(paste("Error: Phenotype", pheno, "does not exist in the dataset `phenos`"))
  }

  # Check if all covariates exist in phenos
  missing_covariates <- setdiff(covariates, colnames(phenos))
  if (length(missing_covariates) > 0) {
    stop(paste("Error: Covariates", paste(missing_covariates, collapse = ", "), "do not exist in the dataset `phenos`"))
  }

  # Filter for phenotype
  pheno_cov_data <- phenos %>%
    select(FID, IID, all_of(pheno), all_of(covariates)) %>% # Use `SAMPLE` as ID and the specified phenotype column
    na.omit()
  colnames(pheno_cov_data)[-c(1,2,3)] <- covariates

  pheno_data <- pheno_cov_data %>%
    select(FID, IID, all_of(pheno))

  # Generate filenames based on phenotype and covariates
  covariates_name <- paste(covariates, collapse = "_")
  pheno_file <- paste0(pheno, "-", covariates_name, "_pheno.txt")
  covariate_file <- paste0(pheno, "-", covariates_name, "_covs.txt")

  # Write phenotype file
  write_tsv(pheno_data, pheno_file)

  # Filter for covariates
  covariate_data <- pheno_cov_data %>%
    select(FID, IID, all_of(covariates)) # Use `SAMPLE` as ID and the specified covariates

  # Rename columns to retain their original names
  colnames(covariate_data)[-c(1,2)] <- covariates

  # Write covariate file
  write_tsv(covariate_data, covariate_file)
})
