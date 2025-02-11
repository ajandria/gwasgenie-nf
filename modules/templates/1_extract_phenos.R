#!/usr/bin/env Rscript

library(readxl)
library(readr)
library(dplyr)

# Load the data

# Define file path
gwas_sheet_path <- "${gwas_sheet}"

# Check file extension and read accordingly
if (grepl("xlsx", gwas_sheet_path)) {
  excel_gwas_sheet=TRUE
  gwas_sheet <- read_excel(gwas_sheet_path) # Read Excel file
  gwas_sheet <- gwas_sheet %>% 
    rename(phenotype = `Trait (on original file)`,
           covariates = `Covariates needed`) %>%
    mutate(covariates = toupper(gsub(' ', '', covariates)),
           phenotype = toupper(phenotype),
           `Results folder` = gsub(' ','_',`Results folder`)) %>%
    filter(`Test association yes/no` == 'Y')

  pcs <- readr::read_delim("${pcs_path}", col_names = FALSE) %>% 
    rename(FAM = X1,
          `Study ID` = X2)
  colnames(pcs)[-c(1,2)] <- c('PC1', 'PC2', 'PC3', 'PC4',
                              'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')

  phenos <- read_excel("${phenos}")  %>% # Clinical data with sample ID and measurements
    left_join(pcs) %>%
    mutate(sex = ifelse(Gender == 'Female',1,ifelse(Gender == 'Male', 0, Gender)))
} else {
  gwas_sheet <- read_tsv(gwas_sheet_path) # Read TSV file
  phenos <- read_excel("${phenos}") %>% # Clinical data with sample ID and measurements
    mutate(sex = ifelse(Gender == 'Female',1,ifelse(Gender == 'Male', 0, Gender))) 
}

colnames(phenos) <- toupper(colnames(phenos))
print(gwas_sheet)
print(phenos)

# Iterate through each row of the gwas_sheet
apply(gwas_sheet, 1, function(row) {
  # Extract the phenotype and covariates
  pheno <- as.character(row["phenotype"]) # Extract phenotype name
  covariates <- unlist(strsplit(row["covariates"], ",")) # Split covariates into separate columns

  # Check if the phenotype exists in phenos
  #if (!(pheno %in% colnames(phenos))) {
  #  stop(paste("Error: Phenotype", pheno, "does not exist in the dataset `phenos`"))
  #}
  if (!(pheno %in% colnames(phenos))) {
    message(paste("Warning: Phenotype", pheno, "does not exist in the dataset `phenos`. Skipping..."))
    return(NULL)  # Skip to the next iteration
  }

  # Check if all covariates exist in phenos
  missing_covariates <- setdiff(covariates, colnames(phenos))
  #if (length(missing_covariates) > 0) {
  #  stop(paste("Error: Covariates", paste(missing_covariates, collapse = ", "), "do not exist in the dataset `phenos`"))
  #}
  if (length(missing_covariates) > 0) {
    message(paste("Warning: Covariates", paste(missing_covariates, collapse = ", "), "do not exist in the dataset `phenos`. Continuing with available data..."))
    covariates <- setdiff(covariates, missing_covariates)  # Remove missing covariates and continue
  }
  
  # Filter for phenotype
  pheno_cov_data <- phenos %>%
    select(FID, IID, all_of(pheno), all_of(covariates)) %>% # Use `SAMPLE` as ID and the specified phenotype column
    na.omit()
  colnames(pheno_cov_data)[-c(1,2,3)] <- covariates

  # Check if there are no rows left after filtering
  if (nrow(pheno_cov_data) == 0) {
    message(paste("Warning: No data available for phenotype", pheno, "after filtering. Skipping..."))
    return(NULL)  # Skip this iteration
  }

  pheno_data <- pheno_cov_data %>%
    select(FID, IID, all_of(pheno))

  if (excel_gwas_sheet == TRUE) {
    if (row["Transformation needed"] == "RINT") {
      pheno_data[[pheno]] <- qnorm((rank(pheno_data[[pheno]], na.last="keep") - 0.5) / sum(!is.na(pheno_data[[pheno]])))
    }
  }

  # Generate filenames based on phenotype and covariates
  covariates_name <- paste(covariates, collapse = "_")
  pheno_file <- paste0(gsub(' ','_',pheno), "-", covariates_name, "_pheno.txt")
  covariate_file <- paste0(gsub(' ','_',pheno), "-", covariates_name, "_covs.txt")

  # Write phenotype file
  colnames(pheno_data) <- gsub(' ','_',colnames(pheno_data))
  write_tsv(pheno_data, pheno_file)

  # Filter for covariates
  covariate_data <- pheno_cov_data %>%
    select(FID, IID, all_of(covariates)) # Use `SAMPLE` as ID and the specified covariates

  # Rename columns to retain their original names
  colnames(covariate_data) <- gsub(' ','_',colnames(covariate_data))
  colnames(covariate_data)[-c(1,2)] <- covariates

  # Write covariate file
  write_tsv(covariate_data, covariate_file)
})
