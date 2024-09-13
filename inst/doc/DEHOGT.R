## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----installation, eval=FALSE-------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE)) {
#      install.packages("BiocManager")
#  }
#  
#  BiocManager::install("DEHOGT")

## ----exampleWorkflow----------------------------------------------------------
## Simulate gene expression data (100 genes, 10 samples)
data <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)

## Randomly assign treatment groups
treatment <- sample(0:1, 10, replace = TRUE)
## Load DEHOGT package
library(DEHOGT)

## Run the function with 2 CPU cores
result <- dehogt_func(data, treatment, num_cores = 2)

## Display results
head(result$pvals)

# Example: Adding covariates and normalization factors
covariates <- matrix(rnorm(1000), nrow = 100, ncol = 10)
norm_factors <- rep(1, 10)

# Run with covariates and normalization factors
result_cov <- dehogt_func(data, treatment, covariates = covariates, norm_factors = norm_factors, num_cores = 2)

## ----sessionInfo, results='asis'----------------------------------------------
sessionInfo()

