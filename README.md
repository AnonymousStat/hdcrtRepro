# hdcrtRepro

`hdcrtRepro` is an R package containing utility functions and benchmark implementations used in the reproducibility materials for the paper:

**Hypothesis testing in high-dimensional censored-transformation models**

This package is prepared for anonymous review and reproducibility evaluation. It provides functions for data generation, simulation setting construction, result summarization, plotting, benchmark tests, and real data analysis utilities.

## Overview

The package supports the reproducibility workflow for simulation studies and empirical analysis in the paper. In particular, it includes functions for:

- generating covariates, regression coefficients, and random errors;
- generating survival outcomes under different censoring mechanisms;
- controlling censoring constants to achieve target censoring rates;
- summarizing empirical Type-I error rates and empirical powers;
- converting simulation result arrays into data frames;
- producing empirical power plots and size tables;
- implementing benchmark tests used in simulation comparisons;
- supporting repeated sample-splitting and pathway-level analysis for the TCGA SKCM data.

## Installation

This package is intended to be used together with the submitted reproducibility materials.

From the package root directory, install it in R by:

```r
devtools::install_github("AnonymousStat/hdcrtRepro")
```

## Requirements

This package requires R and a working C compiler, since some benchmark methods rely on C routines.

Required R packages are listed in the DESCRIPTION file. The main dependencies include packages for the proposed testing procedures, matrix computation, simulation, plotting, and benchmark comparisons.

## Use with reproducibility materials

This package is designed to be used together with the submitted reproducibility materials.

The full workflow is documented in the main README.md file of the reproducibility package. In general:

- Run the corresponding simulation scripts to generate p-value matrices.
- Run the corresponding summary scripts to compute empirical Type-I error rates and empirical powers.
- Run the real data analysis script to reproduce the TCGA SKCM analysis. 
  
## Notes

This is an anonymous version prepared for review. Author-identifying information has been removed where appropriate. After the review process, the package may be updated with author information and maintained in a public repository.

The package is provided for academic and research reproducibility purposes.










