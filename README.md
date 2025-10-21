# Hy-MLC: A Hybrid Identification Workflow using Machine Learning and ABC

**Hy-MLC** is a novel bioinformatics workflow for identifying hybrid generations (e.g., F1, F2, Backcross) from genome-wide SNP data. It overcomes limitations of traditional methods like NewHybrids by combining a Bayesian hybrid index with Machine Learning (Random Forest) and Approximate Bayesian Computation (ABC) to handle thousands of loci efficiently.

> **Reference:** If you use Hy-MLC in your research, please cite: [Your Publication Here]

## Overview

Traditional tools for hybrid classification struggle with large SNP datasets. Hy-MLC addresses this with a four-step process:
1.  **Simulation** of hybrid genotypes from pure parental samples.
2.  **Hybrid Index** calculation for all individuals.
3.  **Feature Selection** using a Random Forest classifier to identify the most informative summary statistics.
4.  **Final Classification** using an Approximate Bayesian Computation (ABC) model for probabilistic assignment.

This workflow is implemented in R.

## Prerequisites

Before you begin, ensure you have the following R packages installed:

```r
# Install required packages
install.packages(c("adegenet", "gghybrid", "randomForest", "abc"))
```

