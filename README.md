# histoplasmosis

# Histoplasmosis Analysis

## Overview

This repository contains scripts and models for analyzing histoplasmosis, a fungal infection caused by *Histoplasma capsulatum*. The tools provided here are designed to facilitate sequential testing and treatment modeling for this disease.

## Contents

- `histo_seq_testing.R`: An R script for conducting sequential testing analysis related to histoplasmosis.
- `histo_sequential_testing.html`: An HTML file to simulate sequential testing results.
- `histo_treat_Markov_model.R`: An R script for modeling treatment strategies for histoplasmosis using a Markov model approach.
- `model_seq_pos_retested_binom.stan`: A Stan model file used for statistical analysis in the sequential testing framework.

## Prerequisites

Before using the scripts in this repository, ensure you have the following software installed:

- [R](https://www.r-project.org/): A programming language for statistical computing.
- [RStudio](https://posit.co/download/rstudio/): An integrated development environment for R.
- [Stan](https://mc-stan.org/): A platform for statistical modeling and high-performance statistical computation.

## Installation

1. **Clone the repository**:

   ```bash
   git clone https://github.com/antoniofl0res/histoplasmosis.git
   cd histoplasmosis
   ```

2. **Install required R packages**:

   Open R or RStudio and run:

   ```R
   install.packages(c("rstan", "ggplot2", "dplyr"))
   ```

## Usage

1. **Sequential Testing Analysis**:

   To perform sequential testing analysis, run the `histo_seq_testing.R` script:

   ```bash
   Rscript histo_seq_testing.R
   ```

   This script will generate outputs detailing the sequential testing procedures and results.

2. **Treatment Modeling**:

   For treatment modeling, execute the `histo_treat_model.R` script:

   ```bash
   Rscript histo_treat_model.R
   ```

   This will run the treatment models and provide insights into various treatment strategies.



*Note: This repository is intended for research and educational purposes. Ensure compliance with local regulations and guidelines when handling and analyzing data related to histoplasmosis.* 
