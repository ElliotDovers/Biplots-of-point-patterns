# Biplots-of-point-patterns
This repo accompanies the manuscript "Biplots of multivariate spatial point processes - an exploratory tool to visualise structure in multiple, correlated point patterns" submitted for review.
The code can be used to replicate both the simulations and real data analyses therein.

## Lansing Woods
Contains code to analyze the Lansing Woods dataset from the spatstat package on R. Including model selection, fitting, and constructing the biplot. Main script is: `analysis.R`

## Northern NSW
Conatains code to analyse the northern New South Wales dataset from the disdat package on R. Includes methods of model selection, fitting, and construction of the biplot. For main model fitting see `mvlgcp_flora_iterated.R` which is a parallelised script for running on a HPC (set job = 1 for an example). For the main analysis, including created biplots from the fitted model, please see `model_explore.R`.

## Simulations
Contains code for simulations to investigate a data-driven strategy for choosing the number of basis functions used within the multivariate log-Gaussian Cox process underpinning the proposed biplot procedure.
