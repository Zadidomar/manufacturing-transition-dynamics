# manufacturing-transition-dynamics
This project estimates unknown parameters of the assumed densities (Gamma and Exponential) of the time needed for the production run and delivery of a product and the unknown point of transition at which data distribution is changed.

## Methods
- Monte Carlo Markov Chain (MCMC)
- Bayesian analysis using prior and posterior densities
- Comparison of accuracy and efficiency between Single Component and Independent Metropolis-Hastings (MH) algorithms
- Synthetic data

## Tools
- R, Rmarkdown

## Key Results
- Single component MH always performed better than independent MH in estimating the point of transition. 
- Independent MH was always much faster than Single component MH in all cases. 
- Single component MH was a little more accurate than independent MH in estimating other parameters.
