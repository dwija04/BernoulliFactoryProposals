# BernoulliFactoryProposals
This repository contains the code for all the examples give in the paper titled "Exact MCMC Sampling for Intractable Proposals" by Dwija Kakkad and Dootika Vats. There are three examples in the paper:

- Truncated Gaussian
- Sensor Network Localization
- Gaussian Process modulated Cox process

There is a folder corresponding to each example. Every example has three main files:
- `***_functions.R` : contains all the background functions required for that example
- `***_single_run.R` : contains code for implementing a single chain to make any plots. The output is an Rdata file of the chain.
- `***_run.R` : contains parallelized replications for ESS, computational time, and other comparisons of different relevant samplers. The output is an Rdata file, organized like an incepted list.
- `***_output.R` :  contains code to create all plots, tables etc. This is completely reprodcuible.

The plots generated are given in a 'plots' folder. The .RData file for a single chain in the Cox modulated GP example is not provided due to its large size and can be obtained by running the 'cox_single_run' file (will take roughly 35 hours).
