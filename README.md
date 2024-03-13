This repo is associated with the manuscript, "Platform-of-1: A Bayesian Adaptive N-of-1 Trial Design For Personalizing Treatment Among Multiple Candidates". It contains code for simulating data from a Platform-of-1 Trial and an accompanying Shiny dashboard for running simulations. 

# How to use: Scripts

There are 4 main files to consider:

- `simulation-pipeline.R`: Main interface for users. Users can specify their their desired parameters for the Platform-of-1 design. This interface creates a grid of all the simulation repetitions and keeps track of which have been simulated. Uses furr library for parallelizing simulations. 
- `simulatePlatformOf1.R`: contains function for running Platform-of-1 trials without interim decisions (for examining how operating characteristics evolve over course of trial)
- `simulatePlatformOf1Decision.R`: contains function for running Platform-of-1 trials with interim decisions
- `helpers.R`: contains functions for Thompson Sampling and probability rebalancing algorithm

To use `simulation-pipeline.R`, a user will need to do the following: 

1. Where to place simulation results
2. Define the parameters for the Platform-of-1 trial
3. Define the parameters for Thompson Sampling, MCMC and interim decisions
4. Run the simulations

# How to use: Shiny dashboard