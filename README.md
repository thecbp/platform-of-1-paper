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

The code is designed to output the results of each individual simulation into a pre-specified directory. The simulations can be compiled together into a single tibble using the `compileSimulations()` function in the `helpers.R` script. 

```
source("helpers.R")

# assuming that the simulations are stored in "simulations/test"
out = compileSimulations("simulations/test")
```

We recommend using the scripts if you plan to run > 1000 simulations. This takes advantage of parallelization and speeds up the process.

# How to use: Shiny dashboard

We have provided a Shiny dashboard for running small scale simulations. The R code for the app is located in the `app` directory. After downloading this repo to your local machine, you can run the app with the following command:

```
library(shiny)
runApp("app")
```