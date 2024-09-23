This repo is associated with the manuscript, "Platform-of-1: a Bayesian adaptive N-of-1 trial design for identifying an optimal treatment".

In this repo, you can find:

- an accompanying Shiny dashboard for running simulations. 
- code for simulating data from a Platform-of-1 Trial

# How to use: Shiny dashboard

We have provided a Shiny dashboard for running small scale simulations. The R code for the app is located in the `app` directory. After downloading this repo to your local machine, you can run the app with the following command:

```
library(shiny)
runApp("app")
```

From here, you can change a few parameters and run your simulations. They will be stored on a folder that the user can specify in the dashboard.

# How to use: Scripts

We have provided scripts for running more large scale simulations. Within the `R` directory, the main files of interest are:

- `pipeline.R`: Script where parameters can be specified to initialize simulation
- `simulateData.R`: contains a function that simulates a Platform-of-1 trial. This version does not act on decision rules and is paired with `pipeline.R`
- `simulateDataDecision.R`: contains a function that simulates a Platform-of-1 trial. This version does act on decision rules.

Other files in the `R` folder hold code and functionality that help with simulating trials.


To use `pipeline.R`, a user will need to do the following: 

1. Specify where to place simulations
2. Define the parameters for the Platform-of-1 trial (observation frequency, standardized effect size, AR1 parameter)
3. Run the simulations

The code is designed to output the results of each individual simulation into a pre-specified directory. From there, they can be compiled into a single dataframe and then analyzed.

We recommend using the scripts if you plan to run > 1000 simulations. Our simulation pipeline takes advantage of parallelization and speeds up the process.
