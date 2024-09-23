if(!require(pacman)) {install.packages("pacman")}
pkgs = c(
  "tidyverse",
  "brms",
  "glue",
  "furrr", 
  "foreach"
)

pacman::p_load(char = pkgs)
source("ui-code.R")
source("../R/helpers.R")               # loads helper functions 
source("../R/configs.R")                # loads simulation configurations
source("../R/simulateDataV2.R")        # loads function for generating data
source("../R/simulateDataStandard.R")        # loads function for generating data
source("../R/ThompsonSamplingV2.R")

ui = fluidPage(
  withMathJax(),
  navbarPage("Platform-of-1 Planner",
             tabPanel("Simulation", simulationTab)
  )
)
