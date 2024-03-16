source("ui-code.R")
source("../helpers.R")
source("../simulatePlatformOf1Shiny.R")

library(tidyverse)

ui = fluidPage(
  withMathJax(),
  navbarPage("Platform-of-1 Planner",
             tabPanel("Simulation", simulationTab),
             tabPanel("How To Use", howtoTab)
  )
)
