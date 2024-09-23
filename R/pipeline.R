################################################################################

# DEFINE PACKAGES NEEDED TO RUN PIPELINE

################################################################################

# Load required packages, and install first if necessary
if(!require(pacman)) {install.packages("pacman")}
pkgs = c(
  "tidyverse",
  "brms",
  "glue",
  "furrr"
)

pacman::p_load(char = pkgs)
source("helpers.R")               # loads helper functions 
source("configs.R")                # loads simulation configurations
source("simulateData.R")        # loads function for generating data
source("code/ThompsonSamplingV2.R")

################################################################################

# DEFINE FILE I/O

################################################################################

# [!!!] Specify where simulation results should be stored
dir.create("simulations")
out_path = "simulations/[[CHANGE-THIS-STRING-TO-DESIRED-NAME]]"

# Create directory if it doesn't already exist
if (!dir.exists(out_path))
  dir.create(out_path)

################################################################################

# DEFINE DESIRED SIMULATION PARAMETERS

################################################################################

# [!!!] Define observation frequency for each treatment period
OBVS_FREQ = c(5, 10, 15)

# [!!!] Define the effect size (mu/sigma) of the effective treatment
EFFECT_SIZE = c(0.0, -0.3, -0.7)

# [!!!] Define AR(p) parameter (x10 for consistent naming in files)
AR_PARAM = c(0.0, 0.3, 0.7)

# [!!!] Define the number of repetitions per simulation configuration
NUM_REPS = 2000

################################################################################

# SETTING OTHER SIMULATION CONFIGURATION PARAMETERS

################################################################################

# Load in standard configuration file from separate file
CONFIG = loadP1Config(NUM_TRTS = 3)

################################################################################

# CREATE GRID OF SIMULATION REPETITIONS

################################################################################

# Check if simulation files have already been generated
done = glue::glue("{out_path}/{list.files(out_path)}")

# Form a grid of all possible scenarios, giving each an unique id
scenario_grid =
  tidyr::expand_grid(
    REP_ID = seq_len(NUM_REPS),
    OBVS_FREQ = OBVS_FREQ,
    EFFECT_SIZE = EFFECT_SIZE, 
    AR_PARAM = AR_PARAM
  ) |> 
  dplyr::mutate(
    
    ALG = "TS",
    
    # Making AR_PARAM a string for consistent file naming
    AR_PARAM_STR = if_else(AR_PARAM == 0.0, "0.0", as.character(AR_PARAM)),
    
    # Assign a file name to each repetition
    FILE = glue::glue("{out_path}/{ALG}-ES{EFFECT_SIZE}-OF{OBVS_FREQ}-PHI{AR_PARAM_STR}-{REP_ID}.rds")
    
  ) |> 
  filter(!(FILE %in% done)) # Filter out simulations that have already been performed

# Set up parallelization for purrr
plan(multisession, workers = availableCores() - 1)

# Run simulations
# Results will be stored in path specified in out.path variable
scenario_grid |> 
  mutate(
    sim = future_pwalk(list(REP_ID, OBVS_FREQ, EFFECT_SIZE, AR_PARAM, FILE), 
                       .options = furrr_options(scheduling = 1, seed = TRUE),
                       .progress = TRUE,
                       function(ri, of, es, arp, f) {
                         
                         # For running simulations without interim decisions
                         s = simulateData(REP_ID = ri,
                                          OBVS_FREQ = of,
                                          EFFECT_SIZE = es,
                                          AR_PARAM = arp,
                                          CONFIG = CONFIG,
                                          ALG = "TS"
                         )
                         
                         saveRDS(s, file = f)
                         
                       }
    ))
