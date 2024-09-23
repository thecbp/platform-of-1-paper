################################################################################

# DEFINE PACKAGES NEEDED TO RUN PIPELINE

################################################################################

# Load required packages, and install first if necessary
if(!require(pacman)) {install.packages("pacman")}
pkgs = c(
  "tidyverse",
  "brms",
  "glue"
)

pacman::p_load(char = pkgs)
source("code/helpers.R")               # loads helper functions 
source("code/configs.R")                # loads simulation configurations
source("code/simulateDataV2Decision.R")        # loads function for generating data
source("code/simulateDataStandardDecision.R")        # loads function for generating data
source("code/ThompsonSamplingV2.R")

################################################################################

# DEFINE FILE I/O

################################################################################

# CHECK OUTPATH, scenariogrid to make sure correct DESIGN is used

# [!!!] Specify where simulation results should be stored
out_path = "simulations/P1-FIXED-decision-2000"

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
AR_PARAM = c(0.0, 0.3)

# [!!!] Define the number of repetitions per simulation configuration
NUM_REPS = 2000

# [!!!] Define which algorithm is being used in this set of simulations
ALG_ = "TS"

# [!!!] Define the number of agents used for consensus Thompson Sampling
# NUM_AGENTS = 2

################################################################################

# SETTING OTHER SIMULATION CONFIGURATION PARAMETERS

################################################################################

# Load in standard configuration file from separate file
CONFIG = loadP1Config(NUM_TRTS = 3)

# Define how many cores that replications will be used
# TSCC allows 100 maximum parallel jobs in an array job
NUM_CORES = 100

################################################################################

# CREATE GRID OF SIMULATION REPETITIONS

################################################################################

# Extract slurm job ID for parallel analysis
job_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(job_id)) { job_id = 1 }
set.seed(job_id)

# Check if simulation files have already been generated
done = glue::glue("{out_path}/{list.files(out_path)}")

# Form a grid of all possible replications needed for simulation study
scenario_grid =
  tidyr::expand_grid(
    REP_ID = seq_len(NUM_REPS),
    DESIGN = c("Fixed"),
    OBVS_FREQ = OBVS_FREQ,
    EFFECT_SIZE = EFFECT_SIZE, 
    AR_PARAM = AR_PARAM
  ) |> 
  left_join(read_csv("./processed/manuscript-fwer-cutoffs.csv"), by = c("DESIGN", "OBVS_FREQ", "AR_PARAM")) |> 
  left_join(read_csv("./processed/manuscript-t2e-cutoffs.csv"), by = c("DESIGN", "OBVS_FREQ", "AR_PARAM", "EFFECT_SIZE")) |> 
  rename(ETAE = min_etaf, ETAF = max_etae)

# Further processing on this scenario grid with specific core information
scenario_grid |> 
  dplyr::mutate(
    
    # Distribute replications among 1000 cores for array job
    JOB_ID = rep(1:NUM_CORES, length.out = nrow(scenario_grid)),
    
    ALG = ALG_,
    
    # Translate the effect size into a string for file naming
    EFFECT_SIZE_CODE = case_when(
      EFFECT_SIZE == 0.0 ~ "N",
      EFFECT_SIZE == -0.3 ~ "S",
      EFFECT_SIZE == -0.5 ~ "M",
      EFFECT_SIZE == -0.7 ~ "L",
    ),
    
    # Making AR_PARAM a string for consistent file naming
    AR_PARAM_STR = if_else(AR_PARAM == 0.0, "0.0", as.character(AR_PARAM)),
    
    # Assign a file name to each repetition
    FILE = glue::glue("{out_path}/{ALG}-ES{EFFECT_SIZE_CODE}-OF{OBVS_FREQ}-PHI{AR_PARAM_STR}-{REP_ID}.rds")
    
  ) |> 
  
  # Filter out simulations that have already been performed and 
  # only include replicates to those designated to this particular array ID
  filter(!(FILE %in% done), 
         JOB_ID == job_id) |> 
  dplyr::mutate(
    
    # Perform the relevant replications for this particular core
    sim = pwalk(list(REP_ID, OBVS_FREQ, EFFECT_SIZE, AR_PARAM, FILE, ETAE, ETAF), 
                function(ri, of, es, arp, f, etae, etaf) {
                  
                  # For running simulations without interim decisions
                  s = simulateDataStandardDecision(REP_ID = ri,
                                                   OBVS_FREQ = of,
                                                   EFFECT_SIZE = es,
                                                   AR_PARAM = arp,
                                                   CONFIG = CONFIG,
                                                   ALG = ALG_,
                                                   ETA_F = etaf,
                                                   ETA_E = etae
                  )
                  
                  saveRDS(s, file = f)
                  
                }
    )
  )
