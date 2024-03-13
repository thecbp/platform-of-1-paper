################################################################################

# DEFINE PACKAGES NEEDED TO RUN PIPELINE

################################################################################

# Load required packages, and install first if necessary
if(!require(pacman)) {install.packages("pacman")}
pkgs = c(
  "tidyverse",
  "brms",
  "furrr", 
  "glue", 
  "brms"
)

pacman::p_load(char = pkgs)
source("helpers.R")          # loads Thompson Sampling, rebalancing
source("simulatePlatformOf1.R") # loads simulatePlatformOf1() function
source("simulatePlatformOf1Decision.R") # loads simulatePlatformOf1Decision() function

################################################################################

# DEFINE FILE I/O

################################################################################

# [!!!] Specify where simulation results should be stored
out_path = "simulations/specify-folder-here"

# Create directory if it doesn't already exist
if (!dir.exists(out_path))
  dir.create(out_path)

################################################################################

# DEFINE DESIRED PARAMETERS FOR PLATFORM-OF-1 DESIGN

################################################################################

# [!!!] Define the number of treatments in the simulation
# should be a single number
choice_num_trts = 3

# [!!!] Define the number of cycles in the initialization phase
# should be single number
choice_num_init_cycles = 1

# [!!!] Define observation frequency made during initialization phase 
# can be single number or vector
choice_num_obvs_init = c(1)

# [!!!] Define observation frequency made during adaptive phase
choice_num_obvs_adaptive = c(1)

# [!!!] Define the maximum length of the trial
# should be a single number
choice_max_duration = 50

# [!!!] Define the effect size of the effective treatment
# can be single number or vector
choice_effect_size = c(0, -2, -5, -10)

# [!!!] Define the autocorrelation of the AR(1)
# can be single number or vector
choice_ar1_phi = c(0, 25, 50, 75)

# [!!!] Define the number of repetitions per simulation configuration
no_reps_per_scenario = 2000

################################################################################

# ADDITIONAL PARAMETERS FOR SIMULATION
# Thompson Sampling, MCMC, interim decision settings

################################################################################

MODEL_ARGS = list()

# [!!!] TS: Define the objective of the optimization function
MODEL_ARGS[["OBJECTIVE"]] = "Minimize"

# [!!!] TS: Define intercept of Thompson Sampling reward model
MODEL_ARGS[["INTERCEPT"]] = 130

# [!!!] TS: Define the standard deviation of Thompson Sampling of the reward model
MODEL_ARGS[["SIGMA"]] = 10

# [!!!] TS: Define priors for TS reward model (brms)
MODEL_ARGS[["PRIORS"]] = c(
  brms::set_prior("normal(0,100)", class = "Intercept"),
  brms::set_prior("normal(0,100)", class = "b")
)

# Pre-build Stan model to speed up simulations
X = diag(choice_num_trts)
X[, 1] = 1
X = X[rep(seq_len(nrow(X)), each = choice_num_obvs_init),]
X = X[rep(seq_len(nrow(X)), times = choice_num_init_cycles),]
Y = rnorm(nrow(X))
DATA = cbind(X, Y = Y) %>% tibble::as_tibble()
colnames(DATA) = c(paste0("X", 1:choice_num_trts), "Y", "period", "trt")

# Store the compiled model for use in simulations
MODEL_ARGS[["COMPILED_FIT"]] = brms::brm(formula = paste0("Y ~", paste0("X", 2:choice_num_trts, collapse = "+")),
                                         data = DATA,
                                         family = gaussian(link = "identity"),
                                         chains = MODEL_ARGS[["NUM_CHAINS"]],
                                         warmup = MODEL_ARGS[["NUM_ITER_WARMUP"]],
                                         iter = MODEL_ARGS[["NUM_ITER_POST"]],
                                         thin = MODEL_ARGS[["THINNING_FACTOR"]],
                                         prior = MODEL_ARGS[["PRIORS"]],
                                         control = list(
                                           max_treedepth = MODEL_ARGS[["MAX_TREE_DEPTH"]],
                                           adapt_delta = MODEL_ARGS[["ADAPT_DELTA"]]
                                         ),
                                         refresh = 0)

# [!!!] MCMC: Number of MCMC chains
MODEL_ARGS[["NUM_CHAINS"]] = 4

# [!!!] MCMC: Number of warm-up iterations
MODEL_ARGS[["NUM_ITER_WARMUP"]] = 5000

# [!!!] MCMC: Number of post warm-up iterations
MODEL_ARGS[["NUM_ITER_POST"]] = 10000

# [!!!] MCMC: Thinning factor
MODEL_ARGS[["THINNING_FACTOR"]] = 2

# [!!!] MCMC: Adapt delta
MODEL_ARGS[["ADAPT_DELTA"]] = 0.95

# [!!!] MCMC: Maximum tree depth
MODEL_ARGS[["MAX_TREE_DEPTH"]] = 13

# [!!!] Interim Decision: clinically relevant value for treatment (used with simulatePlatformOf1Decision())
MODEL_ARGS[["DELTA"]] = 130

# [!!!] Interim Decision: Add a threshold for declaring graduation (used with simulatePlatformOf1Decision())
MODEL_ARGS[["EFFICACY_THRESHOLD"]] = 0.95

# [!!!] Interim Decision: Add a threshold for declaring futility (used with simulatePlatformOf1Decision())
MODEL_ARGS[["FUTILITY_THRESHOLD"]] = 0.20

################################################################################

# RUN THE SIMULATIONS

################################################################################

# Check if simulation files have already been generated
already_simulated = list.files(out_path) %>% str_replace(".rds", "")

# Form a grid of all possible scenarios, giving each an unique id
scenario_grid =
  tidyr::expand_grid(
    REP_ID = seq_len(no_reps_per_scenario),
    NUM_TRTS = choice_num_trts,
    NUM_INIT_CYCLES = choice_num_init_cycles,
    NUM_OBVS_INIT = choice_num_obvs_init,
    NUM_OBVS_ADAPTIVE = choice_num_obvs_adaptive,
    MAX_DURATION = choice_max_duration,
    EFFECT_SIZE = choice_effect_size, 
    AR1_PHI = choice_ar1_phi
  ) %>% 
  mutate(
    EFFECT_SIZE_CODE = case_when(
      EFFECT_SIZE == 0 ~ "null",
      EFFECT_SIZE == -2 ~ "weak",
      EFFECT_SIZE == -5 ~ "mod",
      EFFECT_SIZE == -10 ~ "large",
    ),
    FILE = glue::glue("{EFFECT_SIZE_CODE}_len{NUM_OBVS_INIT}_AR{AR1_PHI}_{REP_ID}.out")
  ) %>% 
  # Filter out simulations that have already been performed
  filter(!(FILE %in% already_simulated)) 

# Set up parallelization
plan(multisession, workers = availableCores() - 1)

# Run simulations, results will be stored in path specified in out.path
scenario_grid %>% 
  mutate(
    sim = future_pmap(list(REP_ID, NUM_TRTS, NUM_INIT_CYCLES, NUM_OBVS_INIT, NUM_OBVS_ADAPTIVE,
                           MAX_DURATION, EFFECT_SIZE, AR1_PHI, FILE), 
                      .options = furrr_options(
                        scheduling = 1,
                        seed = TRUE
                      ),
                      .progress = TRUE,
                      function(REP_ID_, NUM_TRTS_, NUM_INIT_CYCLES_, NUM_OBVS_INIT_, NUM_OBVS_ADAPTIVE_,
                               MAX_DURATION_, EFFECT_SIZE_, AR1_PHI_, FILE_) {
                        
                        PARAM_ARGS = list(
                          REP_ID = REP_ID_,
                          NUM_TRTS = NUM_TRTS_,
                          NUM_INIT_CYCLES = NUM_INIT_CYCLES_,
                          NUM_OBVS_INIT = NUM_OBVS_INIT_,
                          NUM_OBVS_ADAPTIVE = NUM_OBVS_ADAPTIVE_,
                          MAX_DURATION = MAX_DURATION_,
                          EFFECT_SIZE = EFFECT_SIZE_, 
                          AR1_PHI = AR1_PHI_,
                          FILE = FILE_
                        )
                        
                        # For running simulations without interim decisions
                        simulatePlatformOf1(PARAM_ARGS, MODEL_ARGS)
                        
                        # For running simulations with interim decisions
                        # simulatePlatformOf1Decision(PARAM_ARGS, MODEL_ARGS)
                        
                      }
    ))
  