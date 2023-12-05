################################################################################

# DEFINE PACKAGES NEEDED TO RUN PIPELINE

################################################################################

# Load required packages, and install first if necessary
if(!require(pacman)) {install.packages("pacman")}
pkgs = c(
  "tidyverse",
  "brms",
  "furrr"
)

pacman::p_load(char = pkgs)
source("P1-helpers.R") # Loads Thompson Sampling function
source("simulateP1Delayed.R")

################################################################################

# DEFINE FILE I/O

################################################################################

# [!!!] CHOOSE PARAMETERS: OUTPUT FOLDER
out_path = "simulations/delayed"

# Create directory if it doesn't already exist
if (!dir.exists(out_path))
  dir.create(out_path)

################################################################################

# DEFINE DESIRED PARAMETERS FOR CURRENT SET OF SIMULATIONS

################################################################################

# [!!!] Define the number of treatments in the simulation
choice_num_trts = 3

# [!!!] Define the number of cycles in the initialization phase
choice_num_init_cycles = 1

# [!!!] Define the number of observations made during a treatment period (initialization)
choice_num_obvs_init = 1

# [!!!] Define the number of observations made during a treatment period (adaptive)
choice_num_obvs_adaptive = 1

# [!!!] Define the maximum length of the trial
choice_max_duration = 50

# [!!!] Define the effect size of the effective treatment
choice_effect_size = c(0, -2, -5, -10)

# [!!!] Define the autocorrelation of the AR(1)
choice_ar1_phi = c(0)

# [!!!] CHOOSE PARAMETERS: NUMBER OF REPLICATIONS
no_reps_per_scenario = 2000

################################################################################

# ADDITIONAL PREPROCESSING 

################################################################################

MODEL_ARGS = list()

# [!!!] Define the objective of the optimization function
MODEL_ARGS[["OBJECTIVE"]] = "Minimize"

# [!!!] Define intercept of reward model
MODEL_ARGS[["INTERCEPT"]] = 130

# [!!!] Define the variance of the reward model
MODEL_ARGS[["SIGMA"]] = 10

# [!!!] Efficacy decision cutoff
MODEL_ARGS[["DELTA_EFF"]] = -10

# [!!!] Futility decision cutoff
MODEL_ARGS[["DELTA_FUT"]] = -10

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

# Setting priors for the model
MODEL_ARGS[["PRIORS"]] = c(
  brms::set_prior("normal(0,100)", class = "Intercept"),
  brms::set_prior("normal(0,100)", class = "b")
)

# Pre-build Stan models to speed up simulation speed

# For first portion of trials, before optimal treatment is added
# Design matrix of single person getting all of the treatments for n_obvs each
X = diag(choice_num_trts)
X[, 1] = 1
X = X[rep(seq_len(nrow(X)), each = choice_num_obvs_init),]
X = X[rep(seq_len(nrow(X)), times = choice_num_init_cycles),]
Y = rnorm(nrow(X))
DATA = cbind(X, Y = Y) %>% tibble::as_tibble()
colnames(DATA) = c(paste0("X", 1:choice_num_trts), "Y", "period", "trt")

# Build the formula for the model (X1 is reference treatment)
FORMULA = paste0("Y ~", paste0("X", 2:choice_num_trts, collapse = "+"))

MODEL_ARGS[["COMPILED_FIT"]] = brms::brm(formula = FORMULA,
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

# Another model for running when the optimal treatment is added
X = diag(choice_num_trts + 1)
X[, 1] = 1
X = X[rep(seq_len(nrow(X)), each = choice_num_obvs_init),]
X = X[rep(seq_len(nrow(X)), times = choice_num_init_cycles),]
Y = rnorm(nrow(X))
DATA = cbind(X, Y = Y) %>% tibble::as_tibble()
colnames(DATA) = c(paste0("X", 1:(choice_num_trts+1)), "Y", "period", "trt")

# Build the formula for the model (X1 is reference treatment)
FORMULA = paste0("Y ~", paste0("X", 2:(choice_num_trts+1), collapse = "+"))

MODEL_ARGS[["COMPILED_FIT2"]] = brms::brm(formula = FORMULA,
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

################################################################################

# RUNNING THE SIMULATIONS

################################################################################

# Check the currently constructed files
already_simulated = list.files(out_path) %>% str_replace(".rds", "")

# Form a grid of all possible scenarios, giving each an unique id
# Filter out files are already made in case of crashes, need to stop, etc. 
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
    FILE = glue::glue("delayed_{EFFECT_SIZE_CODE}_{REP_ID}.out")
  ) %>% 
  filter(!(FILE %in% already_simulated)) 

plan(multisession, workers = availableCores())
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
                        
                        simulatePlatformOf1Delayed(PARAM_ARGS, MODEL_ARGS)
                        
                      }
    ))
  