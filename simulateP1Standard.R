simulatePlatformOf1Standard = function(PARAM_ARGS, MODEL_ARGS) {
  
  # Initialize output row with 
  OUT = tibble(
    REP_ID = PARAM_ARGS[["REP_ID"]],
    NUM_TRTS = PARAM_ARGS[["NUM_TRTS"]], 
    NUM_INIT_CYCLES = PARAM_ARGS[["NUM_INIT_CYCLES"]],
    NUM_OBVS_INIT = PARAM_ARGS[["NUM_OBVS_INIT"]],
    NUM_OBVS_ADAPTIVE = PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]],
    MAX_DURATION = PARAM_ARGS[["MAX_DURATION"]],
    EFFECT_SIZE = PARAM_ARGS[["EFFECT_SIZE"]],
    AR1_PHI = PARAM_ARGS[["AR1_PHI"]]
  )
  
  ##############################################################################
  # Initialization phase
  ##############################################################################
  
  # Design matrix of single person getting all of the treatments for n_obvs each
  X = diag(PARAM_ARGS[["NUM_TRTS"]])
  X[, 1] = 1
  X = X[rep(seq_len(nrow(X)), each = PARAM_ARGS[["NUM_OBVS_INIT"]]),]
  X = X[rep(seq_len(nrow(X)), times = PARAM_ARGS["NUM_INIT_CYCLES"]),]
  
  BETAS = rep(0, PARAM_ARGS[["NUM_TRTS"]])# Initialize vector of null effects
  BETAS[1] = MODEL_ARGS[["INTERCEPT"]]
  BETAS[2] = PARAM_ARGS[["EFFECT_SIZE"]]
  print(BETAS)
  
  # Generate outcome based on multiple regression, if no serial correlation
  # Otherwise, use an ARIMA model to generate it
  if (PARAM_ARGS[["AR1_PHI"]] == 0) {
    Y = ((X %*% BETAS) %>% c()) + stats::rnorm(nrow(X), 0, sd = MODEL_ARGS[["SIGMA"]])
  } else {
    Y_ar = arima.sim(model = list(ar = PARAM_ARGS[["AR1_PHI"]] / 100), 
                     sd = MODEL_ARGS[["SIGMA"]],
                     n = PARAM_ARGS[["NUM_TRTS"]] * 
                       PARAM_ARGS[["NUM_INIT_CYCLES"]] * 
                       PARAM_ARGS[["NUM_OBVS_INIT"]])
    Y = Y_ar + (X %*% BETAS)
  }
  
  # Initialize ongoing trial data
  TRIAL_DATA = cbind(X, Y = Y) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      period = rep(1:PARAM_ARGS[["NUM_TRTS"]], each = PARAM_ARGS[["NUM_INIT_CYCLES"]] * PARAM_ARGS[["NUM_OBVS_INIT"]]),
      trt = rep(1:PARAM_ARGS[["NUM_TRTS"]], each = PARAM_ARGS[["NUM_INIT_CYCLES"]] * PARAM_ARGS[["NUM_OBVS_INIT"]])
    )
  
  colnames(TRIAL_DATA) = c(paste0("X", 1:PARAM_ARGS[["NUM_TRTS"]]), "Y", "period", "trt")
  
  # Generate first posterior distribution from intiailization phase data
  POSTERIOR = brms::brm(formula = FORMULA,
                        data = TRIAL_DATA,
                        family = gaussian(link = "identity"),
                        fit = MODEL_ARGS[["COMPILED_FIT"]],
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
  
  # Calculate allocation probabilities via on Thompson Sampling
  current_probs = TS(POSTERIOR, 
                     c = (PARAM_ARGS[["NUM_TRTS"]] * PARAM_ARGS[["NUM_INIT_CYCLES"]]) / (2 * PARAM_ARGS[["MAX_DURATION"]]), 
                     objective = MODEL_ARGS[["OBJECTIVE"]]) 
  
  # Record the first set of allocation probabilities
  ALLOC_PROB_TABLE = c(current_probs, PARAM_ARGS[["NUM_TRTS"]] * PARAM_ARGS[["NUM_INIT_CYCLES"]])
  
  # Initialize table to hold futility and efficacy decisions
  DECISION_CUTOFF_TABLE = tibble()
  
  ##############################################################################
  # Adaptive phase
  ##############################################################################
  
  adaptive_start = (PARAM_ARGS[["NUM_TRTS"]] * PARAM_ARGS[["NUM_INIT_CYCLES"]]) + 1
  adaptive_periods = adaptive_start:PARAM_ARGS[["MAX_DURATION"]]
  
  for (period in adaptive_periods) {
    
    # Select the next treatment based on reallocated probabilities
    next_trt = sample(1:PARAM_ARGS[["NUM_TRTS"]], size = 1, prob = current_probs)
    
    # Construct design matrix and "observe" the outcome under this treatment
    X = matrix(0, nrow = PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]], ncol = PARAM_ARGS[["NUM_TRTS"]])
    X[1] = 1        # intercept
    X[next_trt] = 1 # setting next treatment
    
    # Generate outcome based on if we need serial correlation or not
    if (PARAM_ARGS[["AR1_PHI"]]  == 0) {
      Y = (X %*% BETAS) + stats::rnorm(PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]], 0, MODEL_ARGS[["SIGMA"]])
    } else {
      Y_ar = arima.sim(model = list(ar = PARAM_ARGS[["AR1_PHI"]] / 100), sd = MODEL_ARGS[["SIGMA"]], n = PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]])
      Y = Y_ar + (X %*% BETAS)
    }
    
    period_col = matrix(rep(period, PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]]), nrow = PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]])
    assignment_col = matrix(rep(next_trt, PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]]), nrow = PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]])
    
    new_data = cbind(X, Y, period_col, assignment_col)
    
    # Rename this new data to make sure that columns are aligned when combining
    colnames(new_data) = c(paste0("X", 1:PARAM_ARGS[["NUM_TRTS"]]), "Y", "period", "trt")
    
    # Append new data to current set
    TRIAL_DATA = rbind(TRIAL_DATA, new_data)
    
    # Recalculate the posteriors using the updated dataset
    POSTERIOR = update(POSTERIOR, newdata = TRIAL_DATA)
    
    # Based on posterior, calculate the minimum cutoff needed to make either a
    # graduation or futility decision 
    decisions = posterior::as_draws_df(POSTERIOR) %>% 
      select(contains("b_X")) %>% 
      pivot_longer(contains("b_X"), 
                   names_to = "treatment",
                   values_to = "val") %>% 
      nest(data = c("val")) %>% 
      mutate(
        # Smallest efficacy cutoff that results in a graduation decision
        smallest_gamma_eff = map(data, function(postsamp) {
          tibble(gamma = seq(0.001, 0.999, by = 0.001)) %>%
            mutate(effdec_satisfied = map_lgl(gamma, function(g) { mean(postsamp$val < MODEL_ARGS[["DELTA_EFF"]]) > g } )) %>% 
            filter(effdec_satisfied) %>% 
            pull(gamma) %>% 
            max
        }),
        # Largest futility cutoff that results in a futility decision
        largest_gamma_fut = map(data, function(postsamp) {
          tibble(gamma = seq(0.001, 0.999, by = 0.001)) %>%
            mutate(futdec_satisfied = map_lgl(gamma, function(g) { mean(postsamp$val < MODEL_ARGS[["DELTA_EFF"]]) < g } )) %>% 
            filter(futdec_satisfied) %>% 
            pull(gamma) %>% 
            min
        })
      ) %>% 
      select(-data) %>% 
      pivot_longer(smallest_gamma_eff:largest_gamma_fut,
                   names_to = "decision_type",
                   values_to = "gamma") %>% 
      transmute(
        col = glue::glue("{treatment}_{decision_type}"),
        gamma = gamma
      ) %>% 
      pivot_wider(names_from = col, values_from = gamma) %>% 
      mutate( period = period )
    
    # Calculate allocation probabilities via on Thompson Sampling
    current_probs = TS(POSTERIOR, 
                       c = (period) / (2 * PARAM_ARGS[["MAX_DURATION"]]), 
                       objective = MODEL_ARGS[["OBJECTIVE"]]) 
    
    # Append new set of allocation probabilities to current set
    ALLOC_PROB_TABLE = rbind(ALLOC_PROB_TABLE, c(current_probs, period))
    
    # Do the same for the decision cutoffs
    DECISION_CUTOFF_TABLE = bind_rows(DECISION_CUTOFF_TABLE, decisions)
    
  }
  
  ##############################################################################
  # Post simulation processing
  ##############################################################################
  
  ALLOC_PROB_TABLE = tibble::as_tibble(ALLOC_PROB_TABLE)
  names(ALLOC_PROB_TABLE) = c(paste0("prob_X", 1:PARAM_ARGS[["NUM_TRTS"]]), "period")
  
  # Convert allocation probability table to wider format for simulation export
  ALLOC_PROB_TABLE_WIDE = ALLOC_PROB_TABLE %>% 
    pivot_longer(contains("prob"), names_to = "trt", values_to = "prob") %>% 
    transmute(col = paste0(trt, "_period", period), prob) %>% 
    pivot_wider(names_from = col, values_from = prob)
  
  # Convert the decision table into wider format
  DECISON_CUTOFF_TABLE_WIDE = DECISION_CUTOFF_TABLE %>% 
    pivot_longer(-period, 
                 names_to = "cutoff",
                 values_to = "val") %>% 
    transmute(
      col = glue::glue("{cutoff}_{period}"),
      val = val
    ) %>% 
    pivot_wider(names_from = col, values_from = val)
  
  # Convert the actual treatment assignment into wide
  ASSIGNMENTS_WIDE = TRIAL_DATA %>% 
    transmute(
      col = glue::glue("assigned_{period}"),
      val = trt
    ) %>% 
    distinct() %>% 
    pivot_wider(names_from = col, values_from = val)
  
  
  # Combine everything into single output row
  OUT = bind_cols(OUT, ALLOC_PROB_TABLE_WIDE)
  OUT = bind_cols(OUT, DECISON_CUTOFF_TABLE_WIDE)
  OUT = bind_cols(OUT, ASSIGNMENTS_WIDE)
  
  saveRDS(OUT, file = sprintf("%s/%s.rds", out_path, PARAM_ARGS[["FILE"]]))
  
  return(OUT)
  
}