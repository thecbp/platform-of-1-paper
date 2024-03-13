simulatePlatformOf1Decision = function(PARAM_ARGS, MODEL_ARGS) {
  
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
  
  # Initialize table to hold data on distance of posterior mean to set point
  PP_SPDIFF_TABLE = tibble()
  
  CLIN_INTERVAL_TABLE = tibble()
  
  FUTILITY_TABLE = tibble(
    treatment = paste0("X", 1:PARAM_ARGS[["NUM_TRTS"]]),
    futile = FALSE, 
    futility_declared = NA
  )
  
  ##############################################################################
  # Adaptive phase
  ##############################################################################
  
  adaptive_start = (PARAM_ARGS[["NUM_TRTS"]] * PARAM_ARGS[["NUM_INIT_CYCLES"]]) + 1
  adaptive_periods = adaptive_start:PARAM_ARGS[["MAX_DURATION"]]
  starting_treatments = current_treatments = 1:PARAM_ARGS[["NUM_TRTS"]]
  
  for (period in adaptive_periods) {
    
    # Select the next treatment based on reallocated probabilities
    next_trt = sample(current_treatments, size = 1, prob = current_probs)
    
    # Construct design matrix and "observe" the outcome under this treatment
    X = matrix(0, nrow = PARAM_ARGS[["NUM_OBVS_ADAPTIVE"]], ncol = PARAM_ARGS[["NUM_TRTS"]])
    X[, 1] = 1        # intercept
    X[, next_trt] = 1 # setting next treatment
    
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
    
    # doing the check for distance of posterior mean to set point
    SPDIFF = posterior::as_draws_df(POSTERIOR) %>% 
      select(contains("b_Intercept"), contains("b_X")) %>% 
      rename(b_X1 = b_Intercept) %>% 
      mutate(across(-b_X1, ~ . + b_X1)) %>% 
      pivot_longer(everything(), names_to = "treatment", values_to = "val") %>% 
      group_by(treatment) %>% 
      summarize(mu = mean(val)) %>% 
      mutate(
        spdiff = abs(mu - MODEL_ARGS[["SET_POINT"]])
      ) 
    
    PP_SPDIFF_ROW = SPDIFF %>% 
      pivot_wider(values_from = c("mu", "spdiff"), names_from = treatment) %>% 
      mutate(
        trt_min_dist = which.min(SPDIFF$spdiff),
        period = period
      )
    
    # calculating how much of the posterior distribution falls within a 
    # also record whether or not a decision was made
    CLIN_INTERVAL = posterior::as_draws_df(POSTERIOR) %>% 
      select(contains("b_Intercept"), contains("b_X")) %>% 
      rename(b_X1 = b_Intercept) %>% 
      mutate(across(-b_X1, ~ . + b_X1)) %>% 
      pivot_longer(everything(), names_to = "treatment", values_to = "val") %>% 
      group_by(treatment) %>% 
      summarize(
        in_interval = mean(val <= MODEL_ARGS[["DELTA"]])
      ) %>% 
      mutate(
        futility_met = in_interval <= MODEL_ARGS[["FUTILITY_THRESHOLD"]],
        graduation_met = in_interval >= MODEL_ARGS[["EFFICACY_THRESHOLD"]]
      )
    
    CLIN_INTERVAL_ROW = CLIN_INTERVAL %>% 
      pivot_wider(values_from = c("in_interval", "futility_met", "graduation_met"), 
                  names_from = treatment) %>% 
      mutate( period = period )
    
    # Check if any treatment has gone below the futility threshold
    if (any(CLIN_INTERVAL$futility_met[current_treatments]) & period >= MODEL_ARGS[["DECISION_START"]]) { 
      
      is_futile = which(CLIN_INTERVAL$futility_met)[1]
      current_treatments = setdiff(current_treatments, is_futile)
      
      FUTILITY_TABLE$futile[is_futile] = TRUE
      FUTILITY_TABLE$futility_declared[is_futile] = period
      
    }
    
    # Calculate allocation probabilities via on Thompson Sampling
    current_probs = TS(POSTERIOR, 
                       c = (period) / (2 * PARAM_ARGS[["MAX_DURATION"]]), 
                       objective = MODEL_ARGS[["OBJECTIVE"]]) 
    current_probs = rebalance_probabilities(current_probs,
                                            current = current_treatments, 
                                            start = starting_treatments)
    
    # Append new set of allocation probabilities to current set
    # NOTE TO SELF: this doesnt work for this simulation set!
    ALLOC_PROB_TABLE = rbind(ALLOC_PROB_TABLE, c(current_probs, period))
    
    # Do the same for the distances from the set point
    PP_SPDIFF_TABLE = bind_rows(PP_SPDIFF_TABLE, PP_SPDIFF_ROW)
    
    CLIN_INTERVAL_TABLE = bind_rows(CLIN_INTERVAL_TABLE, CLIN_INTERVAL_ROW)
    
    # Check if any treatment has surpassed the graduation threshold
    # Record period and treatment that has graduated
    if (any(CLIN_INTERVAL$graduation_met) & 
        length(which(CLIN_INTERVAL$graduation_met)) == 1 & 
        period >= MODEL_ARGS[["DECISION_START"]]) {
      
      DECISION = tibble(
        finish = period, 
        chosen = which(CLIN_INTERVAL$graduation_met)
      )
      
      break
      
    }
    
    # If there's only one treatment left, just graduate with it
    if (length(current_treatments) == 1 & 
        period >= MODEL_ARGS[["DECISION_START"]]) {
      
      DECISION = tibble(
        finish = period, 
        chosen = current_treatments
      )
      
      break
      
    }
    
  }
  
  ##############################################################################
  # Post simulation processing
  ##############################################################################
  
  # if no decision reached by maximum duration, choose the treatment whose mean is closest to the set point
  if (!exists("DECISION")) {
    DECISION = tibble(
      finish = PARAM_ARGS[["MAX_DURATION"]],
      chosen = which.min(SPDIFF$spdiff)
    )
  }
  
  ALLOC_PROB_TABLE = tibble::as_tibble(ALLOC_PROB_TABLE)
  names(ALLOC_PROB_TABLE) = c(paste0("prob_X", 1:PARAM_ARGS[["NUM_TRTS"]]), "period")
  
  # Convert allocation probability table to wider format for simulation export
  ALLOC_PROB_TABLE_WIDE = ALLOC_PROB_TABLE %>% 
    pivot_longer(contains("prob"), names_to = "trt", values_to = "prob") %>% 
    transmute(col = paste0(trt, "_period", period), prob) %>% 
    pivot_wider(names_from = col, values_from = prob)
  
  # Convert the setpoint differences into wider format
  SPDIFF_TABLE_WIDE = PP_SPDIFF_TABLE %>% 
    select(period, trt_min_dist, contains("spdiff")) %>% 
    pivot_wider(names_from = period, 
                values_from = c("trt_min_dist", "spdiff_b_X1", "spdiff_b_X2", "spdiff_b_X3"))
  
  CLIN_INTERVAL_TABLE_WIDE = CLIN_INTERVAL_TABLE %>% 
    pivot_wider(names_from = period, 
                values_from = everything())
  
  FUTILITY_WIDE = FUTILITY_TABLE %>% 
    pivot_wider(names_from = treatment, values_from = c("futile", "futility_declared"))
  
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
  OUT = bind_cols(OUT, SPDIFF_TABLE_WIDE)
  OUT = bind_cols(OUT, CLIN_INTERVAL_TABLE_WIDE)
  OUT = bind_cols(OUT, ASSIGNMENTS_WIDE)
  OUT = bind_cols(OUT, FUTILITY_WIDE)
  OUT = bind_cols(OUT, DECISION)
  
  saveRDS(OUT, file = sprintf("%s/%s.rds", out_path, PARAM_ARGS[["FILE"]]))
  
  return(OUT)
  
}