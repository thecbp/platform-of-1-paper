simulateDataV2Decision = function(REP_ID, 
                                  OBVS_FREQ, 
                                  EFFECT_SIZE, 
                                  AR_PARAM, 
                                  CONFIG, 
                                  ALG, 
                                  ETA_F, 
                                  ETA_E,
                                  NUM_AGENTS = NULL) {
  
  # Start timer for simulation
  tic = lubridate::now()
  
  # Initialize output row with given simulation parameters
  OUT = tibble(
    REP_ID = REP_ID,
    OBVS_FREQ = OBVS_FREQ, 
    EFFECT_SIZE = EFFECT_SIZE, 
    AR_PARAM = AR_PARAM,
    NUM_AGENTS = NUM_AGENTS,
    CONFIG = CONFIG[["NAME"]],
    ALG = ALG, 
    DESIGN = "Adaptive"
  )
  
  ##############################################################################
  # Initialization phase
  ##############################################################################
  
  # Design matrix of single person getting all of the treatments for n_obvs each
  X = diag(CONFIG[["NUM_TRTS"]])
  X[, 1] = 1
  X = X[rep(seq_len(nrow(X)), each = OBVS_FREQ),]
  
  # Generate outcome based on given AR model
  # For initialization, all 3 treatments are observed, so we must multiply to account for this
  Y = generateOutcomes(X, EFFECT_SIZE, CONFIG[["NUM_TRTS"]] * OBVS_FREQ, AR_PARAM, CONFIG)
  
  DATA = cbind(X, Y = Y) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      period = rep(1:CONFIG[["NUM_TRTS"]], each = OBVS_FREQ),
      trt = rep(1:CONFIG[["NUM_TRTS"]], each = OBVS_FREQ),
      obs = row_number()
    )
  
  colnames(DATA) = c(paste0("X", 1:CONFIG[["NUM_TRTS"]]), "Y", "period", "trt", "obs")
  
  # Update the brms model to use the new data
  FIT = update(CONFIG[["COMPILED_FIT"]], newdata = DATA)
  
  # Calculate allocation probabilities via Thompson Sampling (or other algorithm)
  stab_parameter = CONFIG[["NUM_TRTS"]] / (2 * CONFIG[["MAX_DURATION"]])
  current_probs = adaptiveRandomizationV2(FIT, DATA, ALG = "TS", STAB = stab_parameter) 
  
  # Initialize tibbles to hold metrics for simulation
  ALLOC_PROB_TABLE = c(current_probs, CONFIG[["NUM_TRTS"]])
  PP_EFFECTIVE_TABLE = tibble()
  ESTIMATES_TABLE = tibble() # bias and credible interval of parameter
  
  # Initialize pool of treatments
  current_treatments = 1:CONFIG[["NUM_TRTS"]]
  
  # Initialize efficacy and futility columns
  OUT[["GRADUATED_X2"]] = FALSE
  OUT[["GRADUATED_X3"]] = FALSE
  OUT[["TRIAL_FINISH"]] = FALSE
  OUT[["FUTILE_X2"]] = FALSE
  OUT[["FUTILE_X3"]] = FALSE
  OUT[["PERIOD_REMOVED_X2"]] = NA
  OUT[["PERIOD_REMOVED_X3"]] = NA
 
  ##############################################################################
  # Adaptive phase
  ##############################################################################
  
  adaptive_periods = (CONFIG[["NUM_TRTS"]] + 1):CONFIG[["MAX_DURATION"]]
  
  for (period in adaptive_periods) {

    ############################################################################
    # Generate data for an adaptive treatment assignment
    ############################################################################
    
    # Select the next treatment
    next_trt = sample(current_treatments, size = 1, prob = current_probs[current_treatments])
    
    # Construct design matrix and "observe" the outcome under this treatment
    X = matrix(0, nrow = OBVS_FREQ, ncol = CONFIG[["NUM_TRTS"]])
    X[, 1] = 1        # intercept
    X[, next_trt] = 1 # setting next treatment
    Y = generateOutcomes(X, EFFECT_SIZE, OBVS_FREQ, AR_PARAM, CONFIG)
    
    period_cols = matrix(rep(period, OBVS_FREQ), nrow = OBVS_FREQ)
    assignment_cols = matrix(rep(next_trt, OBVS_FREQ), nrow = OBVS_FREQ)
    # agent_cols = rep(1:NUM_AGENTS, length.out = OBVS_FREQ)
    obs_cols = (nrow(DATA) + 1):(nrow(DATA) + OBVS_FREQ)

    new_data = cbind(X, Y, period_cols, assignment_cols, obs_cols)
    colnames(new_data) = c(paste0("X", 1:CONFIG[["NUM_TRTS"]]), "Y", "period", "trt", "obs")
    
    # Append new data to current set
    DATA = rbind(DATA, new_data)
    
    ############################################################################
    # Calculate metrics of interest from the posterior distribution
    ############################################################################
    
    # Calculate the posterior probability that a parameter > DELTA (since "Maximize")
    PP_EFFECTIVE = getPosteriorSamples(FIT, DATA, ALG) |> 
      select(contains("b_Intercept"), contains("b_X")) |> 
      rename(b_X1 = b_Intercept) |> 
      mutate(across(-b_X1, ~ . + b_X1)) |> 
      pivot_longer(everything(), names_to = "treatment", values_to = "sample") |> 
      group_by(treatment) |> 
      summarize(
        # prob_eff = mean(sample > CONFIG[["DELTA"]]), # For maximizing
        prob_eff = mean(sample < CONFIG[["DELTA"]])    # For minimizing
      ) |> 
      mutate(
        futility_met = prob_eff <= ETA_F,
        graduation_met = prob_eff >= ETA_E
      )
    
    # Transform df into table for storage
    PP_EFFECTIVE_ROW = PP_EFFECTIVE |> 
      select(-futility_met, -graduation_met) |> 
      pivot_wider(values_from = "prob_eff", 
                  names_from = treatment) |> 
      mutate( period = period )
    
    # Calculate the posterior mean and the 95% credible interval for the regression parameters
    ESTIMATES_ROW = getPosteriorSamples(FIT, DATA, ALG) |> 
      select(contains("b_Intercept"), contains("b_X")) |> 
      rename(b_X1 = b_Intercept) |> 
      pivot_longer(everything(), names_to = "treatment", values_to = "sample") |> 
      group_by(treatment) |> 
      summarize(
        post_mean = mean(sample),
        q025 = quantile(sample, 0.025),
        q975 = quantile(sample, 0.975)
      ) |> 
      pivot_wider(values_from = c("post_mean", "q025", "q975"), 
                  names_from = treatment) |> 
      mutate( period = period )
    
    # Append results of simulation to the ongoing set of results
    ALLOC_PROB_TABLE = rbind(ALLOC_PROB_TABLE, c(current_probs, period))
    PP_EFFECTIVE_TABLE = bind_rows(PP_EFFECTIVE_TABLE, PP_EFFECTIVE_ROW)
    ESTIMATES_TABLE = bind_rows(ESTIMATES_TABLE, ESTIMATES_ROW)
    
    # Perform futility and graduation checks
    grad_checks = PP_EFFECTIVE |> pull(graduation_met)
    if (grad_checks[next_trt] & next_trt != 1) {
      
      OUT[[glue::glue("GRADUATED_X{next_trt}")]] = TRUE
      OUT[["TRIAL_FINISH"]] = period
      break # end adaptive phase due to efficacy
      
    }
    
    futile_checks = PP_EFFECTIVE |> pull(futility_met)
    if (futile_checks[next_trt] & next_trt != 1) {
      
      OUT[[glue::glue("FUTILE_X{next_trt}")]] = TRUE
      OUT[[glue::glue("PERIOD_REMOVED_X{next_trt}")]] = period
      
      # Update current set of arms
      current_treatments = setdiff(current_treatments, next_trt)
      
    }
    
    if (length(current_treatments) == 1) { 
      OUT[["TRIAL_FINISH"]] = period
      break  
    }
    
    # Rebalance allocation probabilities according to updated posterior
    stab_parameter = period / (2 * CONFIG[["MAX_DURATION"]])
    current_probs = adaptiveRandomizationV2(FIT = FIT, DATA = DATA, ALG = ALG, STAB = stab_parameter)
    
  }
  
  ##############################################################################
  # Post simulation processing
  # Take all output tables and expand into one giant row
  ##############################################################################
  
  ALLOC_PROB_TABLE = tibble::as_tibble(ALLOC_PROB_TABLE)
  names(ALLOC_PROB_TABLE) = c(paste0("prob_X", 1:CONFIG[["NUM_TRTS"]]), "period")
  
  # Convert allocation probability table to wider format for simulation export
  ALLOC_PROB_TABLE_WIDE = ALLOC_PROB_TABLE |> 
    pivot_longer(contains("prob"), names_to = "trt", values_to = "prob") |> 
    transmute(col = paste0(trt, "_period", period), prob) |> 
    pivot_wider(names_from = col, values_from = prob)
  
  PP_EFFECTIVE_TABLE_WIDE = PP_EFFECTIVE_TABLE |> 
    pivot_wider(names_from = period, 
                values_from = contains("b_"))
  
  ESTIMATES_TABLE_WIDE = ESTIMATES_TABLE |> 
    pivot_wider(names_from = period, values_from = contains("b_"))
   
  # Convert the actual treatment assignment into wide
  ASSIGNMENTS_WIDE = DATA |> 
    transmute(
      col = glue::glue("assigned_{period}"),
      val = trt
    ) |> 
    distinct() |> 
    pivot_wider(names_from = col, values_from = val)
  
  # Combine everything into single output row
  OUT = bind_cols(OUT, ALLOC_PROB_TABLE_WIDE)
  OUT = bind_cols(OUT, PP_EFFECTIVE_TABLE_WIDE)
  OUT = bind_cols(OUT, ASSIGNMENTS_WIDE)
  OUT = bind_cols(OUT, ESTIMATES_TABLE_WIDE)
  
  # Finally, log the time it took to finish the 
  toc = lubridate::now()
  OUT[["DURATION"]] = toc - tic
  
  OUT

}