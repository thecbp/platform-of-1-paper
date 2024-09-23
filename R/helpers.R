
compileSimulations = function(path) {
  
  # List all files in the given directory
  files = paste0(path,list.files(path))
  
  # Pre-allocate list for simulations
  sims = vector(mode = "list", length = length(files))
  
  for (i in 1:length(files)) {
    t = readRDS(files[i])
    sims[[i]] = t
  }
  
  bind_rows(sims) 
  
}

precompileBrmsModel = function(NUM_TRTS) {
  # Precompile Stan model for simulation repetitions
  
  # Build design matrix for
  X = diag(NUM_TRTS) 
  X[, 1] = 1
  X = X[rep(seq_len(nrow(X)), each = 1),]
  X = X[rep(seq_len(nrow(X)), times = 1),]
  Y = rnorm(nrow(X))
  DATA = cbind(X, Y = Y) |> tibble::as_tibble()
  colnames(DATA) = c(paste0("X", 1:NUM_TRTS), "Y", "period", "trt")
  
  brms::brm(formula = paste0("Y~", paste0("X", 2:NUM_TRTS, collapse = "+")),
            data = DATA,
            family = gaussian(link = "identity"),
            prior = c(
              brms::set_prior("normal(0,100)", class = "Intercept"),
              brms::set_prior("normal(0,100)", class = "b")
            ),
            refresh = 0)
  
}

precompileAR1BrmsModel = function(NUM_TRTS) {
  # Precompile Stan model for simulation repetitions
  
  # Build design matrix for
  X = diag(NUM_TRTS) 
  X[, 1] = 1
  X = X[rep(seq_len(nrow(X)), each = 1),]
  X = X[rep(seq_len(nrow(X)), times = 1),]
  Y = rnorm(nrow(X))
  DATA = cbind(X, Y = Y) |> 
    tibble::as_tibble() |> 
    mutate(
      obs = row_number()
    )
  colnames(DATA) = c(paste0("X", 1:NUM_TRTS), "Y", "obs")
  
  brms::brm(formula = paste0("Y ~ ",                             # outcome
                             paste0("X", 2:3, collapse = " + "), # treatments
                             " + ar(time = obs, p = 1)"),        # AR1 term
            data = DATA,
            family = gaussian(link = "identity"),
            prior = c(
              brms::set_prior("normal(0,100)", class = "Intercept"),
              brms::set_prior("normal(0,100)", class = "b")
            ),
            refresh = 0)
  
}

generateOutcomes = function(X, EFFECT_SIZE, OBVS_FREQ, AR_PARAM, CONFIG) {
  
  # Initialize regression coefficients
  BETAS = rep(0, CONFIG[["NUM_TRTS"]])       # Initialize vector of null effects
  BETAS[1] = CONFIG[["INTERCEPT"]]
  BETAS[2] = EFFECT_SIZE * CONFIG[["SIGMA"]] # Set second treatment as effective treatment
  
  if (AR_PARAM == 0.0) {
    
    Y = ((X %*% BETAS) |> c()) + stats::rnorm(nrow(X), 0, sd = CONFIG[["SIGMA"]])
    
  } else if (AR_PARAM > 0.0) {
    
    Y_ar = arima.sim(n = OBVS_FREQ,
                     model = list(ar = c(AR_PARAM)), 
                     sd = CONFIG[["SIGMA"]])
    Y = Y_ar + (X %*% BETAS)
    
  } else if (AR_PARAM == "SAR") {
    
    # Could be an interesting result to include as a supplemental
    model = Arima(ts(rnorm(100),freq=4), 
                   order = c(AR_PARAM, 0, 0), 
                   seasonal = c(1,1,1),
                   fixed = c(phi=0.5, theta=-0.4, Phi=0.3, Theta=-0.2))
    Y = simulate(model, nsim = OBVS_FREQ) |> c()
    
  }
  
  Y
  
}

adaptiveRandomizationV2 = function(FIT, DATA, ALG, STAB = NULL, NUM_AGENTS = NULL) {
  
  # Use standard Thompson Sampling to adaptively allocate the probabilites
  if (ALG == "TS" || ALG == "TS2" || ALG == "TOPTWO" || ALG == "TOPTWO2") {
    
    probs = ThompsonSamplingV2(FIT, DATA, objective = CONFIG[["OBJECTIVE"]], STAB)
    
  }
  
  # Use the consensus Thompson Sampling algorithm (ver. 1)
  # One consensus posterior as weighted average of shard posteriors
  if (ALG == "CTS1") {
    
    probs = ConsensusThompsonSamplingV1(FIT, DATA, NUM_AGENTS, objective = CONFIG[["OBJECTIVE"]])
    
  }
  
  if (ALG == "CTS2") {
    
    probs = ConsensusThompsonSamplingV2(FIT, DATA, NUM_AGENTS, objective = CONFIG[["OBJECTIVE"]])
    
  }
  
  probs
}

selectNextTreatment = function(FIT, DATA, ALG, OBJECTIVE) {
  
  if (OBJECTIVE == "Maximize") { m = max }
  else { m = min }
  
  # Calculate allocation probabilities via Thompson Sampling (or other algorithm)
  current_probs = adaptiveRandomizationV2(FIT, DATA, ALG = "TS") 
  
  # My original probability matching implementation of Thompson Sampling
  if (ALG == "TS") {
    
    sample(1:CONFIG[["NUM_TRTS"]], size = 1, prob = current_probs)
    
  # A new TS implementation that samples from the outcome posterior  
  } else if (ALG == "TS2") {
    
    # Update posterior with given data
    newfit = update(FIT, newdata = DATA)
    
    # Calculate and sample from outcome posterior distributions
    samples = posterior::as_draws_df(newfit) |> 
      select(contains("b_")) |> 
      rename(b_X1 = b_Intercept) |> 
      mutate(across(-b_X1, ~ . + b_X1)) |> 
      summarise(across(everything(), ~ sample(., 1))) |> 
      unlist()
    
    # Return arm with optimized, sampled reward
    which(samples == m(samples)) 
  
  # Top-Two Thompson Sampling  
  } else if (ALG == "TOPTWO") {
    
    # Update posterior with given data
    newfit = update(FIT, newdata = DATA)
    
    # Calculate and sample from outcome posterior distributions
    samples = posterior::as_draws_df(newfit) |> 
      select(contains("b_")) |> 
      rename(b_X1 = b_Intercept) |> 
      mutate(across(-b_X1, ~ . + b_X1)) |> 
      summarise(across(everything(), ~ sample(., 1))) |> 
      unlist()
    
    # Get the top two arms based on the sample outcomes
    toptwo = order(samples, decreasing = (OBJECTIVE == "Maximize"))[1:2]
    
    # Sample evenly from the top two contenders
    sample(toptwo, size = 1)
    
  } else if (ALG == "TOPTWO2") {
    
    # Update posterior with given data
    newfit = update(FIT, newdata = DATA)
    
    # Calculate and sample from outcome posterior distributions
    samples = posterior::as_draws_df(newfit) |> 
      select(contains("b_")) |> 
      rename(b_X1 = b_Intercept) |> 
      mutate(across(-b_X1, ~ . + b_X1)) |> 
      summarise(across(everything(), ~ sample(., 1))) |> 
      unlist()
    
    toptwo = order(samples, decreasing = (OBJECTIVE == "Maximize"))[1:2]
    
    # Take the probability of being optimal of each arm and use this to sample
    toptwo_probs = current_probs[toptwo]
    prob_param = max(DATA$period) / CONFIG[["MAX_DURATION"]]
    toptwo_probs = toptwo_probs^prob_param / sum(toptwo_probs^prob_param)
    
    # Sample evenly from the top two contenders
    sample(toptwo, size = 1, prob = toptwo_probs)
    
  }
  
  
}

getPosteriorSamples = function(FIT, DATA, ALG) {
  
  # Use standard Thompson Sampling to adaptively allocate the probabilites
  if (ALG == "TS"|| ALG == "TS2" || ALG == "TOPTWO" || ALG == "TOPTWO2") {
    
    newfit = update(FIT, newdata = DATA)
    
    samples = posterior::as_draws_df(newfit)
    
  }
  
  if (ALG %in% c("CTS1", "CTS2")) {
    
    adaptive_data = DATA |> filter(agent != 0)  
    
    local_posteriors = adaptive_data |> 
      as_tibble() |> 
      group_by(agent) |>  
      summarize(
        weight = n() / nrow(adaptive_data)
      ) |> 
      mutate(
        post_mat = map2(agent, weight, function(a, w) {
          
          # Get in initialization data and data assigned to agent
          shard = DATA |> filter((agent == a | agent == 0))
          
          # Generate new posterior
          new_posterior = update(FIT, newdata = shard)
          
          # Get weighted local posterior
          as_draws_matrix(new_posterior) * w
          
        })
      )
    
    consensus_posterior = Reduce(`+`, local_posteriors$post_mat)
    
    # Now run standard Thompson Sampling algorithm on consensus posterior
    samples = as.data.frame(consensus_posterior)
    
  }
  
  samples
}

