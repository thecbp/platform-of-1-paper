ThompsonSamplingV2 = function(FIT, DATA, objective = "Maximize", STAB = NULL) {
  
  # Set the correct objective function to optimize for
  if (objective == "Maximize") { m = max }
  else { m = min }
  
  updated_fit = update(FIT, newdata = DATA)
  
  betas = as_draws_df(updated_fit) %>%
    dplyr::select(tidyselect::contains("b_"))
  
  n_trts = ncol(betas)
  
  # Design matrix representing outcome based on treatments
  X = diag(n_trts)
  X[,1] = 1
  colnames(X) = paste0("X", 1:n_trts)
  
  # Find which arm produced optimal outcome
  posterior_outcomes = X %*% t(betas)
  optimal_arms = apply(posterior_outcomes, 2, function(x) { which(x == m(x)) })
  
  # Allocation = proportion of times a treatment was optimal in posterior sample
  probs = c(table(optimal_arms) / nrow(betas))
  
  if (!is.null(STAB)) {
    probs = probs^STAB / sum(probs^STAB)
  }
  
  # Check if any arms were not optimal in any sample, then re-add
  missing = setdiff(1:n_trts, names(probs))
  
  out = rep(0, times = n_trts)
  names(out) = 1:n_trts
  
  for (i in 1:n_trts) {
    
    if (i %in% missing) { out[i] = 0 }
    else { out[i] = probs[[as.character(i)]]  }
  }
  
  out
  
}
