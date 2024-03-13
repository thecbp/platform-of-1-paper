TS = function(posterior, c, objective = "Maximize") {
  
  # Set the correct objective function to optimize for
  if (objective == "Maximize") { m = max }
  else { m = min }
  
  betas = as.data.frame(posterior) %>%
    dplyr::select(tidyselect::contains("b_"))
  
  n_trts = ncol(betas)
  
  # Design matrix representing outcome based on treatments
  X = diag(n_trts)
  X[,1] = 1
  colnames(X) = paste0("X", 1:n_trts)
  
  # Find which arm produced optimal outcome
  posterior_outcomes = X %*% t(betas)
  optimal_arms = apply(posterior_outcomes, 2, function(x) {which(x == m(x))})
  
  # Allocation = proportion of times a treatment was optimal in posterior sample
  probs = c(table(optimal_arms) / nrow(betas))
  
  # Stabilize the probabilities
  stabilized_probs = probs^c / sum(probs^c)
  
  # Check if any arms were not optimal in any sample, then re-add
  missing = setdiff(1:n_trts, names(stabilized_probs))
  
  out = rep(0, times = n_trts)
  names(out) = 1:n_trts
  
  for (i in 1:n_trts) {
    
    if (i %in% missing) { out[i] = 0 }
    else { out[i] = stabilized_probs[[as.character(i)]]  }
  }
  
  out
  
}

rebalance_probabilities = function(probs, current, start) { 
  
  # Return original probabilities if current treatment set matches the start
  if (all(current == start))
    return(probs)
  
  # If there's only one treatment left, just return 1
  if (length(current) == 1)
    return(1)
  
  # find missing arm
  futile_arm = setdiff(start, current)
  
  # Otherwise, rebalance the probabilities
  # Calculate the redistributions in probability such that the relative odds are invariant
  # Note to self: does not generalize to more than 2 treatments
  A = matrix(c(1, 1, (rev(probs[-futile_arm]) * c(1, -1))),
             nrow = length(current), byrow = T
  )
  b = matrix(c(probs[futile_arm], 0), nrow = length(current))
  
  deltas = c(solve(A, b))
  
  return(probs[-futile_arm] + deltas)
}