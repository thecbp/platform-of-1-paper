# Reducing duration in order to cut down on computation time
loadP1Config = function(NUM_TRTS) {
  
  list(
    NAME = "P1-Standard",
    NUM_TRTS = NUM_TRTS,
    MAX_DURATION = 50,
    OBJECTIVE = "Minimize",
    INTERCEPT = 130,
    SIGMA = 10,
    PRIORS = c(
      brms::set_prior("normal(0,100)", class = "Intercept"),
      brms::set_prior("normal(0,100)", class = "b")
    ), 
    COMPILED_FIT = precompileBrmsModel(NUM_TRTS),
    NUM_CHAINS = 4, 
    NUM_ITER_WARMUP = 5000,
    NUM_ITER_POST = 10000,
    THINNING_FACTOR = 2,
    ADAPT_DELTA = 0.95,
    DELTA = 130.0
  )
  
}
