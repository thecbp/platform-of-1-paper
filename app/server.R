server <- function(input, output, session) {

  ############################

  # SIMULATION FUNCTIONALITY #

  ############################

  simulations = reactive({

    if (input$simulate > 0) { 
      
      PARAM_ARGS = list(
        NUM_TRTS = input$n_trts,
        NUM_INIT_CYCLES = input$n_burn_cycles,
        NUM_OBVS_INIT = input$burn_n_obvs_per_period,
        NUM_OBVS_ADAPTIVE = input$adaptive_n_obvs_per_period,
        MAX_DURATION = input$maximum_duration,
        EFFECT_SIZE = input$effect_size, 
        AR1_PHI = input$serial_correlation,
        OUT_PATH = input$out_path
      )
    
      # Create directory if it doesn't already exist
      if (!dir.exists(PARAM_ARGS[["OUT_PATH"]]))
        dir.create(PARAM_ARGS[["OUT_PATH"]])
      
      prior_text = list(
        "Intercept" = paste0("normal(", input$intercept_prior_mean,
                             ",", input$intercept_prior_variance, ")"),
        "b" = paste0("normal(", input$treatment_prior_mean,
                     ",", input$treatment_prior_variance, ")")
      )
      
      MODEL_ARGS = list(
        OBJECTIVE = input$objective, 
        INTERCEPT = input$intercept, 
        SIGMA = input$within_person_noise,
        DELTA = input$delta,
        PRIORS = c(
          brms::set_prior(prior_text[["Intercept"]], class = "Intercept"),
          brms::set_prior(prior_text[["b"]], class = "b")
        ),
        NUM_CHAINS = input$n_chains,
        NUM_ITER_WARMUP = input$warm_up,
        NUM_ITER_POST = input$samples_per_chain,
        THINNING_FACTOR = 2, 
        ADAPT_DELTA = input$adapt_delta,
        MAX_TREE_DEPTH = input$max_treedepth
      )
      
      withProgress(message = "Pre-compiling TS Model...", {
        
        # Pre-build Stan model to speed up simulations
        X = diag(input$n_trts)
        X[, 1] = 1
        X = X[rep(seq_len(nrow(X)), each = 1),]
        X = X[rep(seq_len(nrow(X)), times = 1),]
        Y = rnorm(nrow(X))
        DATA = cbind(X, Y = Y) %>% tibble::as_tibble()
        colnames(DATA) = c(paste0("X", 1:input$n_trts), "Y", "period", "trt")
        
        # Store the compiled model for use in simulations
        MODEL_ARGS[["COMPILED_FIT"]] = brms::brm(formula = paste0("Y ~", paste0("X", 2:input$n_trts, collapse = "+")),
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
      })
      
      sims = vector('list', length = input$n_sims)

      # Simulate trials with the specified true treatment effect
      withProgress(message = "Simulating Platform-of-1 trials...", {

        for (i in seq_len(input$n_sims)) {
  
          PARAM_ARGS[["REP_ID"]] = i
          PARAM_ARGS[["FILE"]] = glue::glue("sim-{i}.rds")
          
          sim = simulatePlatformOf1Shiny(PARAM_ARGS, MODEL_ARGS)
          
          sims[[i]] = sim

          # Increment progress on the progress bar
          incProgress(1 / input$n_sims)
        }

      }) 
      
      sims = sims %>% bind_rows()

      return(sims)

    }
  })

  optimal_arm = reactive({

    treatment_effect_input_names = paste0("trt-effect-", 1:input$n_trts)
    effects = c()

    # Add mean changes to intercept to get total effect for non-reference arms
    for (i in 1:input$n_trts) {
      if (i == 1)  {
        effects = c(effects, input[[treatment_effect_input_names[i]]])
      } else {
        effects = c(effects, input[[treatment_effect_input_names[1]]] + input[[treatment_effect_input_names[i]]])
      }
    }

    if (input$objective == "Maximize") {
      optimal_index = which(effects == max(effects))
    } else {
      optimal_index = which(effects == min(effects))
    }

    # Give back the name of the optimal arm
    if (length(optimal_index) > 1) {
      paste0("X1")
    } else {
      paste0("X", optimal_index)
    }

  })

  howto = reactive({
    readLines("howto.txt")
  })
  
  output$simulations = renderTable({
    simulations() 
  })
  
  output$epp = renderPlot({
    
    epp = simulations() %>% 
      select(EFFECT_SIZE, AR1_PHI, NUM_OBVS_INIT, contains("assigned")) %>% 
      pivot_longer(contains("assigned"), names_to = "col", values_to = "assigned") %>% 
      mutate(
        period = str_extract(col, "\\d+") %>% as.numeric,
        ES = factor(EFFECT_SIZE, levels = c(0, -2, -5, -10), 
                    labels = c("Null", "Small", "Moderate", "Large")),
        PHI = factor(AR1_PHI, levels = c(0, 25, 50, 75),
                     labels = c("0%", "25%", "50%", "75%")),
      ) %>% 
      group_by(period, ES, PHI, NUM_OBVS_INIT) %>% 
      summarize(n = n(), EPP = mean(assigned == 2)) 
    
    p  = epp %>% 
      ggplot(aes(x = period, y = EPP, color = ES, group = ES)) + 
      geom_hline(aes(yintercept = 1/3), alpha = 0.5, linetype = "dotted") +
      geom_line() +
      theme_minimal() + 
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5)) +
      labs(
        title = "Expected proportion of periods on optimal treatment by effect size"
      ) + scale_color_brewer(palette = "Dark2")
    
    p
    
  })
  
  output$fwer_power = renderPlot({
    
    pwr = simulations() %>% 
      select(EFFECT_SIZE, AR1_PHI, NUM_OBVS_INIT, contains("interval2")) %>% 
      pivot_longer(contains("interval2"), names_to = "col", values_to = "prob") %>% 
      mutate( 
        period = str_extract(col, "_\\d+") %>% str_replace("_", "") %>%  as.numeric,
        treatment = str_extract(col, "X\\d+")
      ) %>% 
      group_by(period, EFFECT_SIZE, AR1_PHI, NUM_OBVS_INIT, treatment) %>% 
      summarize(
        pwr_90 = mean(prob > 0.90),
        pwr_95 = mean(prob > 0.95),
        pwr_975 = mean(prob > 0.975)
      ) %>% 
      pivot_longer(contains("pwr"),
                   names_to = "gamma",
                   values_to = "power")
    
    p = pwr %>% 
      filter(treatment == "X2") %>% 
      ggplot(aes(x = period, y = power, color = gamma)) +
      geom_line(linewidth = 1) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_color_brewer(
        palette = "Dark2",
        labels = c("0.9", "0.95", "0.975"),
        name = "Graduation Threshold") +
      labs(
        x = "Treatment Period",
        y = "Power"
      ) 
    
    p
    
  })
  

  # Text for helping guide the user to use the app
  output$step1 = renderUI({ withMathJax(howto()[1]) })
  output$step2 = renderUI({ withMathJax(howto()[2]) })
  output$step3 = renderUI({ withMathJax(howto()[3]) })
  output$step4 = renderUI({ withMathJax(howto()[4]) })

}
