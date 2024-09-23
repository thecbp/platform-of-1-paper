server <- function(input, output, session) {

  simulations = reactive({

    # Create directory if it doesn't already exist
    if (!dir.exists(input[['OUTPATH']]))
      dir.create(input[['OUTPATH']])

    withProgress(message = "Compiling Stan model...", {
      CONFIG = loadP1Config(NUM_TRTS = 3)
      CONFIG[["MAX_DURATION"]] = input[["MAXIMUM_DURATION"]] 
    })
    
    # Set up dataframe to log all iterations that need to be run
    ITERS =
      tidyr::expand_grid(
        REP_ID = seq_len(input[["NUM_SIMS"]]),
        OBVS_FREQ = input[["OBVS_FREQ"]],
        EFFECT_SIZE = input[["EFFECT_SIZE"]],
        AR_PARAM = input[["AR_PARAM"]]
      ) |>
      dplyr::mutate(
        ALG = "TS",
        AR_PARAM_STR = if_else(AR_PARAM == 0.0, "0.0", as.character(AR_PARAM)),
        FILE = glue::glue("{input[['OUTPATH']]}/{ALG}-ES{EFFECT_SIZE}-OF{OBVS_FREQ}-PHI{AR_PARAM_STR}-{REP_ID}.rds")
      )

    # Simulate trials with the specified true treatment effect
    withProgress(message = "Simulating Platform-of-1 trials...", {

      # Set up parallelization for purrr
      plan(multisession, workers = availableCores() - 1)

      future_pwalk(
        list(
          ITERS[["REP_ID"]],
          ITERS[["OBVS_FREQ"]],
          ITERS[["EFFECT_SIZE"]],
          ITERS[["AR_PARAM"]],
          ITERS[["FILE"]]
        ),
        .options = furrr_options(scheduling = 1, seed = TRUE),
        function(ri, of, es, arp, f) {

          # For running simulations without interim decisions
          s = simulateData(REP_ID = ri,
                           OBVS_FREQ = of,
                           EFFECT_SIZE = es,
                           AR_PARAM = arp,
                           CONFIG = CONFIG,
                           ALG = "TS"
          )

          saveRDS(s, file = f)

        })
    })
    
    withProgress(message = "Compiling simulated Platform-of-1 trials...", {

      sims = foreach(file = list.files(input[['OUTPATH']])) %do% {
        readRDS(paste0(input[['OUTPATH']], "/", file))
      } |>
        bind_rows()

    })
    
    return(sims)

  }) |> 
    bindEvent(input[["simulate"]])
  
  output$sims = renderTable({
    simulations()
  })
  
  output$epp = renderPlot({
    
    epp = simulations() |> 
      select(REP_ID, OBVS_FREQ, EFFECT_SIZE, AR_PARAM, contains("assigned")) |> 
      pivot_longer(contains("assigned"), 
                   names_to = "col", 
                   values_to = "assigned") |> 
      mutate(
        period = str_extract(col, "\\d+") |> as.numeric()
      ) |> 
      group_by(OBVS_FREQ, EFFECT_SIZE, AR_PARAM, period) |> 
      summarize(
        EPP = mean(assigned == 2)
      )
    
    epp |> 
      filter(period > input[["NUM_TRTS"]]) |>  
      ggplot(aes(x = period, y = EPP)) + 
      geom_line() +
      labs(
        title = "EPP as a function of period",
        y = "EPP",
        x = "Period"
      ) +
      theme_minimal() + 
      theme(legend.position = "bottom") + 
      theme(panel.border = element_rect(fill = "transparent", # Needed to add the border
                                        color = "black", linewidth = 1.0))
      
  })
  
  output$power = renderPlot({
    
    power = simulations() |> 
      select(REP_ID, OBVS_FREQ, EFFECT_SIZE, AR_PARAM, starts_with("b_X")) |> 
      pivot_longer(starts_with("b_X"), 
                   names_to = "col", 
                   values_to = "prob_eff") |> 
      mutate(
        arm = str_extract(col, "X\\d+") |> factor(levels =paste0("X", 1:input[["NUM_TRTS"]])),
        period = str_extract(col, "_\\d+") |> 
          str_replace("_", "") |> 
          as.numeric()
      ) |> 
      filter(arm == "X2", period > input[["NUM_TRTS"]]) |> 
      group_by(OBVS_FREQ, AR_PARAM, EFFECT_SIZE, period) |> 
      summarize(
        power = mean(prob_eff > input[["ETAE"]])
      )
    
    power |> 
      ggplot(aes(x = period, y = power)) + 
      geom_line() +
      geom_hline(aes(yintercept = 1 / input[["NUM_TRTS"]])) +
      labs(
        title = "Power as a function of period",
        y = "Power",
        x = "Period"
      ) +
      theme_minimal() + 
      theme(legend.position = "bottom") + 
      theme(panel.border = element_rect(fill = "transparent", # Needed to add the border
                                        color = "black", linewidth = 1.0))
    
  })


}
