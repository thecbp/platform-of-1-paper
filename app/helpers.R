# Generates posterior samples for each sample size
generate_trial_posteriors = function(sim, i) {

  rep_posteriors = list()

  burnin_cutoff = sim[["trial_params"]][["n_trts"]] *
    sim[["trial_params"]][["n_burn_cycles"]]

  burnin_model = sim[["burnin_posteriors"]]

  posterior = as.data.frame(burnin_model)
  posterior$period = burnin_cutoff
  posterior$sim = i

  rep_posteriors[[burnin_cutoff]] = posterior

  for (p in (burnin_cutoff+1):sim[["trial_params"]][["max_duration"]]) {
    current_data = sim[["data"]] %>% filter(period <= p)
    current_posteriors = update(burnin_model, newdata = current_data)
    posterior = as.data.frame(current_posteriors)
    posterior$period = p
    rep_posteriors[[p]] = posterior
  }

  rep_posteriors
}

create_power_plot = function(decisions, nsims, duration, ntrts, optimal, fwer = F) {

  gammas_E = seq(0.50, 0.99, 0.01) # efficacy gammas

  # Strings to categorize futlity + efficacy decisions by treatment arm
  decision_strings = c(
    paste0("X", 2:n_trts, "_fut"),
    paste0("X", 2:n_trts, "_eff")
  )

  data = tibble(
    sim = 1:nsims
  ) %>%
    mutate(
      decisions = map(sim, function(i) {
        decisions[[i]]
        }),
      decision_tibble = map(decisions, function(d) {
        tibble(
          gammas = gammas_E,
          decisions = d
        ) %>%
          mutate(
            decision_type = map(decisions, function(d) { decision_strings })
          ) %>%
          unnest(cols = c("decisions")) %>%
          mutate(
            period = rep(1:duration, times = length(gammas_E))
          ) %>%
          unnest(cols = c("decisions", "decision_type"))
      })
    ) %>%
    select(sim, decision_tibble) %>%
    unnest(c("decision_tibble")) %>%
    drop_na() %>%
    mutate(decisions = unlist(decisions)) %>%
    pivot_wider(names_from = "decision_type",
                values_from = "decisions") %>%
    mutate(
      power = !!sym(paste0(optimal, "_eff"))
    ) %>%
    rowwise() %>%
    mutate(
      fwe = any(c_across(ends_with("_eff"))),
    ) %>%
    ungroup()

  if (fwer) {

    data %>%
      group_by(gammas, period) %>%
      summarize(
        FWER = mean(fwe)
      ) %>%
      ggplot(aes(x = period, y = gammas, fill = FWER)) +
      geom_tile() +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

  } else {

    data %>%
      group_by(gammas, period) %>%
      summarize(
        Power = mean(power)
      ) %>%
      ggplot(aes(x = period, y = gammas, fill = Power)) +
      geom_tile() +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

  }



}
