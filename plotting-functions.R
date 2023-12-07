Figure2 = function() {
  
  # Describes how EPP changes by treatment period (1 obs per period) at 0 auto-correlation
  
  epp %>% 
    filter(PHI == "0%", treatment == "X2") %>% 
    ggplot(aes(x = period, y = EPP, color = ES, group = ES)) + 
    geom_hline(aes(yintercept = 1/3), alpha = 0.5, linetype = "dotted") +
    geom_line() +
    theme_minimal() + 
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) +
    labs(
      title = "Expected proportion of periods on optimal treatment by effect size"
    )
  
}

Figure2b = function() {
  
  # Describes how EPP changes by treatment period (1 obs per period) 
  # and looks at the effect of auto-correlation on it
  
  epp %>% 
    filter(treatment == "X2", ES != "Null") %>% 
    ggplot(aes(x = period, y = EPP, color = PHI, group = PHI)) + 
    geom_hline(aes(yintercept = 1/3), alpha = 0.5, linetype = "dotted") +
    geom_line() +
    theme_minimal() + 
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ES) +
    labs(
      title = "Expected proportion of periods on optimal treatment by\neffect size and auto-correlation"
    )
  
}
