simulationTab = sidebarLayout(
  sidebarPanel(
    wellPanel(
      h3("Trial Parameters"),
      textInput(inputId = "OUTPATH",
                label = "Directory for simulations",
                value = "test"),
      numericInput(inputId = "NUM_TRTS",
                   label = "Number of treatments",
                   value = 3),
      numericInput(inputId = "OBVS_FREQ",
                   label = "Observation frequency",
                   value = 5),
      numericInput(inputId = "MAXIMUM_DURATION",
                   label = "Trial length",
                   value = 10),
      numericInput(inputId = "NUM_SIMS",
                   label = "Number of simulations",
                   value = 10),
      numericInput(inputId = "DELTA",
                   label = "Clinically significant value (delta)",
                   value = 100),
      sliderInput(inputId = "ETAE",
                  label = "Efficacy cutoff",
                  min = 0, max = 0.99, value = 0.8, step = 0.01)
    ),
    wellPanel(
      h3("Outcome Model Parameters"),
      selectInput(inputId = "OBJECTIVE",
                  label = "Maximize or minimize reward?",
                  choices = c("Maximize", "Minimize")),
      numericInput(inputId = "INTERCEPT",
                   label = "Model Intercept",
                   value = 100),
      numericInput(inputId = "EFFECT_SIZE",
                   label = "Standardized Effect Size of Optimal Arm",
                   value = 0),
      sliderInput(inputId = "AR_PARAM",
                  label = "Degree of serial correlation",
                  min = 0, max = 0.95, value = 0, step = 0.05)
    ),
    wellPanel(
      h3("Prior Parameters"),
      fluidRow(
        column(6, numericInput(inputId = "intercept_prior_mean",
                               label = paste0("Intercept Prior Mean"),
                               value = 0)),
        column(6, numericInput(inputId = "intercept_prior_variance",
                               label = paste0("Intercept Prior Variance"),
                               value = 100))
        
      ),
      fluidRow(
        column(6, numericInput(inputId = "treatment_prior_mean",
                               label = paste0("Treatment Prior Mean"),
                               value = 0)),
        column(6, numericInput(inputId = "treatment_prior_variance",
                               label = paste0("Treatment Prior Variance"),
                               value = 100))
      )
    ),
    fluidRow(
      column(12, actionButton(inputId = "simulate",
                              label = "Simulate",
                              width = "100%"))
    )
  ),
  mainPanel(
    conditionalPanel(condition = "input.simulate == 0",
                     wellPanel(
                       h2("No simulations created yet"),
                       p("Set desired parameters and press the Simulate button")
                     )
    ),
    conditionalPanel(condition = "input.simulate > 0",
                     wellPanel(
                       tableOutput("sims"),
                       plotOutput("epp"),
                       plotOutput("power")
                     )
    )
  )
)
