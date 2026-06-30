# ============================================================
# Yeast-Toxin Model: Shiny App
#
# Model:
#   dY/dt = r*Y*(1 - Y/K) - d*P*Y      (logistic growth, toxin-induced death)
#   dP/dt = g*Y - z*P                  (production by Y, decay)
#
# Nontrivial equilibrium:
#   Y* = r*K*z / (r*z + d*g*K)
#   P* = (g/z) * Y*
#
# Required packages: shiny, deSolve, ggplot2
# ============================================================

library(shiny)
library(deSolve)
library(ggplot2)

# ------------------------------------------------------------
# Model definitions (kept outside server: pure functions, no
# reactive context needed)
# ------------------------------------------------------------

# ODE system for deSolve::ode()
yeast_toxin_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dY <- r * Y * (1 - Y / K) - d * P * Y
    dP <- g * Y - z * P
    list(c(dY, dP))
  })
}

# Closed-form nontrivial equilibrium
compute_equilibrium <- function(r, K, g, d, z) {
  Y_star <- (r * K * z) / (r * z + d * g * K)
  P_star <- (g / z) * Y_star
  list(Y_star = Y_star, P_star = P_star)
}

# Fixed range used for the g-sweep equilibrium plots
G_SWEEP_MIN <- 0
G_SWEEP_MAX <- 1
G_SWEEP_N   <- 200

# ------------------------------------------------------------
# UI
# ------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Yeast-Toxin Model: Logistic Growth with Self-Produced Toxicity"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Model Parameters"),
      sliderInput("r", "Growth rate (r):", min = 0.01, max = 2, value = 0.5, step = 0.01),
      sliderInput("K", "Carrying capacity (K):", min = 10, max = 500, value = 100, step = 10),
      sliderInput("g", "Production rate (g):", min = 0, max = 2, value = 0.05, step = 0.01),
      sliderInput("d", "Toxicity susceptability (d):", min = 0, max = 1, value = 0.1, step = 0.001),
      sliderInput("z", "Product decay percentage (z):", min = 0, max = 1, value = 0.1, step = 0.01),
      
      hr(),
      h4("Initial Conditions"),
      numericInput("Y0", "Initial yeast population (Y0):", value = 10, min = 0),
      numericInput("P0", "Initial product concentration (P0):", value = 0, min = 0),
      
      hr(),
      h4("Simulation Settings"),
      sliderInput("tmax", "Simulation time:", min = 10, max = 500, value = 150, step = 10),
      
      hr(),
      h4("Current Equilibrium"),
      verbatimTextOutput("eqText")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
        "Time Series",
        plotOutput("timeSeriesPlot", height = "500px"),
        p(
          "Solid lines: simulated trajectories. Dashed lines: predicted equilibrium ",
          "Y* and P* for the current parameter set, shown as a convergence check.",
          style = "color: #888888; margin-left: 40px; margin-right: 40px; margin-top: 10px;"
        )
      ),
        tabPanel(
          "Equilibrium vs g",
          fluidRow(
            column(6, plotOutput("eqYPlot", height = "450px")),
            column(6, plotOutput("eqPPlot", height = "450px"))
          ),
          p("Closed-form equilibrium Y* and P* swept across a fixed range of g, with ",
            "r, K, d, z held at their current slider values. The vertical dashed line ",
            "and point mark the current g.",
            style = "color: #888888; margin-left: 40px; margin-right: 40px; margin-top: 10px;")
        ),
        tabPanel(
          "Phase Plane",
          plotOutput("phasePlanePlot", height = "500px"),
          p("Trajectory of (Y(t), P(t)) for the current parameters. Triangle marks the ",
            "initial condition, star marks the predicted equilibrium.",
            style = "color: #888888; margin-left: 40px; margin-right: 40px; margin-top: 10px;"),
        )
      ),
      hr(),
      h4("Model Diagram"),
      tags$img(src = "model_diagram.png", style = "max-width: 70%; height: auto;")
    )
  )
)

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------

server <- function(input, output, session) {
  
  # Reactive: run the ODE simulation whenever parameters/ICs/tmax change
  simResult <- reactive({
    state <- c(Y = input$Y0, P = input$P0)
    parameters <- c(r = input$r, K = input$K, g = input$g, d = input$d, z = input$z)
    times <- seq(0, input$tmax, length.out = 500)
    
    out <- ode(y = state, times = times, func = yeast_toxin_model, parms = parameters)
    as.data.frame(out)
  })
  
  # Reactive: equilibrium at the current slider values
  currentEq <- reactive({
    compute_equilibrium(r = input$r, K = input$K, g = input$g, d = input$d, z = input$z)
  })
  
  output$eqText <- renderText({
    eq <- currentEq()
    paste0("Y* = ", round(eq$Y_star, 2), "\nP* = ", round(eq$P_star, 2))
  })
  
  # --- Time series plot, with equilibrium convergence check ---
  output$timeSeriesPlot <- renderPlot({
    df <- simResult()
    eq <- currentEq()
    
    ggplot(df, aes(x = time)) +
      geom_line(aes(y = Y, color = "Yeast cells (Y)"), linewidth = 1) +
      geom_line(aes(y = P, color = "Product (P)"), linewidth = 1) +
      geom_hline(yintercept = eq$Y_star, linetype = "dashed", color = "#1b9e77", linewidth = 0.7) +
      geom_hline(yintercept = eq$P_star, linetype = "dashed", color = "#d95f02", linewidth = 0.7) +
      coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
      scale_color_manual(values = c("Yeast cells (Y)" = "#1b9e77", "Product (P)" = "#d95f02")) +
      labs(
        x = "Time", y = "Population / Concentration", color = NULL,
        title = "Time Series with Predicted Equilibria (dashed)"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
  })
  
  # --- Equilibrium vs g: Y* panel ---
  output$eqYPlot <- renderPlot({
    g_seq <- seq(G_SWEEP_MIN, G_SWEEP_MAX, length.out = G_SWEEP_N)
    eq_df <- data.frame(
      g = g_seq,
      Y_star = sapply(g_seq, function(gg) {
        compute_equilibrium(input$r, input$K, gg, input$d, input$z)$Y_star
      })
    )
    
    ggplot(eq_df, aes(x = g, y = Y_star)) +
      geom_line(color = "#1b9e77", linewidth = 1) +
      geom_vline(xintercept = input$g, linetype = "dashed", color = "gray40") +
      coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
      geom_point(
        data = data.frame(g = input$g, Y_star = currentEq()$Y_star),
        aes(x = g, y = Y_star), color = "#1b9e77", size = 3
      ) +
      labs(
        x = "Production rate (g)", y = "Y* (equilibrium yeast population)",
        title = "Yeast Equilibrium vs Production Rate"
      ) +
      theme_minimal(base_size = 14)
  })
  
  # --- Equilibrium vs g: P* panel ---
  output$eqPPlot <- renderPlot({
    g_seq <- seq(G_SWEEP_MIN, G_SWEEP_MAX, length.out = G_SWEEP_N)
    eq_df <- data.frame(
      g = g_seq,
      P_star = sapply(g_seq, function(gg) {
        compute_equilibrium(input$r, input$K, gg, input$d, input$z)$P_star
      })
    )
    
    ggplot(eq_df, aes(x = g, y = P_star)) +
      geom_line(color = "#d95f02", linewidth = 1) +
      geom_vline(xintercept = input$g, linetype = "dashed", color = "gray40") +
      coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
      geom_point(
        data = data.frame(g = input$g, P_star = currentEq()$P_star),
        aes(x = g, y = P_star), color = "#d95f02", size = 3
      ) +
      labs(
        x = "Production rate (g)", y = "P* (equilibrium toxin concentration)",
        title = "Toxin Equilibrium vs Production Rate"
      ) +
      theme_minimal(base_size = 14)
  })
  
  # --- Phase plane plot ---
  output$phasePlanePlot <- renderPlot({
    df <- simResult()
    eq <- currentEq()
    
    ggplot(df, aes(x = Y, y = P)) +
      geom_path(
        color = "#7570b3", linewidth = 1,
        arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")
      ) +
      geom_point(data = df[1, ], aes(x = Y, y = P), color = "black", size = 3, shape = 17) +
      geom_point(aes(x = eq$Y_star, y = eq$P_star), color = "red", size = 4, shape = 8) +
      coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
      labs(
        x = "Yeast population (Y)", y = "Toxin concentration (P)",
        title = "Phase Plane: Trajectory (triangle = start, star = equilibrium)"
      ) +
      theme_minimal(base_size = 14)
  })
}

shinyApp(ui = ui, server = server)