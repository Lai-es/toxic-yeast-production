# ============================================================
# Yeast-Toxin Model: Shiny App
#
# Model:
#   dY/dt = r*Y - d*P*Y        (linear growth, toxin-induced death)
#   dP/dt = g*Y - z*P          (production by Y, decay)
#
# Nontrivial equilibrium:
#   Y* = r*z / (d*g)
#   P* = r / d
#
# Required packages: shiny, deSolve, ggplot2
# install.packages(c("shiny", "deSolve", "ggplot2"))
# ============================================================

library(shiny)
library(deSolve)
library(ggplot2)

# ============================================================
# Model functions (pure, no reactive context — safe to call
# from the console for testing)
# ============================================================

# ODE system passed to deSolve::ode()
yeast_toxin_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dY <- r * Y - d * P * Y
    dP <- g * Y - z * P
    list(c(dY, dP))
  })
}

# Closed-form nontrivial equilibrium
# Returns a list with Y_star and P_star
compute_equilibrium <- function(r, g, d, z) {
  list(
    Y_star = (r * z) / (d * g),
    P_star = r / d
  )
}

# ============================================================
# Equilibrium sweep helper
#
# Computes Y* and P* across a sequence of values for one
# chosen parameter, holding all others fixed.
#
# Arguments:
#   param_name  : one of "r", "g", "d", "z" (string)
#   param_range : numeric vector of values to sweep over
#   fixed       : named list with the four fixed parameter
#                 values, e.g. list(r=0.5, g=0.5, d=0.5, z=0.5)
#
# Returns a data frame with columns:
#   param_value, Y_star, P_star
# ============================================================
sweep_equilibrium <- function(param_name, param_range, fixed) {
  results <- lapply(param_range, function(val) {
    # Build the full parameter set by overriding the swept parameter
    params <- fixed
    params[[param_name]] <- val
    eq <- compute_equilibrium(r = params$r, g = params$g,
                              d = params$d, z = params$z)
    data.frame(param_value = val,
               Y_star      = eq$Y_star,
               P_star      = eq$P_star)
  })
  do.call(rbind, results)
}

# ============================================================
# Plot helper: render one equilibrium-vs-parameter curve
#
# Arguments:
#   sweep_df    : data frame from sweep_equilibrium()
#   y_col       : "Y_star" or "P_star"
#   y_label     : y-axis label string
#   title       : plot title string
#   line_color  : hex color for the curve and highlight point
#   current_x   : current slider value (draws vline + point)
#   current_y   : equilibrium value at current_x
#   x_label     : x-axis label string
# ============================================================
plot_equilibrium_curve <- function(sweep_df, y_col, y_label, title,
                                   line_color, current_x, current_y,
                                   x_label) {
  x_max <- max(sweep_df$param_value) * 1.03
  y_max <- max(sweep_df[[y_col]])    * 1.03
  
  ggplot(sweep_df, aes(x = param_value, y = .data[[y_col]])) +
    geom_line(color = line_color, linewidth = 1) +
    geom_vline(xintercept = current_x, linetype = "dashed",
               color = "gray40", linewidth = 0.6) +
    geom_point(data = data.frame(param_value = current_x,
                                 val         = current_y),
               aes(x = param_value, y = val),
               color = line_color, size = 3, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max), expand = FALSE) +
    labs(x = x_label, y = y_label, title = title) +
    theme_minimal(base_size = 14)
}

# ============================================================
# Shared sweep range for all parameters
# All sliders currently share [0.01, 1], so one set of
# constants covers all four panels.
# ============================================================
SWEEP_MIN <- 0.01
SWEEP_MAX <- 1
SWEEP_N   <- 200

# Caption style reused across all tab descriptions
CAPTION_STYLE <- paste0(
  "color: #888888;",
  "margin-left: 40px;",
  "margin-right: 40px;",
  "margin-top: 10px;"
)

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  titlePanel("Yeast-Toxin Model: Linear Growth with Self-Produced Toxicity"),
  
  sidebarLayout(
    
    # ------ Sidebar: diagram + parameter sliders + ICs ------
    sidebarPanel(
      h4("Model Diagram"),
      tags$img(src = "model_diagram.png",
               style = "max-width: 100%; height: auto;"),
      
      hr(),
      h4("Model Parameters"),
      sliderInput("r", "Replication rate (r):",
                  min = 0.01, max = 1, value = 0.5, step = 0.01),
      sliderInput("g", "Product generation rate (g):",
                  min = 0.01, max = 1, value = 0.5, step = 0.01),
      sliderInput("d", "Yeast death rate (d):",
                  min = 0.01, max = 1, value = 0.5, step = 0.01),
      sliderInput("z", "Product decay rate (z):",
                  min = 0.01, max = 1, value = 0.5, step = 0.01),
      
      hr(),
      h4("Initial Conditions"),
      numericInput("Y0", "Initial yeast population:", value = 10, min = 0),
      numericInput("P0", "Initial product amount:",   value = 0,  min = 0)
    ),
    
    # ------ Main panel: equilibrium display + tab plots -----
    mainPanel(
      
      # Always-visible equilibrium values
      wellPanel(
        h4("Current Equilibrium"),
        fluidRow(
          column(6, verbatimTextOutput("eqYText")),
          column(6, verbatimTextOutput("eqPText"))
        )
      ),
      
      tabsetPanel(
        
        # ---- Tab 1: Time series ----
        tabPanel(
          "Time Series",
          plotOutput("timeSeriesPlot", height = "500px"),
          p("Solid lines: simulated trajectories. Dashed lines: predicted",
            "equilibrium Y* and P* for the current parameter set,",
            "shown as a convergence check.",
            style = CAPTION_STYLE),
          sliderInput("tmax", "Simulation time:",
                      min = 10, max = 300, value = 150,
                      step = 10, width = "100%")
        ),
        
        # ---- Tab 2: Phase plane ----
        tabPanel(
          "Phase Plane",
          plotOutput("phasePlanePlot", height = "500px"),
          p("Trajectory of (Y(t), P(t)) for the current parameters.",
            "Triangle marks the initial condition,",
            "star marks the predicted equilibrium.",
            style = CAPTION_STYLE)
        ),
        
        # ---- Tab 3: Equilibria vs Parameters (nested tabs) ----
        tabPanel(
          "Equilibria vs Parameters",
          tabsetPanel(
            
            tabPanel(
              "Product generation rate (g)",
              br(),
              fluidRow(
                column(6, plotOutput("eqYPlotG", height = "450px")),
                column(6, plotOutput("eqPPlotG", height = "450px"))
              ),
              p("Closed-form equilibrium Y* and P* swept across a fixed range of g,",
                "with r, d, z held at their current slider values.",
                "The vertical dashed line and point mark the current g.",
                style = CAPTION_STYLE)
            ),
            
            tabPanel(
              "Replication rate (r)",
              br(),
              fluidRow(
                column(6, plotOutput("eqYPlotR", height = "450px")),
                column(6, plotOutput("eqPPlotR", height = "450px"))
              ),
              p("Closed-form equilibrium Y* and P* swept across a fixed range of r,",
                "with g, d, z held at their current slider values.",
                "The vertical dashed line and point mark the current r.",
                style = CAPTION_STYLE)
            ),
            
            tabPanel(
              "Yeast death rate (d)",
              br(),
              fluidRow(
                column(6, plotOutput("eqYPlotD", height = "450px")),
                column(6, plotOutput("eqPPlotD", height = "450px"))
              ),
              p("Closed-form equilibrium Y* and P* swept across a fixed range of d,",
                "with r, g, z held at their current slider values.",
                "The vertical dashed line and point mark the current d.",
                style = CAPTION_STYLE)
            ),
            
            tabPanel(
              "Product decay rate (z)",
              br(),
              fluidRow(
                column(6, plotOutput("eqYPlotZ", height = "450px")),
                column(6, plotOutput("eqPPlotZ", height = "450px"))
              ),
              p("Closed-form equilibrium Y* and P* swept across a fixed range of z,",
                "with r, g, d held at their current slider values.",
                "The vertical dashed line and point mark the current z.",
                style = CAPTION_STYLE)
            )
            
          ) # end inner tabsetPanel
        )   # end "Equilibria vs Parameters" tabPanel
        
      ) # end outer tabsetPanel
    )   # end mainPanel
  )     # end sidebarLayout
)       # end fluidPage

# ============================================================
# Server
# ============================================================

server <- function(input, output, session) {
  
  # -- Reactive: ODE simulation (re-runs on any parameter change) --
  simResult <- reactive({
    state      <- c(Y = input$Y0, P = input$P0)
    parameters <- c(r = input$r, g = input$g, d = input$d, z = input$z)
    times      <- seq(0, input$tmax, length.out = 500)
    out <- ode(y = state, times = times,
               func = yeast_toxin_model, parms = parameters)
    as.data.frame(out)
  })
  
  # -- Reactive: closed-form equilibrium at current slider values --
  currentEq <- reactive({
    compute_equilibrium(r = input$r, g = input$g,
                        d = input$d, z = input$z)
  })
  
  # -- Reactive: fixed parameters used inside each sweep
  #    (the four sliders minus whichever is being swept) --
  fixedParams <- reactive({
    list(r = input$r, g = input$g, d = input$d, z = input$z)
  })
  
  # ---- Equilibrium display (wellPanel) ----
  output$eqYText <- renderText({
    paste0("Y* = ", round(currentEq()$Y_star, 3))
  })
  output$eqPText <- renderText({
    paste0("P* = ", round(currentEq()$P_star, 3))
  })
  
  # ============================================================
  # Time series plot
  # ============================================================
  output$timeSeriesPlot <- renderPlot({
    df <- simResult()
    eq <- currentEq()
    
    ggplot(df, aes(x = time)) +
      geom_line(aes(y = Y, color = "Yeast cells (Y)"),   linewidth = 1) +
      geom_line(aes(y = P, color = "Product amount (P)"), linewidth = 1) +
      geom_hline(yintercept = eq$Y_star, linetype = "dashed",
                 color = "#1b9e77", linewidth = 0.7) +
      geom_hline(yintercept = eq$P_star, linetype = "dashed",
                 color = "#d95f02", linewidth = 0.7) +
      coord_cartesian(xlim = c(0, max(df$time) * 1.01), ylim = c(0, NA), expand = FALSE) +
      scale_color_manual(
        values = c("Yeast cells (Y)" = "#1b9e77",
                   "Product amount (P)" = "#d95f02")
      ) +
      labs(x = "Time", y = "Yeast population / Product amount",
           color = NULL,
           title = "Time Series with Predicted Equilibria (dashed)") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
  })
  
  # ============================================================
  # Phase plane plot
  # ============================================================
  output$phasePlanePlot <- renderPlot({
    df <- simResult()
    eq <- currentEq()
    
    ggplot(df, aes(x = Y, y = P)) +
      geom_path(color = "#7570b3", linewidth = 1,
                arrow = arrow(length = unit(0.2, "cm"),
                              ends = "last", type = "closed")) +
      geom_point(data = df[1, ], aes(x = Y, y = P),
                 color = "black", size = 3, shape = 17) +
      geom_point(aes(x = eq$Y_star, y = eq$P_star),
                 color = "red", size = 4, shape = 8) +
      coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = FALSE) +
      labs(x = "Yeast population (Y)", y = "Product amount (P)",
           title = "Phase Plane: Trajectory (triangle = start, star = equilibrium)") +
      theme_minimal(base_size = 14)
  })
  
  # ============================================================
  # Equilibria vs Parameters: one pair of renderPlot blocks per
  # parameter, all using sweep_equilibrium() + plot_equilibrium_curve()
  # ============================================================
  
  # Helper that builds the sweep data frame for the currently
  # active tab, so we don't repeat the seq() + sweep call inline
  make_sweep <- function(param_name) {
    sweep_equilibrium(
      param_name  = param_name,
      param_range = seq(SWEEP_MIN, SWEEP_MAX, length.out = SWEEP_N),
      fixed       = fixedParams()
    )
  }
  
  # ---- g sweep ----
  output$eqYPlotG <- renderPlot({
    df <- make_sweep("g")
    plot_equilibrium_curve(df,
                           y_col     = "Y_star",
                           y_label   = "Y* (equilibrium yeast population)",
                           title     = "Yeast Equilibrium vs Product Generation Rate",
                           line_color = "#1b9e77",
                           current_x = input$g,
                           current_y = currentEq()$Y_star,
                           x_label   = "Product generation rate (g)"
    )
  })
  output$eqPPlotG <- renderPlot({
    df <- make_sweep("g")
    plot_equilibrium_curve(df,
                           y_col     = "P_star",
                           y_label   = "P* (equilibrium product amount)",
                           title     = "Product Equilibrium vs Product Generation Rate",
                           line_color = "#d95f02",
                           current_x = input$g,
                           current_y = currentEq()$P_star,
                           x_label   = "Product generation rate (g)"
    )
  })
  
  # ---- r sweep ----
  output$eqYPlotR <- renderPlot({
    df <- make_sweep("r")
    plot_equilibrium_curve(df,
                           y_col     = "Y_star",
                           y_label   = "Y* (equilibrium yeast population)",
                           title     = "Yeast Equilibrium vs Replication Rate",
                           line_color = "#1b9e77",
                           current_x = input$r,
                           current_y = currentEq()$Y_star,
                           x_label   = "Replication rate (r)"
    )
  })
  output$eqPPlotR <- renderPlot({
    df <- make_sweep("r")
    plot_equilibrium_curve(df,
                           y_col     = "P_star",
                           y_label   = "P* (equilibrium product amount)",
                           title     = "Product Equilibrium vs Replication Rate",
                           line_color = "#d95f02",
                           current_x = input$r,
                           current_y = currentEq()$P_star,
                           x_label   = "Replication rate (r)"
    )
  })
  
  # ---- d sweep ----
  output$eqYPlotD <- renderPlot({
    df <- make_sweep("d")
    plot_equilibrium_curve(df,
                           y_col     = "Y_star",
                           y_label   = "Y* (equilibrium yeast population)",
                           title     = "Yeast Equilibrium vs Yeast Death Rate",
                           line_color = "#1b9e77",
                           current_x = input$d,
                           current_y = currentEq()$Y_star,
                           x_label   = "Yeast death rate (d)"
    )
  })
  output$eqPPlotD <- renderPlot({
    df <- make_sweep("d")
    plot_equilibrium_curve(df,
                           y_col     = "P_star",
                           y_label   = "P* (equilibrium product amount)",
                           title     = "Product Equilibrium vs Yeast Death Rate",
                           line_color = "#d95f02",
                           current_x = input$d,
                           current_y = currentEq()$P_star,
                           x_label   = "Yeast death rate (d)"
    )
  })
  
  # ---- z sweep ----
  output$eqYPlotZ <- renderPlot({
    df <- make_sweep("z")
    plot_equilibrium_curve(df,
                           y_col     = "Y_star",
                           y_label   = "Y* (equilibrium yeast population)",
                           title     = "Yeast Equilibrium vs Product Decay Rate",
                           line_color = "#1b9e77",
                           current_x = input$z,
                           current_y = currentEq()$Y_star,
                           x_label   = "Product decay rate (z)"
    )
  })
  output$eqPPlotZ <- renderPlot({
    df <- make_sweep("z")
    plot_equilibrium_curve(df,
                           y_col     = "P_star",
                           y_label   = "P* (equilibrium product amount)",
                           title     = "Product Equilibrium vs Product Decay Rate",
                           line_color = "#d95f02",
                           current_x = input$z,
                           current_y = currentEq()$P_star,
                           x_label   = "Product decay rate (z)"
    )
  })
  
}

shinyApp(ui = ui, server = server)