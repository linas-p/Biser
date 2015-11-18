library(shiny)

# Define UI for slider demo application
shinyUI(fluidPage(

#  Application title
            titlePanel("Sliders"),

# Sidebar with sliders that demonstrate various available
# options
            sidebarLayout(
                sidebarPanel(
# Simple integer interval
                    actionButton("go", "Calculate"),
                    numericInput("dt", "dt:", 1e-5),

                    numericInput("n", "Mesh on each interval:", 10),
                    numericInput("vmax_1", "Vmax_1:", 4e-5),
                    numericInput("vmax_2", "Vmax_2:", 1e-5),
                    numericInput("km_1", "Km_1:", 8e-5),
                    numericInput("km_2", "Km_2:", 2e-5),
                    numericInput("k1", "K1:", 4e-3),
                    numericInput("k2", "K2:", 2e-3),

                    numericInput("d", "Interval length(cm):", 1),
                    sliderInput("decimal", "Percentage of first layer:",
                                min = 1, max = 99, value = 10, step= 1),

                    numericInput("D_gr", "laeyr1 D_p:", 2.1e-6),
                    numericInput("D_pr", "laeyr1 D_g:", 2.1e-6),
                    numericInput("D_o2r", "laeyr1 D_o2:", 0.66e-5),

                    numericInput("D_gg", "laeyr2 D_p:", 6.3e-6),
                    numericInput("D_pg", "laeyr2 D_g:", 6.3e-6),
                    numericInput("D_o2g", "laeyr2 D_o2:", 2e-5),

                    numericInput("p_0", "p_0:", 0),
                    numericInput("g_0", "g_0:", 1e-5),
                    numericInput("o2_0", "o2_0:", 1.2e-6),
                    numericInput("ox1_0", "ox1_0:", 1e-9),
                    numericInput("ox2_0", "ox2_0:", 1e-9),
                    numericInput("red1_0", "red1_0:", 0),
                    numericInput("red2_0", "red2_0:", 0)

                ),

# Show a table summarizing the values entered
                mainPanel(
                    plotOutput("System1"),
                    plotOutput("System2")
                )
            )
        ))
