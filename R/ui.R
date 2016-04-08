library(shiny)

# Define UI for slider demo application
shinyUI(fluidPage(

#  Application title
            titlePanel("Biser model"),

# Sidebar with sliders that demonstrate various available
# options
            sidebarLayout(
                sidebarPanel(
# Simple integer interval
                    actionButton("go", "Calculate"),
                    
                    numericInput("g_0", "g_0(10E-3):", 1),
                    numericInput("o2_0", "o2_0(10E-4):", 2.5),
                    
                    numericInput("km_1", "Km_1(10E-3):", 6.8),
                    numericInput("vmax_1", "Vmax_1(10E-5):", 4),
                    helpText("\n----------------------------------\n"),
                    

                    numericInput("D_gm", "laeyr1 D_{g, m}(10E-6):", 2.2),
                    numericInput("D_pm", "laeyr1 D_{p, m}(10E-6):", 2.2),
                    numericInput("D_o2m", "laeyr1 D_{o2, m}(10E-5):", 0.8),

                    numericInput("D_gd", "laeyr2 D_{g, d}(10E-6):", 6.6),
                    numericInput("D_pd", "laeyr2 D_{p, d}(10E-6):", 6.6),
                    numericInput("D_o2d", "laeyr2 D_{o2, d}(10E-5):", 2.4),

                    helpText("\n----------------------------------\n"),
                    
                    
                    sliderInput("alpha", label = h5("Proportion coef. alpha:"),
                                min = 0, max = 1, value = 0.50),
                    numericInput("T", "Time to simulate (s):", 100),
                    numericInput("d_m", "Micro-reactor interval length(cm):", 0.10),
                    numericInput("d_d", "Diffusion interval length(cm):", 0.02),
                    numericInput("d_b", "Baundary interval length(cm):", 0.03),
                    sliderInput("n", label = h5("Mesh on each interval(on each 3 layers!):"),
                                min = 10, max = 100, value = 30)
                    #numericInput("n", "Mesh on each interval:", 30)
                    

                ),

# Show a table summarizing the values entered
                mainPanel(
                    textOutput("text1"),
                    plotOutput("System2"),
                    plotOutput("System1"),
                    textOutput("text2")
                )
            )
        ))
