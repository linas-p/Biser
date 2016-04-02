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
                    numericInput("T", "Time to simulate (s):", 100),
                    numericInput("n", "Mesh on each interval(on each 3 layers!):", 30),
                    numericInput("km_1", "Km_1(10E-3):", 6.8),
                    numericInput("vmax_1", "Vmax_1(10E-5):", 4),
                    numericInput("alpha", "Proportion coef.:", 0.5),
                    
                    numericInput("d_m", "Micro-reactor interval length(cm):", 0.10),
                    numericInput("d_d", "Diffusion interval length(cm):", 0.02),
                    numericInput("d_b", "Baundary interval length(cm):", 0.03),
                    

                    numericInput("D_gm", "laeyr1 D_{g, m}(10E-6):", 2.2),
                    numericInput("D_pm", "laeyr1 D_{p, m}(10E-6):", 2.2),
                    numericInput("D_o2m", "laeyr1 D_{o2, m}(10E-5):", 0.8),

                    numericInput("D_gd", "laeyr2 D_{g, d}(10E-6):", 6.6),
                    numericInput("D_pd", "laeyr2 D_{p, d}(10E-6):", 6.6),
                    numericInput("D_o2d", "laeyr2 D_{o2, d}(10E-5):", 2.4),

                    numericInput("g_0", "g_0(10E-3):", 1),
                    numericInput("o2_0", "o2_0(10E-4):", 2.5)

                ),

# Show a table summarizing the values entered
                mainPanel(
                	  textOutput("text1"),
                    plotOutput("System1"),
                    plotOutput("System2"),
                	  textOutput("text2")
                )
            )
        ))
