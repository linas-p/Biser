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
                    numericInput("T", "Time to simulate (s):", 25),
                    numericInput("n", "Mesh on each interval(on each 3 layers!):", 10),
                    numericInput("km_1", "Km_1(10E-3):", 9.6),
                    numericInput("km_2", "Km_2(10E-4):", 5),
                    numericInput("vmax_1", "Vmax_1(10E-4):", 1.9),
                    numericInput("vmax_2", "Vmax_2(10E-4):", 3.9),
                    
                    numericInput("rho", "Proportion coef.:", 0.56),
                    
                    numericInput("d_m", "Micro-reactor interval length(cm):", 0.025),
                    numericInput("d_d", "Diffusion interval length(cm):", 0.005),
                    numericInput("d_b", "Baundary interval length(cm):", 0.054),
                    

                    numericInput("D_lm", "laeyr1 D_{g, m}(10E-6):", 2.2),
                    numericInput("D_pm", "laeyr1 D_{p, m}(10E-6):", 2.2),
                    numericInput("D_o2m", "laeyr1 D_{o2, m}(10E-5):", 0.8),

                    numericInput("D_ld", "laeyr2 D_{g, d}(10E-6):", 6.7),
                    numericInput("D_pd", "laeyr2 D_{p, d}(10E-6):", 6.7),
                    numericInput("D_o2d", "laeyr2 D_{o2, d}(10E-5):", 2.4),

                    numericInput("l_0", "l_0(10E-3):", 2)

                ),

# Show a table summarizing the values entered
                mainPanel(
                	  textOutput("text1"),
                      plotOutput("System1"),
                      plotOutput("System2"),
                	  plotOutput("System4"),
                	  plotOutput("System3"),
                	  textOutput("text2")
                )
            )
        ))
