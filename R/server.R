library(shiny)
library(Rcpp)
dyn.load("../calculator.so");

result <- c();
params<-c(1,2);
result<- .Call("calculate", params);
dt<- 0.001;

#time<- seq(0, 201, 0.01)

# Define server logic for slider examples
shinyServer(function(input, output) {
    output$Plot <- renderPlot({
        randomVals()
    })

    randomVals <- eventReactive(input$go, {
        result<- .Call("calculate", params);
        
        plot(result$P);
    })
})