library(shiny)
library(Rcpp)

dyn.load("../calculator.so");
source("initial.r");
#result<- .Call("calculate", params);
#print(params);

# Define server logic for slider examples
shinyServer(function(input, output) {
    output$Plot <- renderPlot({
        randomVals()
    })

    randomVals <- eventReactive(input$go, {

        dx1<- input$d*(input$decimal/100)/input$n;
        dx2<- input$d*((100-input$decimal)/100)/input$n;

        deltat<- cumsum(c(0, rep(dx1, input$n), rep(dx2, input$n)));
#dt<- (min(dx1, dx2))^2/(2);
        dt <- input$dt;
        dx_l <- c(rep(dx1, input$n-1), 0.05, rep(dx2, input$n-1));

        params<-c(input$k1, input$k2, input$km_1, input$km_2,
        input$vmax_1, input$vmax_2, dt, input$n,
        input$p_0, input$g_0, input$o2_0, input$ox1_0, input$ox2_0, input$red1_0, input$red2_0,
        1, input$D_gr, input$D_o2r, input$D_pr, input$d/100,
        0, input$D_gg, input$D_o2g, input$D_pg, (100-input$d)/100 );
        print(params);

        result<- .Call("calculate", params);

        plot(deltat, result$O2, xlab = "x(cm)", col="blue");
#lines(deltat, result$P, col="red");
        lines(deltat, result$G, col="green");


    })
})
