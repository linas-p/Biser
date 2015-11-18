library(shiny)
library(Rcpp)

dyn.load("../calculator.so");
source("initial.r");
result<- .Call("calculate", params);

# Define server logic for slider examples
shinyServer(function(input, output) {

    output$System2 <- renderPlot({
        randomVals2()

    })

    output$System1 <- renderPlot({
        randomVals1()
    })

    randomVals2 <- eventReactive(input$go, {
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

        result<- .Call("calculate", params);

        plot(deltat, result$O2/max(result$O2), type="l", xlab = "x(cm)", ylab="normalized (0,1)", col="blue");
        lines(deltat, result$Ox2/max(result$Ox2), col="green");
#lines(deltat, result$Red2, col="brown");
        legend("bottomright",legend=c("O2", "Ox2"),
        text.col=c("blue", "green", "brown"), col=c("blue", "green"))
    })

    randomVals1 <- eventReactive(input$go, {

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

        result<- .Call("calculate", params);

        plot(deltat, result$G/max(result$G), type="l", xlab = "x(cm)", ylab="normalized (0,1)", col="blue");
        lines(deltat, result$P/max(result$P), col="red");
        lines(deltat, result$Ox1/max(result$Ox1), col="green");
        lines(deltat, result$Red1/max(result$Red1), col="brown");
        legend("bottomright",legend=c("G","P", "Ox1", "Red1"),
        text.col=c("blue","red", "green", "brown"), col=c("blue","red", "green", "brown"))



    })
})
