library(shiny)
library(Rcpp)

#dyn.unload("../calculator.so")
dyn.load("../calculator.so");
source("initial.r");
result<- .Call("calculate", params);

# Define server logic for slider examples
shinyServer(function(input, output) {
    
    
    getDt <- function(){
        dx <- min(input$d_m, input$d_d, input$d_b);
        dt <- dx^2/(2);
        print(paste("dt ", dt));
        return(dt);
    }
    
     output$text1 <- renderText({
         paste( "Estimated average time for calculations " , round(input$T/getDt()/300000*2, 2) , " (s)!" );
     	           
     })

    output$System2 <- renderPlot({
        randomVals2()

    })

    output$System1 <- renderPlot({
        randomVals1()
    })

    randomVals2 <- eventReactive(input$go, {
        dx1<- input$d_m/input$n;
        dx2<- input$d_d/input$n;
        dx3<- input$d_b/input$n;
        
        deltat<- cumsum(c(0, rep(dx1, input$n), rep(dx2, input$n), rep(dx3, input$n)));

        params<- c(
            input$km_1 * 1e-3, input$km_1 * 1e-3,
            input$vmax_1 * 1e-5, input$vmax_1 * 1e-5,
            getDt(), input$n,
            0., input$g_0 * 1e-3, input$o2_0 * 1e-4,
            input$alpha,
            1, input$D_gm * 1e-6, input$D_pm * 1e-6, input$D_o2m * 1e-5, input$d_m,
            0, input$D_gd * 1e-6, input$D_pd * 1e-6, input$D_o2d * 1e-5, input$d_d,
            0, input$D_gd * 1e-6, input$D_pd * 1e-6, input$D_o2d * 1e-5, input$d_b,
            input$T
        );
        
        result<- .Call("calculate", params);

        plot(deltat, result$G, type="l", xlab = "x(cm)", ylab="raw curves", col="blue");
        
        lines(deltat, result$P, col="red");
        lines(deltat, result$O2, col="green");
        
        
        legend("bottomright",legend=c("G", "P", "O_2"),
               text.col=c("blue","red", "green"), col=c("blue","red", "green"))
        
    })

    randomVals1 <- eventReactive(input$go, {

        dx1<- input$d_m/input$n;
        dx2<- input$d_d/input$n;
        dx3<- input$d_b/input$n;
        
        deltat<- cumsum(c(0, rep(dx1, input$n), rep(dx2, input$n), rep(dx3, input$n)));
     
        params<- c(
            input$km_1 * 1e-3, input$km_1 * 1e-3,
            input$vmax_1 * 1e-5, input$vmax_1 * 1e-5,
            getDt(), input$n,
            0., input$g_0 * 1e-3, input$o2_0 * 1e-4,
            input$alpha,
            1, input$D_gm * 1e-6, input$D_pm * 1e-6, input$D_o2m * 1e-5, input$d_m,
            0, input$D_gd * 1e-6, input$D_pd * 1e-6, input$D_o2d * 1e-5, input$d_d,
            0, input$D_gd * 1e-6, input$D_pd * 1e-6, input$D_o2d * 1e-5, input$d_b,
            input$T
        );
        
        result<- .Call("calculate", params);

        plot(deltat, result$G/input$g_0, ylim= c(0,1), type="l", xlab = "x(cm)", ylab="normalized (0,1)", col="blue");
        
        lines(deltat, result$P/max(result$P), col="red");
        lines(deltat, result$O2/input$o2_0, col="green");
        
        
        legend("bottomright",legend=c("G", "P", "O_2"),
        text.col=c("blue","red", "green"), col=c("blue","red", "green"))



    })
})
