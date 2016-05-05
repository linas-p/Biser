library(shiny)
library(Rcpp)


dyn.load("calculator.so");
source("initial.r");
result<- .Call("calculate", params);
dyn.unload("calculator.so")

# Define server logic for slider examples
shinyServer(function(input, output) {
  
    getDt <- function(){
        dx <- min(input$d_m, input$d_d, input$d_b);
        dt <- dx^2/(2);
        print(paste("dt ", dt));
        return(dt);
    }
    
     output$text1 <- renderText({
         paste( "Estimated average time for calculations " , round(input$T/getDt()/300000*2/(30/input$n), 2) , " (s)!\n" , 
                "Time step dt: ", getDt());
     	           
     })
     
     output$text2 <- renderText({
       print(paste( "--- "));
       
     })



    output$System1 <- renderPlot({
        randomVals1()
    })
    
    output$System2 <- renderPlot({
        randomVals2()
        
    })

    output$System4 <- renderPlot({
        randomVals4()
        
    })
    
    output$System3 <- renderPlot({
        randomVals3()
        
    })
    
    randomVals2 <- eventReactive(input$go, {

        load(file = "results")

        plot(result$dt, result$L, type="o", ylim = c(0, max(result$L)), xlab = "x(cm)", ylab="Concentations(M)", col="blue");
        
        lines(result$dt, result$P, type="o", col="red");
        lines(result$dt, result$O2, type="o", col="green");
        
        legend("bottomright",legend=c("L", "P", "O_2"),
        text.col=c("blue","red", "green"), col=c("blue","red", "green"));


    })
    
    randomVals3 <- eventReactive(input$go, {
        
        load(file = "results")
        
        plot(result$T, result$CL, ylim = c(0, max(result$CL)), type="o", xlab = "t(s)", ylab="C(t)", col="blue");
        
        lines(result$T, result$CP, type="o", col="red");
        lines(result$T, result$CO2, type="o", col="green");
        
        legend("bottomright",legend=c("C_L(t)", "C_P(t)", "C_O2(t)"),
               text.col=c("blue","red", "green"), col=c("blue","red", "green"));
        
        
    })
    
    randomVals4 <- eventReactive(input$go, {
        
        load(file = "results")
        
        plot(result$T, (result$CL - min(result$CL))/(max(result$CL) - min(result$CL)),  ylim=c(0, 1), type="o", xlab = "t(s)", ylab="normalized C(t)", col="blue");
        
        lines(result$T, (result$CP - min(result$CP))/(max(result$CP) - min(result$CP)), type="o", col="red");
        lines(result$T, (result$CO2 - min(result$CO2))/(max(result$CO2) - min(result$CO2)), type="o", col="green");
        
        legend("bottomright",legend=c("C_L(t)", "C_P(t)", "C_O2(t)"),
               text.col=c("blue","red", "green"), col=c("blue","red", "green"));
        
        
    })

    randomVals1 <- eventReactive(input$go, {
        dyn.load("calculator.so");
        dx1<- input$d_m/input$n;
        dx2<- input$d_d/input$n;
        dx3<- input$d_b/input$n;

        deltat<- cumsum(c(0, rep(dx1, input$n), 0, rep(dx2, input$n), rep(dx3, input$n)));

        params<- c(
            input$km_1 * 1e-3, input$km_2 * 1e-4,
            input$vmax_1 * 1e-4, input$vmax_2 * 1e-4,
            getDt(), input$n,
            0., input$l_0 * 1e-3, 2.5 * 1e-4,
            input$rho,
            1, input$D_lm * 1e-6, input$D_pm * 1e-6, input$D_o2m * 1e-5, input$d_m,
            0, input$D_ld * 1e-6, input$D_pd * 1e-6, input$D_o2d * 1e-5, input$d_d,
            0, input$D_ld * 1e-6, input$D_pd * 1e-6, input$D_o2d * 1e-5, input$d_b,
            input$T
        );

        result<- .Call("calculate", params);
        
        plot(deltat, result$L/(input$l_0 * 1e-3), type="o", ylim= c(0,1),  xlab = "x(cm)", ylab="normalized concentration", col="blue");
        
        lines(deltat, result$P/max(result$P), type="o", col="red");
        lines(deltat, result$O2/(2.5 * 1e-4), type="o", col="green");
        
        
        legend("bottomright",legend=c("L", "P", "O_2"),
        text.col=c("blue","red", "green"), col=c("blue","red", "green"))
        dyn.unload("calculator.so")
        result$dt <- deltat;
        save(result, file = "results");
        
        print(paste("< - ", ls()));
        
    })
})
