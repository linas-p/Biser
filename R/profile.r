#################
require(Rcpp);
library("parallel");

dyn.load("../calculator.so");
source("./initial.r");

gen_s0 <- function(index, optimi) { 
  Total <- 0.03;
  
  params <- c(
    100 * 1e-3, 0, # param: k_m,-
    optimi[index,1], 0,  # param: v_max,-
    1e-11, 240,# param: dt, N
    0, optimi[index, 2], 0,# param: -,s_0,-
    0.6,# param: phi
    1, 200 * 1e-6, 0, 200 * 1e-6, 250 * 1e-6,# param: is mm, d_m, -, d_m, h_0 = r_0
    0, 600 * 1e-6, 0, 600 * 1e-6, 60 * 1e-6,# param: is mm, d_d, -, d_d, h = r_1 - r_0
    0, 0, 0, 0, 0, # -,-,-,-
    5e-05 # time to simulate to 
  );
  result<- .Call("calculate", params);
  result$params <- optimi[index,]
  result$index <- index
  return(result);
}

rate <- 100 * 1e-3;# param: v_max
s0 <- 100 * 1e-3;# param: s_0
Gr <-  expand.grid(Vm = rate, S0 = s0);

res_grid <- mclapply(1:dim(Gr)[1], gen_s0, mc.cores = 8, optimi = Gr);

dyn.unload("../calculator.so")


c <- rainbow(8);
c <- c('black', c[1], c[4], c[6], c[8])

a <- res_grid[[1]]
nr0m <- 121
nr0p <- 122
nr1 <- 242
nr0m <- 241
nr0p <- 242
nr1 <- 482

setEPS()
postscript("profile.eps")

par(fig = c(0,1,0,1))
plot(a$points[1:nr0m], a$L[1:nr0m], type = 'l', xlim = c(0, max(a$points)), ylim = c(0, max(a$L)), 
     axes = FALSE, col = c[1], xlab = '', ylab = '', cex.lab=1.5, cex.main = 1.5, lwd = 4)
id <- c(1, as.integer(((1:10)/10) *  nr0m))
lines(a$points[nr0p:nr1], a$L[nr0p:nr1], col = c[1], lwd = 4)
lines(a$points[1:nr0m], a$L[1:nr0m+nr1*1], col = c[2], lty = 2, lwd =4)
id <- c(1, as.integer(((1:15)/15) *  nr0m))
lines(a$points[nr0p:nr1], a$L[nr0p:nr1+nr1*1], col = c[2], lty = 2, lwd = 4)
id <- c(1, as.integer(((1:5)/5) *  nr0m))
lines(a$points[1:nr0m], a$L[1:nr0m+nr1*200], col = c[3], lty = 3, lwd = 4)
id <- c(1, as.integer(((1:10)/10) *  nr0m))
lines(a$points[nr0p:nr1], a$L[nr0p:nr1+nr1*200], col = c[3], lty = 3, lwd = 4)
id <- c(1, as.integer(((1:5)/5) *  nr0m))
lines(a$points[1:nr0m], a$L[1:nr0m+nr1*500], col = c[4], lty = 4, lwd = 4)
id <- c(1, as.integer(((1:10)/10) *  nr0m))
id <- sample(nr0m, 3)
lines(a$points[nr0p:nr1], a$L[nr0p:nr1+nr1*500], col = c[4], lty = 4, lwd = 4)
lines(a$points[1:nr0m], a$L[1:nr0m+nr1*1001], col = c[5], lty = 5, lwd = 4)
id <- c(1, as.integer(((1:10)/10) *  nr0m))
lines(a$points[nr0p:nr1], a$L[nr0p:nr1+nr1*1001], col = c[5], lty = 5, lwd = 4)
arrows(275 / 1000000, 15/ 1000000, 255/ 1000000, 5/ 1000000, lwd = 2)
arrows(220 / 1000000, 25/ 1000, 240/ 1000000, 20/ 1000, lwd = 2)
arrows(150 / 1000000, 35/ 1000, 175/ 1000000, 30/ 1000, lwd = 2)
arrows(30 / 1000000, 35/ 1000, 55/ 1000000, 30/ 1000, lwd = 2)
arrows(30 / 1000000, 60/ 1000, 55/ 1000000, 55/ 1000, lwd = 2)
arrows(230 / 1000000, 95/ 1000, 250/ 1000000, 95/ 1000, lwd = 2)
arrows(275 / 1000000, 75/ 1000, 260/ 1000000, 80/ 1000, lwd = 2)
text(280 / 1000000, 20/ 1000000, "1", cex = 1.5)
text(215 / 1000000, 25/ 1000, "2", cex = 1.5)
text(145 / 1000000, 35/ 1000, "3", cex = 1.5)
text(25 / 1000000, 35/ 1000, "4",  cex = 1.5)
text(25 / 1000000, 65/ 1000, "5",  cex = 1.5)
text(225 / 1000000, 95/ 1000, "3", cex = 1.5)
text(280 / 1000000, 70/ 1000, "2", cex = 1.5)
pp <- c(seq(0, 0.00025, 0.00005), 0.00028, 0.00031);
axis(1, at=pp, labels=pp* 1000000, las=1, col.ticks="black", col.axis="black", col = "black", cex.axis=1.5)   
ss <- seq(0, 1e-01, 2e-2)
axis(2, at=ss, labels=ss* 1000, las=1, col.ticks="black", col.axis="black", col = "black", cex.axis=1.5)   
mtext(TeX("$r$, \ \ $\\mu m$"), 1, line=3, at=0.00015, cex =1.5, cex.axis=1.5)
mtext(TeX("$s$, \ \ $\\mu M$"), 2, line=2, at=5e-02, cex =1.5, cex.axis=1.5)
abline(v = 0.00025, col = 'gray', lty = 2, cex = 1.25, lwd = 2)
dev.off()



