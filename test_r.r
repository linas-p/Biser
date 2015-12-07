require(Rcpp);
#nm -D calculator_r.so | grep calculate
dyn.load("calculator.so");
source("initial.r");

result<- .Call("calculate", params);

