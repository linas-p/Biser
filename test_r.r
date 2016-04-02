require(Rcpp);
#nm -D calculator_r.so | grep calculate
dyn.load("calculator.so");
source("R/initial.r");

result<- .Call("calculate", params);
dyn.unload("calculator.so")

require(Rcpp);
dyn.load("../calculator.so");
source("initial.r");

result<- .Call("calculate", params);
dyn.unload("../calculator.so")
