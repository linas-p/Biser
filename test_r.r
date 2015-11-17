require(Rcpp);
#nm -D calculator_r.so | grep calculate
dyn.load("calculator.so");
params<-c(1,2);

result<- .Call("calculate", params);

