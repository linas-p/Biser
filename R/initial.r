N <- 4;
N_0 <- 1;
N_R0m <- N + 1;
N_R0p <- N + 2;

N_R1 <- 2*N + 2;
N_R <- 3*N + 2;

layers<- 3;
grid_size <- N * layers + 2;

R_0 <- 0.025;
R_1 <- 0.03;
R<- 0.084;

dx_m <- R_0 / N;
dx_d <- (R_1-R_0) / N;
dx_b <- (R - R_1) / N;

points <- cumsum(c(0, rep(dx_m, N), 0, rep(dx_d, N), rep(dx_b, N)));

dt<- (min(dx_m, dx_d, dx_b))^2/(2);

P_0<- 0;
O2_0 <- 2.5*1e-4;
L_0 <- 2e-3;

KM1 <- 9.6 * 1e-3;
KM2 <- 5 * 1e-4;
VMAX1 <- 1.9 * 1e-4;
VMAX2 <- 3.9 * 1e-4;

DL_m <- 2.2*1e-6;
DP_m <- 2.2*1e-6;
DO2_m <- 0.8*1e-5;

DL_d <- 6.7 * 1e-6;
DP_d <- 6.7 * 1e-6;
DO2_d <- 2.4 * 1e-5;

rho <- 0.56;

params<- c(
KM1, KM2,
VMAX1, VMAX2,
dt, N,
P_0, L_0, O2_0,
rho,
1, DL_m, DO2_m, DP_m, 0.025,
0, DL_d, DO2_d, DP_d, 0.005,
0, DL_d, DO2_d, DP_d, 0.034,
dt*10000
);



LaplacePolar <- function(vec, p, dr, r){
  #tmp <- ((vec[p+1] - 2*vec[p] + vec[p-1])/(dr^2) + (2/r)*(vec[p+1]-vec[p])/(dr));
  tmp <- ((vec[p+1] - 2*vec[p] + vec[p-1])/(dr^2) + (1/r)*(vec[p+1]-vec[p-1])/(dr));
  
  return(tmp);
}

LaplacePolar0 <- function(vec, dr){
  
  tmp <- 2*(vec[2] - vec[1])/(dr^2);
  
  return(tmp);
}

MM <- function(vec, p, vmax, km){
  return(vmax * vec[p]/(km+vec[p]));
}

MM2 <- function(v1, v2){
  return(v1*v2/(v1+2*v2));
}
