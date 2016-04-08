N <- 4;
N_0 <- 1;
N_R0m <- N + 1;
N_R0p <- N + 2;

N_R1 <- 2*N + 2;
N_R <- 3*N + 2;

layers<- 3;
grid_size <- N * layers + 2;

R_0 <- 0.10;
R_1 <- 0.12;
R<- 0.15;

dx_m <- R_0 / N;
dx_d <- (R_1-R_0) / N;
dx_b <- (R - R_1) / N;

points <- cumsum(c(0, rep(dx_m, N), 0, rep(dx_d, N), rep(dx_b, N)));

dt<- (min(dx_m, dx_d, dx_b))^2/(2);

P_0<- 0;
O2_0 <- 2.5*1e-4;
G_0 <- 1e-3;

KM1 <- 6.8*1e-3;
VMAX1 <- 4*1e-5;

DG_m <- 2.2*1e-6;
DP_m <- 2.2*1e-6;
DO2_m <- 0.8*1e-5;

DG_d <- DG_m*3;
DP_d <- DP_m*3;
DO2_d <- DO2_m*3;

alpha <- 0.5;

params<- c(
    KM1, KM1,
    VMAX1, VMAX1,
    dt, N,
    P_0, G_0, O2_0,
    alpha,
    1, DG_m, DO2_m, DP_m, 0.10,
    0, DG_d, DO2_d, DP_d, 0.02,
    0, DG_d, DO2_d, DP_d, 0.03,
    dt*100000, 0
);
