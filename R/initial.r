N <- 10;
N_l <- 11;

layers<- 2;
grid_size <- N*layers + 1;
R_0 <- 0.1;
R<- 1;
dx1<- R_0/N;
dx2<- (R-R_0)/N;

deltat<- cumsum(c(0, rep(dx1, N), rep(dx2, N)));
dt<- (min(dx1, dx2))^2/(2);
dx_l <- c(rep(dx1, N-1), 0.05,rep(dx2, N-1));

E1RED_0 <- E2RED_0 <- 0;
P_0<- 0;
O2_0 <- 1.2*1e-6;
G_0 <- 1e-5;
E1OX_0 = E2OX_0 <- 1e-9;
KM1 <- 8*1e-5;
KM2 <- 2*1e-5;
#VMAX1 <- 4*1e-5;
#VMAX2 <- 1*1e-5;
KCAT1 <- 100;
KCAT2 <- 100;
K1 <- 4*1e-3;
K2 <- 2*1e-3;

DG_r <- 2.1*1e-6;
DP_r <- 2.1*1e-6;
DO2_r <- 0.66*1e-5;

DG_g <- 6.3*1e-6;
DP_g <- 6.3*1e-6;
DO2_g <- 2*1e-5;
TIME <- 10

params<-c(K1, K2, KM1, KM2, KCAT1, KCAT2, dt, N,
          P_0, G_0, O2_0, E1OX_0, E2OX_0, E1RED_0, E2RED_0,
          1, DG_r, DO2_r, DP_r, R_0,
          0, DG_g, DO2_g, DP_g, (R-R_0),
          TIME);
