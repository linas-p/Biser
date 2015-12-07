source("initial.r");

last_g <- rep(0, grid_size);
last_g[N_l:grid_size] <- G_0;
last_p <- rep(0, grid_size);
last_o2 <- rep(O2_0, grid_size);
last_1ox <- rep(E1OX_0, grid_size);
last_2ox <- rep(E2OX_0, grid_size);
last_1red <- rep(0, grid_size);
last_2red <- rep(0, grid_size);

current_g <- rep(0, grid_size);
current_p <- rep(0, grid_size);
current_o2 <- rep(0, grid_size);
current_1ox <- rep(0, grid_size);
current_2ox <- rep(0, grid_size);
current_1red <- rep(0, grid_size);
current_2red <- rep(0, grid_size);



for(time in 0:100) {
#pin<- function()
#    {

    for(k in 2:(N_l-1)) {
        kinetics_partg <- dt*VMAX1*last_g[k]/(KM1+last_g[k]);
        kinetics_parto2 <- dt*VMAX2*last_o2[k]/(KM2+last_o2[k]);


        current_g[k] <- last_g[k] + got_dif(last_g, k, DG_r, dt, dx1) - kinetics_partg;
        current_p[k] <- last_p[k] + got_dif(last_p, k, DP_r, dt, dx1) + kinetics_partg;
        current_1ox[k] <- last_1ox[k] + dt*2*K1*last_1red[k]  - kinetics_partg;
        current_1red[k] <- last_1red[k] - dt*2*K1*last_1red[k]+ kinetics_partg;

        current_2ox[k] <- last_2ox[k] - dt*4*K2*last_2ox[k]   + kinetics_parto2;
        current_2red[k] <- last_2red[k] + dt*4*K2*last_2ox[k] - kinetics_parto2;
        current_o2[k] <- last_o2[k] + got_dif(last_o2, k, DO2_r, dt, dx1) - kinetics_parto2;
    }

    for(k in (N_l+1):(grid_size-1)) {
        current_g[k] <- last_g[k] + got_dif(last_g, k, DG_g, dt, dx2);
        current_p[k] <- last_p[k] + got_dif(last_p, k, DP_g, dt, dx2);
        current_o2[k] <- last_o2[k] + got_dif(last_o2, k, DO2_g, dt, dx2);
    }

    current_g[N_l] = got_mid(current_g[N_l+1], DG_r, current_g[N_l-1], DG_g);
    current_p[N_l] = got_mid(current_p[N_l+1], DP_r, current_p[N_l-1], DP_g);
    current_o2[N_l] = got_mid(current_o2[N_l+1], DO2_r, current_o2[N_l-1], DO2_g);


#current_g[grid_size] <- G_0;
    current_g[1] <- current_g[2];
    current_p[1] <- current_p[2];
    current_o2[1] <- current_o2[2];

    current_g[grid_size] <- current_g[grid_size-1];
    current_p[grid_size] <- current_p[grid_size-1];
    current_o2[grid_size] <- current_o2[grid_size-1];

# Where put this?
# current_1red <- 2*K2*current_2ox/K1;

    last_g <- current_g;
    last_p <- current_p;
    last_o2 <- current_o2;
    last_1red <- current_1red;
    last_1ox <- current_1ox;
    last_2red <- current_2red;
    last_2ox <- current_2ox;


}

plot(deltat, current_g, type='l');
lines(deltat, current_p, type='l');
lines(deltat, current_o2, type='l');

current_g
current_p
current_1ox
current_2ox
current_1red
current_2red
current_o2


output_i <- read.table("~/Documents/de/numerical/testing/cpp/simul/Biosensor-Calculator-Library/output.dat", quote="\"", comment.char="");
plot(output_i$V2, output_i$V1);
lines(output_e$V2, output_e$V1);


got_mid <- function(a_1, D1, a_2, D2) {
    a <- (D2 * dx1 * a_1 + D1 * dx2 * a_2) /  (D2 * dx1 + D1 * dx2);
    return(a);
}

got_dif <- function(a1, k, D, dt, dx) {
    a <- dt * D * (a1[k + 1] - 2 * a1[k] + a1[k - 1]) / (dx^2)
    return(a);
}
