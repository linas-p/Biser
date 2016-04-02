source("initial.r");

last_g <- rep(0, grid_size);
last_g[(N_R0+1):grid_size] <- G_0;
last_p <- rep(0, grid_size);
last_o2 <- rep(O2_0, grid_size);
last_o2[N_0:N_R0] <- alpha * O2_0;



current_g <- rep(-1, grid_size);
current_p <- rep(-1, grid_size);
current_o2 <- rep(-1, grid_size);

delta <- 1/((R^3-R_1^3)/(3*R_1^2));


for(time in 0:1000000) {
    for(k in 2:(N_R0-1)) {
        kinetics_partg <- dt * VMAX1 * last_g[k] / (KM1 + last_g[k]);
        current_g[k] <- last_g[k] + dt * DG_m * ((last_g[k+1] - 2*last_g[k] +
                last_g[k-1])/(dx_m^2) +
                (2/k)*(last_g[k+1]-last_g[k])
                /(dx_m^2)) - kinetics_partg;
        current_p[k] <- last_p[k] + dt * DP_m * ((last_p[k+1] - 2*last_p[k] +
                last_p[k-1])/(dx_m^2) +
                (2/k)*(last_p[k+1]-last_p[k])
                /(dx_m^2)) + kinetics_partg;
        current_o2[k] <- last_o2[k] + dt * DO2_m * ((last_o2[k+1] - 2*last_o2[k] +
                                                      last_o2[k-1])/(dx_m^2) +
                                                     (2/k)*(last_o2[k+1]-last_o2[k])
                                                 /(dx_m^2)) - kinetics_partg;
    }
    
    kinetics_partg <- dt * VMAX1 * last_g[0] / (KM1 + last_g[0]);
    current_g[1] <- last_g[1] + dt * DG_m * 2 * (last_g[2] - last_g[1])/(dx_m^2)
    - kinetics_partg;
    current_p[1] <- last_p[1] + dt * DP_m * 2 * (last_p[2] - last_p[1])/(dx_m^2)
    + kinetics_partg;
    current_o2[1] <- last_o2[1] + dt * DO2_m * 2 * (last_o2[2] - last_o2[1])/(dx_m^2)
    - kinetics_partg;

    current_g[N_R0] = (DG_d * dx_m * last_g[N_R0+1] +
                       DG_m * dx_d * last_g[N_R0-1]) /
                      (DG_d * dx_m + DG_m * dx_d);
    current_p[N_R0] = (DP_d * dx_m * last_p[N_R0+1] +
                       DP_m * dx_d * last_p[N_R0-1]) /
                      (DP_d * dx_m + DP_m * dx_d);
    current_o2[N_R0] = (DO2_d * dx_m * last_o2[N_R0+1] +
                           DO2_m * dx_d * last_o2[N_R0-1]) /
        (DO2_d * dx_m + DO2_m * dx_d);

    for(k in (N_R0+1):(N_R1-1)) {
        current_g[k] <- last_g[k] + dt * DG_d * ((last_g[k+1] -
                2*last_g[k] + last_g[k-1]) /
                (dx_d^2) + (2/k)*(last_g[k+1]-last_g[k])/(dx_d^2));
        current_p[k] <- last_p[k] + dt * DP_d * ((last_p[k+1] -
                2*last_p[k] + last_p[k-1]) /
                (dx_d^2) + (2/k)*(last_p[k+1]-last_p[k])/(dx_d^2));
        current_o2[k] <- last_o2[k] + dt * DO2_d * ((last_o2[k+1] -
                                                      2*last_o2[k] + last_o2[k-1]) /
                                                     (dx_d^2) + (2/k)*(last_o2[k+1]-last_o2[k])/(dx_d^2));
    }

    current_g[N_R1] = (last_g[N_R1+1] + last_g[N_R1-1]) / 2;
    current_p[N_R1] = (last_p[N_R1+1] + last_p[N_R1-1]) / 2;
    current_o2[N_R1] = (last_o2[N_R1+1] + last_o2[N_R1-1]) / 2;
    


    for(k in (N_R1+1):N_R) {
        current_g[k] <- last_g[k] - dt * delta *(last_g[N_R1] -
                last_g[N_R1-1])/dx_d;
        current_p[k] <- last_p[k] - dt * delta *(last_p[N_R1] -
                last_p[N_R1-1])/dx_d;
        current_o2[k] <- last_o2[k] - dt * delta *(last_o2[N_R1] -
                                                     last_o2[N_R1-1])/dx_d;
    }


    last_g <- current_g;
    last_p <- current_p;
    last_o2 <- current_o2;
    
    if(!(time %% 100000)) {
        cat(".");
        plot(points, current_g, type='o', ylim = c(0, G_0));
        lines(points, current_o2, type='o');
    }
}

#cbind(points, current_g)
cbind(points, current_o2)

plot(points, current_o2, type='o');

