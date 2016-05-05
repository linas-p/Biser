source("initial.r");
source("utils.r");

last_l <- rep(0, grid_size);
last_l[(N_R0p):grid_size] <- L_0;

last_p <- rep(0, grid_size);

last_o2 <- rep(O2_0, grid_size);
last_o2[N_0:N_R0m] <- rho * O2_0;

current_l <- rep(-1, grid_size);
current_p <- rep(-1, grid_size);
current_o2 <- rep(-1, grid_size);

delta <- 1/((R^3-R_1^3)/(3*R_1^2));

tt<- c(); Ct_g <- c(); Ct_p <- c(); Ct_o2 <- c();

for(time in 0:99999) {

    kinetics_partl  <- MM(last_l[N_0], VMAX1 , KM1);
    kinetics_parto2 <- MM(last_o2[N_0], VMAX2 , KM2);
    kinetics <- MM2(kinetics_partl, kinetics_parto2);

    current_l[N_0] <- last_l[N_0]   + dt * (DL_m * LaplacePolar0(last_l[N_0], last_l[N_0 + 1], dx_m) - 2 * kinetics);
    current_p[N_0] <- last_p[N_0]   + dt * (DP_m * LaplacePolar0(last_p[N_0], last_p[N_0 + 1], dx_m) + 2 * kinetics);
    current_o2[N_0] <- last_o2[N_0] + dt * (DO2_m * LaplacePolar0(last_o2[N_0], last_o2[N_0 + 1], dx_m) - kinetics);

    
    for(k in 2:(N_R0m-1)) {
        kinetics_partl <- MM(last_l[k], VMAX1 , KM1);
		kinetics_parto2 <- MM(last_o2[k], VMAX2 , KM2);
		kinetics <- MM2(kinetics_partl, kinetics_parto2);
		
        current_l[k] <- last_l[k] + dt * (DL_m * LaplacePolar(last_l[k-1], last_l[k], last_l[k + 1], dx_m, points[k]) - 2 * kinetics);
        current_p[k] <- last_p[k] + dt * (DP_m * LaplacePolar(last_p[k-1], last_p[k], last_p[k + 1], dx_m, points[k]) + 2 * kinetics);
        current_o2[k] <- last_o2[k] + dt * (DO2_m * LaplacePolar(last_o2[k-1], last_o2[k], last_o2[k+1], dx_m, points[k]) - kinetics);
    }

    current_l[N_R0m] = rho * (DL_d * dx_m  * last_l[N_R0p+1] +
                       DL_m * dx_d * last_l[N_R0m-1]) /
                      (DL_d * dx_m + rho * DL_m * dx_d);
    current_p[N_R0m] = rho * (DP_d * dx_m  * last_p[N_R0p+1] +
                      DP_m * dx_d * last_p[N_R0m-1]) /
                      (DP_d * dx_m + rho * DP_m * dx_d);
    current_o2[N_R0m] = rho * (DO2_d * dx_m  * last_o2[N_R0p+1] +
                        DO2_m * dx_d * last_o2[N_R0m-1]) /
                       (DO2_d * dx_m + rho * DO2_m * dx_d);
    current_l[N_R0p] = (DL_d * dx_m * last_l[N_R0p+1] +
                       DL_m * dx_d * last_l[N_R0m-1]) /
                      (DL_d * dx_m + rho * DL_m * dx_d);
    current_p[N_R0p] = (DP_d * dx_m * last_p[N_R0p+1] +
                      DP_m * dx_d  * last_p[N_R0m-1]) /
                      (DP_d * dx_m + rho * DP_m * dx_d);
    current_o2[N_R0p] = (DO2_d * dx_m * last_o2[N_R0p+1] +
                        DO2_m * dx_d  * last_o2[N_R0m-1]) /
                       (DO2_d * dx_m + rho * DO2_m * dx_d);
    
    
    for(k in (N_R0p+1):(N_R1-1)) {
        current_l[k] <- last_l[k] + dt * DL_d * LaplacePolar(last_l[k-1], last_l[k], last_l[k+1], dx_d, points[k]);
        current_p[k] <- last_p[k] + dt * DP_d * LaplacePolar(last_p[k-1], last_p[k], last_p[k+1], dx_d, points[k]);
        current_o2[k] <- last_o2[k] + dt * DO2_d * LaplacePolar(last_o2[k-1], last_o2[k], last_o2[k+1], dx_d, points[k]);
    }

    #current_l[N_R1] = (DL_d * dx_d * last_l[N_R1+1] +
    #                     DL_m * dx_b * last_l[N_R1-1]) /
    #  (DL_d * dx_d + DL_m * dx_b);
    #current_p[N_R1] = (DP_d * dx_d * last_p[N_R1+1] +
    #                     DP_m * dx_b * last_p[N_R1-1]) /
    #  (DP_d * dx_d + DP_m * dx_b);
    #current_o2[N_R1] = (DO2_d * dx_d * last_o2[N_R1+1] +
    #                      DO2_m * dx_b * last_o2[N_R1-1]) /
    #  (DO2_d * dx_d + DO2_m * dx_b);



    for(k in N_R1:N_R) {
        current_l[k] <- last_l[k] - dt * delta * DL_d * (last_l[N_R1] - last_l[N_R1-1])/dx_d;
        current_p[k] <- last_p[k] - dt * delta * DP_d * (last_p[N_R1] - last_p[N_R1-1])/dx_d;
        current_o2[k] <- last_o2[k] - dt * delta * DO2_d * (last_o2[N_R1] - last_o2[N_R1-1])/dx_d;
    }


    last_l <- current_l;
    last_p <- current_p;
    last_o2 <- current_o2;

    if(!(time %% 10000)) {
        tt <- c(tt, time*dt);
        Ct_g <- c(Ct_g, 3 /(R^3 - R_0^3) * ( sum(diff(points[6:10]) * (current_l[6:9] + current_l[7:10])/2 * ((points[6:9] + points[7:10])/2)^2) + (R^3 - R_1^3)/3*current_l[10]));
        Ct_p <- c(Ct_p, 3 /(R^3 - R_0^3) * ( sum(diff(points[6:10]) * (current_p[6:9] + current_p[7:10])/2 * ((points[6:9] + points[7:10])/2)^2) + (R^3 - R_1^3)/3*current_p[10]));
        Ct_o2 <- c(Ct_o2, 3 /(R^3 - R_0^3) * ( sum(diff(points[6:10]) * (current_o2[6:9] + current_o2[7:10])/2 * ((points[6:9] + points[7:10])/2)^2) + (R^3 - R_1^3)/3*current_o2[10]));
        cat(".");
        plot(points, current_l, type='o', ylim = c(0, L_0));
        lines(points, current_o2, type='o');
    }
}

cbind(points, current_o2, result$O2)
cbind(points, current_l, result$L)
cbind(points, current_p, result$P)

plot(points, current_l)
lines(points, result$L)

plot(points, current_p)
lines(points, result$P)

plot(points, current_o2)
lines(points, result$O2)

plot(tt, Ct_g);
plot(tt, Ct_p);
plot(tt, Ct_o2);
