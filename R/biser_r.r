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


for(time in 1:100000) {

    kinetics_partg <- MM(last_g, 1, VMAX1 , KM1);
    current_g[1] <- last_g[1] + dt * (DG_m * LaplacePolar0(last_g, dx_m) - kinetics_partg);
    current_p[1] <- last_p[1] + dt * (DP_m * LaplacePolar0(last_p, dx_m) + kinetics_partg);
    current_o2[1] <- last_o2[1] + dt * (DO2_m * LaplacePolar0(last_o2, dx_m) - kinetics_partg);

    
    for(k in 2:(N_R0-1)) {
        kinetics_partg <- MM(last_g, k, VMAX1 , KM1);
        current_g[k] <- last_g[k] + dt * (DG_m * LaplacePolar(last_g, k, dx_m, points[k]) - kinetics_partg);
        current_p[k] <- last_p[k] + dt * (DP_m * LaplacePolar(last_p, k, dx_m, points[k]) + kinetics_partg);
        current_o2[k] <- last_o2[k] + dt * (DO2_m * LaplacePolar(last_o2, k, dx_m, points[k]) - kinetics_partg);
    }



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
        current_g[k] <- last_g[k] + dt * DG_d * LaplacePolar(last_g, k, dx_d, points[k]);
        current_p[k] <- last_p[k] + dt * DP_d * LaplacePolar(last_p, k, dx_d, points[k]);
        current_o2[k] <- last_o2[k] + dt * DO2_d * LaplacePolar(last_o2, k, dx_d, points[k]);
    }

    current_g[N_R1] = (DG_d * dx_d * last_g[N_R1+1] +
                         DG_m * dx_b * last_g[N_R1-1]) /
      (DG_d * dx_d + DG_m * dx_b);
    current_p[N_R1] = (DP_d * dx_d * last_p[N_R1+1] +
                         DP_m * dx_b * last_p[N_R1-1]) /
      (DP_d * dx_d + DP_m * dx_b);
    current_o2[N_R1] = (DO2_d * dx_d * last_o2[N_R1+1] +
                          DO2_m * dx_b * last_o2[N_R1-1]) /
      (DO2_d * dx_d + DO2_m * dx_b);



    for(k in (N_R1+1):N_R) {
        current_g[k] <- last_g[k] - dt * delta * DG_d * (last_g[N_R1] - last_g[N_R1-1])/dx_d;
        current_p[k] <- last_p[k] - dt * delta * DP_d * (last_p[N_R1] - last_p[N_R1-1])/dx_d;
        current_o2[k] <- last_o2[k] - dt * delta * DO2_d * (last_o2[N_R1] - last_o2[N_R1-1])/dx_d;
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

cbind(points, current_o2, result$O2)
cbind(points, current_g, result$G)
cbind(points, current_p, result$P)

plot(points, current_g)
lines(points, result$G)

plot(points, current_p)
lines(points, result$P)

plot(points, current_o2)
lines(points, result$O2)

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