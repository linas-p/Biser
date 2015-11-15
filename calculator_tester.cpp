/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */

#include <stdio.h>
#include <stdlib.h>
#include <BiserLikeModel/biosensor_information.h>
#include <BiserLikeModel/explicit_calculator.h>

using namespace BiserLikeModel;

void callback_crunched(void *ptr, int time) {
    printf("%ds simulated\n", time);
}

int main() {
    struct bio_params *bio_info;

    bio_info = new bio_params;
    bio_info->explicit_scheme = 0;
    bio_info->subs_inh = 1;
    bio_info->prod_inh = 1;
    // [s^-1]
    bio_info->k2 = 2*1e-3;
    bio_info->k1 = 4*1e-3;

    // [mol/l] -> [mol/cm^3]
    bio_info->km = 0.01 * 1e-3;
    bio_info->km1 = 8 * 1e-5;
    bio_info->km2 = 2 * 1e-5;
    // [mol/l] -> [mol/cm^3]
    bio_info->ks = 0.001 * 1e-3;
    // [mol/l] -> [mol/cm^3]
    bio_info->kp = 0.001 * 1e-3;
    bio_info->vmax1 = 4 * 1e-5;
    bio_info->vmax2 = 1e-5;

    // [s]
    bio_info->dt = 1e-2;
    bio_info->n = 10;
    bio_info->resp_t_meth = MIN_TIME;
    // [s]
    bio_info->min_t = 100;
    // [s]
    bio_info->resp_t = 0;
    bio_info->out_file_name = "output.dat";
    bio_info->ne = 1;
    // [mol/l] -> [mol/cm^3] 4e-5;
    bio_info->s0 = 0.04 * 1e-3;
    // [mol/l] -> [mol/cm^3]
    bio_info->p0 = 0 * 1e-3;

    bio_info->pr_0 = 0 * 1e-3;
    bio_info->g_0 = 1e-5;
    bio_info->o2_0 = 1.2 * 1e-6;
    bio_info->e1ox_0 = 1e-9;
    bio_info->e2ox_0 = 1e-9;
    bio_info->e1red_0 = 0*1e-9;
    bio_info->e2red_0 = 0*1e-9;

    bio_info->layer_count = 2;
    bio_info->layers = new layer_params[bio_info->layer_count];

    // Užpildoma sluoksnių informacija
    // 0
    bio_info->layers[0].enz_layer = 1;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Ds = 100 * 1e-8;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Dp = 100 * 1e-8;

    bio_info->layers[0].Dg = 2.1 * 1e-6;
    bio_info->layers[0].Do2 = 0.666 * 1e-5;
    bio_info->layers[0].Dpr = 2.1 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[0].d = 0.1;
    // [mol/l] -> [mol/cm^3]
    bio_info->layers[0].e0 = 0.01 * 1e-3;

    // 1
    bio_info->layers[1].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Ds = 200 * 1e-8;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Dp = 200 * 1e-8;

    bio_info->layers[1].Dg = 6.3 * 1e-6;
    bio_info->layers[1].Do2 = 2 * 1e-5;
    bio_info->layers[1].Dpr = 6.3 * 1e-6;

    // [um] -> [cm]
    bio_info->layers[1].d = 0.9;
    // [mol/l] -> [mol/cm^3]
    bio_info->layers[1].e0 = 0 * 1e-3;

    calculate_explicitly(bio_info, NULL, &callback_crunched);

    free(bio_info->layers);
    free(bio_info);

    return 0;
}
