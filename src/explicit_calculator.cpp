/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <BiserLikeModel/utils.h>
#include <BiserLikeModel/constants.h>
#include <BiserLikeModel/explicit_calculator.h>
#include <BiserLikeModel/biosensor_information.h>

namespace BiserLikeModel {
void calculate_explicitly(struct bio_params *bio_info, void *ptr, \
                          void (*callback_crunched)(void *, int)) {
    int a;

    // Srovės tankis
    double i, last_i = 0;

    // Kintamasis rodo kaip pakito srovės tankis nuo praėjusios iteracijos
    double di;

    // Žingsnių pagal erdvę masyvas
    double *space_steps;

    // Tinklo taškų skaičius per visus biojutiklio sluoksnius
    int point_count;

    // Iteracija pagal laiką
    int t = 0;

    // Simuliacijos laikas sekundėmis
    double execution_time;

    // Kintamasis nurodo ar jau pasiektas atsako laikas
    int response_time_reached;

    // Kintamasis nurodo ties kuriuo sluoksniu esame
    int layer;

    // Kintamasis nurodo, kad esame ties sluoksnio kraštu
    // (ties sluoksnių sandūra)
    int is_boundary;

    // Kinetikos dedamoji
    double kinetics_partg;
    double kinetics_parto2;

    // Rezultatų saugojimui skirtas failas
    FILE *output_file;

    // Sukuriami lokalūs kintamieji dėl optimizavimo
    double k1                    = bio_info->k1;
    double k2                    = bio_info->k2;
    double km1                    = bio_info->km1;
    double km2                    = bio_info->km2;
    double dt                    = bio_info->dt;
    int n                        = bio_info->n;
    enum resp_method resp_t_meth = bio_info->resp_t_meth;
    double min_t                 = bio_info->min_t;
    double resp_t                = bio_info->resp_t;
    char *out_file_name          = bio_info->out_file_name;
    double e1ox_0                = bio_info->e1ox_0;
    double e2ox_0                = bio_info->e2ox_0;
    double o2_0                  = bio_info->o2_0;
    double g_0                  = bio_info->g_0;

    int layer_count              = bio_info->layer_count;
    int enz_layer;

    double Dg, Dg0, Dg1;
    double Dpr, Dpr0, Dpr1;
    double Do2, Do20, Do21;

    double dx, dx0, dx1;
    double v_max1                  = bio_info->vmax1;
    double v_max2                  = bio_info->vmax2;

    // Sukuriamas rezultatų saugojimui skirtas failas
    output_file = fopen(out_file_name, "w");
    fclose(output_file);

    // Apskaičiuojamas tinklo taškų skaičius per visus biojutiklio sluoksnius
    point_count = layer_count * n + 1;
    printf("start ini %d \n", point_count);
    fflush(stdout);

    // Medžiagų koncentracijų masyvams išskiriama atmintis
    double *current_g = new double[point_count],
    *last_g = new double[point_count];
    double *current_pr = new double[point_count],
    *last_pr = new double[point_count];
    double *current_o2 = new double[point_count],
    *last_o2 = new double[point_count];
    double *current_1red = new double[point_count],
    *last_1red = new double[point_count];
    double *current_2red = new double[point_count],
    *last_2red = new double[point_count];
    double *current_1ox = new double[point_count],
    *last_1ox = new double[point_count];
    double *current_2ox = new double[point_count],
    *last_2ox = new double[point_count];


    // Priskiriamos pradinės ir kai kurios kraštinės sąlygos
    fill_array(last_1ox, point_count, e1ox_0, 0);  //  not nessisary hole array
    fill_array(last_2ox, point_count, e2ox_0, 0);
    fill_array(last_1red, point_count, 0, 0);
    fill_array(last_2red, point_count, 0, 0);
    fill_array(last_pr, point_count, 0, 0);
    fill_array(last_o2, point_count, o2_0, 0);
    fill_array(last_g, point_count, 0, 0);
    fill_array(last_g, point_count, g_0, n);

    // Kiekvienam sluoksniui apskaičiuojami žingsniai pagal erdvę
    space_steps = new double[layer_count];
    for (a = 0; a < layer_count; a++)
        space_steps[a] = bio_info->layers[a].d / n;



    printf("start\n");
    do {
        // Iteruojama per biojutiklio sluoksnius,
        // skaičiuojamos medžiagų koncentracijos
        layer = 0;
        // Surenkami pirmojo sluoksnio parametrai
        enz_layer = bio_info->layers[layer].enz_layer;
        Dg = bio_info->layers[layer].Dg;
        Dpr = bio_info->layers[layer].Dpr;
        Do2 = bio_info->layers[layer].Do2;

        dx = space_steps[layer];

        for (a = 1; a < point_count - 1; a++) {
            // Nustatome ar tai nėra sluoksnių sandūra
            is_boundary = !(a % n);

            // Reikšmės sluoksnių sandūrose bus skaičiuojamos vėliau pagal
            // derinimo sąlygas
            if (is_boundary) {
                // Nustatome kuriame sluoksnyje esame
                layer++;
                // Surenkami kito sluoksnio parametrai
                enz_layer = bio_info->layers[layer].enz_layer;
                Dg = bio_info->layers[layer].Dg;
                Dpr = bio_info->layers[layer].Dpr;
                Do2 = bio_info->layers[layer].Do2;

                dx = space_steps[layer];
            } else {
                // Įskaičiuojama difuzijos įtaka
                current_g[a] = dt * Dg * \
                               (last_g[a + 1] - 2 * last_g[a] + \
                                last_g[a - 1]) / (dx * dx) + \
                               last_g[a];
                current_pr[a] = dt * Dpr * \
                                (last_pr[a + 1] - 2 * last_pr[a] + \
                                 last_pr[a - 1]) / (dx * dx) + \
                                last_pr[a];



                current_o2[a] = dt * Do2 * \
                                (last_o2[a + 1] - 2 * last_o2[a] \
                                 + last_o2[a - 1]) / (dx * dx) + \
                                last_o2[a];




                // Jeigu sluoksnis yra fermentinis,
                // tuomet prisideda ir kinetikos dalis
                if (enz_layer) {
                    current_1ox[a] = dt * 2 * k1 * last_1red[a] + last_1ox[a];
                    current_1red[a] = -dt * 2 * k1 * last_1red[a] + \
                                      last_1red[a];
                    current_2ox[a] = -dt * 4 * k2 * last_2ox[a] + last_2ox[a];
                    current_2red[a] = dt * 4 * k2 * last_2ox[a] + last_2red[a];

                    kinetics_partg = dt * v_max1 * last_g[a] / \
                                     (last_g[a]+km1);
                    kinetics_parto2 = dt * v_max2 * last_o2[a] / \
                                      (last_o2[a]+km2);

                    current_g[a] -=    kinetics_partg;
                    current_pr[a] +=   kinetics_partg;
                    current_1ox[a] -=  kinetics_partg;
                    current_1red[a] += kinetics_partg;

                    current_2ox[a] +=   kinetics_parto2;
                    current_2red[a] -=  kinetics_parto2;
                    current_o2[a] -=    kinetics_parto2;
                }
            }
        }

        // Sluoksnių sandūroms pritaikomos derinimo sąlygos
        for (layer = 0; layer < layer_count - 1; layer++) {
            // Apskaičiuojame kuriame taške yra layer
            // ir layer + 1 sluoksnių sandūra
            a = n * (layer + 1);
            Dg0 = bio_info->layers[layer].Dg;
            Dpr0 = bio_info->layers[layer].Dpr;
            Do20 = bio_info->layers[layer].Do2;

            dx0 = space_steps[layer];

            Dg1 = bio_info->layers[layer + 1].Dg;
            Dpr1 = bio_info->layers[layer + 1].Dpr;
            Do21 = bio_info->layers[layer + 1].Do2;


            dx1 = space_steps[layer + 1];

            current_g[a] = (Dg1 * dx0 * current_g[a + 1] + \
                            Dg0 * dx1 * current_g[a - 1]) / \
                           (Dg1 * dx0 + Dg0 * dx1);
            current_pr[a] = (Dpr1 * dx0 * current_pr[a + 1] + \
                             Dpr0 * dx1 * current_pr[a - 1]) / \
                            (Dpr1 * dx0 + Dpr0 * dx1);
            current_o2[a] = (Do21 * dx0 * current_o2[a + 1] + \
                             Do20 * dx1 * current_o2[a - 1]) / \
                            (Do21 * dx0 + Do20 * dx1);
        }

        // Kraštinė substrato nepratekėjimo sąlyga
        current_g[0] = current_g[1];
        current_pr[0] = current_pr[1];
        current_o2[0] = current_o2[1];

        current_g[point_count - 1] = current_g[point_count - 2];
        current_pr[point_count - 1] = current_pr[point_count - 2];
        current_o2[point_count - 1] = current_o2[point_count - 2];


        // Skaičiuojamas srovės tankis
        i = bio_info->ne * F * bio_info->layers[0].Dp * \
            (current_pr[1] - current_pr[0]) / space_steps[0];
        di = fabs(i - last_i);
        last_i = i;

        // TODO(linas): Combine oxygen and glucose.
        // condition_assing(current_1red, current_2ox,
        // point_count - 1, 2*k2/k1);



        // Masyvai sukeičiami vietomis
        swap_arrays(&current_g, &last_g);
        swap_arrays(&current_pr, &last_pr);
        swap_arrays(&current_1ox, &last_1ox);
        swap_arrays(&current_2ox, &last_2ox);
        swap_arrays(&current_1red, &last_1red);
        swap_arrays(&current_2red, &last_2red);
        swap_arrays(&current_o2, &last_o2);


        // Apskaičiuojamas laikas
        t++;
        execution_time = t * dt;


        // Spausdinami rezultatai
        if ((t % INTERVAL) == 0) {
            output_file = fopen(out_file_name, "a");
            fprintf(output_file, "%e %e \n", i, execution_time);
            fclose(output_file);
            if (callback_crunched != NULL)
                callback_crunched(ptr, execution_time);
        }

        // Nustatoma ar tęsti simuliaciją
        switch (resp_t_meth) {
        case MIN_TIME:
            if (execution_time < min_t) {
                response_time_reached = 0;
                break;
            }
            // Jeigu jau pasiekė minimalų laiką,
            // tuomet tikrinama pagal DEFAULT_TIME sąlygas
        case DEFAULT_TIME:
            if (i > 1e-30)
                response_time_reached = ((execution_time / i) * (di / dt) \
                                         <= EPSILON);
            else
                response_time_reached = 0;
            break;
        case FIXED_TIME:
            response_time_reached = (execution_time >= resp_t);
            break;
        }
    } while (!response_time_reached);

    // Atspausdinamas paskutinis taškas
    output_file = fopen(out_file_name, "a");
    fprintf(output_file, "%e %e\n", i, execution_time);
    fclose(output_file);
    if (callback_crunched != NULL)
        callback_crunched(ptr, execution_time);

    // Atlaisvinama atmintis
    free(current_g);
    free(last_g);
    free(current_pr);
    free(last_pr);
    free(current_1ox);
    free(last_1ox);
    free(current_2ox);
    free(last_2ox);
    free(current_1red);
    free(last_1red);
    free(current_2red);
    free(last_2red);
    free(current_o2);
    free(last_o2);

    free(space_steps);
}
}  //  namespace BiserLikeModel
