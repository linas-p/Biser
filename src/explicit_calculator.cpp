/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <BiserLikeModel/utils.h>
#include <BiserLikeModel/constants.h>
#include <BiserLikeModel/explicit_calculator.h>
#include <BiserLikeModel/biosensor_information.h>

namespace BiserLikeModel {
void calculate_explicitly(struct bio_params *bio_info, void *ptr, \
                          void (*callback_crunched)(void *, int),
                          std::vector<double> * P,\
                          std::vector<double> * G,\
                          std::vector<double> * O2,\
                          std::vector<double> * Ox1,\
                          std::vector<double> * Ox2,\
                          std::vector<double> * Red1,\
                          std::vector<double> * Red2
                         ) {
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
    bool write_to_file = bio_info->write_to_file;

    // Sukuriami lokalūs kintamieji dėl optimizavimo
    double k1                    = bio_info->k1;
    double k2                    = bio_info->k2;
    double kcat1                    = bio_info->kcat1;
    double kcat2                    = bio_info->kcat2;
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

    double dr, dr0, dr1;
    //double v_max1                  = bio_info->vmax1;
    //double v_max2                  = bio_info->vmax2;

    // Sukuriamas rezultatų saugojimui skirtas failas
    if(write_to_file) {
        output_file = fopen(out_file_name, "w");
        fclose(output_file);
    }
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


	std::clock_t start = std::clock();
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

        dr = space_steps[layer];

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

                dr = space_steps[layer];
            } else {
                // Įskaičiuojama difuzijos įtaka
                current_g[a] = dt * Dg * ((last_g[a + 1] - 2 * last_g[a] + last_g[a - 1]) / (dr * dr) +
					  (last_g[a + 1] - last_g[a - 1]) / (a * 2 * dr * dr))+
                               last_g[a];
                current_pr[a] = dt * Dpr * ((last_pr[a + 1] - 2 * last_pr[a] + last_pr[a - 1]) / (dr * dr) +
					  (last_pr[a + 1] - last_pr[a - 1]) / (a * 2 * dr * dr))+
                                last_pr[a];

                current_o2[a] = dt * Do2 * ((last_o2[a + 1] - 2 * last_o2[a] + last_o2[a - 1]) / (dr * dr) +
					  (last_o2[a + 1] - last_o2[a - 1]) / (a * 2 * dr * dr))+
                                last_o2[a];


                // Jeigu sluoksnis yra fermentinis,
                // tuomet prisideda ir kinetikos dalis
                if (enz_layer) {
                    current_1ox[a] = dt * 2 * k1 * last_1red[a] + last_1ox[a];
                    current_1red[a] = -dt * 2 * k1 * last_1red[a] + \
                                      last_1red[a];
                    current_2ox[a] = -dt * 2 * k2 * last_1red[a] + last_2ox[a];
                    current_2red[a] = dt * 2 * k2 * last_1red[a] + last_2red[a];

                    kinetics_partg = dt * e1ox_0 * kcat1 * last_g[a] / \
                                     (last_g[a]+km1);
                    kinetics_parto2 = dt * e2ox_0 * kcat2 * last_o2[a] / \
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

            dr0 = space_steps[layer];

            Dg1 = bio_info->layers[layer + 1].Dg;
            Dpr1 = bio_info->layers[layer + 1].Dpr;
            Do21 = bio_info->layers[layer + 1].Do2;


            dr1 = space_steps[layer + 1];

            current_g[a] = (Dg1 * dr0 * current_g[a + 1] + \
                            Dg0 * dr1 * current_g[a - 1]) / \
                           (Dg1 * dr0 + Dg0 * dr1);
            current_pr[a] = (Dpr1 * dr0 * current_pr[a + 1] + \
                             Dpr0 * dr1 * current_pr[a - 1]) / \
                            (Dpr1 * dr0 + Dpr0 * dr1);
            current_o2[a] = (Do21 * dr0 * current_o2[a + 1] + \
                             Do20 * dr1 * current_o2[a - 1]) / \
                            (Do21 * dr0 + Do20 * dr1);
        }

        // Kraštinė substrato nepratekėjimo sąlyga
        current_g[0] = last_g[0] -  dt * Dg * 4 * (last_g[1] - last_g[0]) / (dr * dr) + k1 * last_g[0]/(km1 + last_g[0]);
        current_pr[0] = last_pr[0] +  dt * Dg * 4 * (last_pr[1] - last_pr[0]) / (dr * dr) + k1 * last_pr[0]/(km1 + last_pr[0]);
        current_o2[0] = last_o2[0] -  dt * Dg * 4 * (last_o2[1] - last_o2[0]) / (dr * dr) + k2 * last_o2[0]/(km2 + last_o2[0]);

        current_g[point_count - 1] = current_g[point_count - 2];
        current_pr[point_count - 1] = current_pr[point_count - 2];
        current_o2[point_count - 1] = current_o2[point_count - 2];


        // Skaičiuojamas srovės tankis
        i = bio_info->ne * F * bio_info->layers[0].Dpr * \
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
            printf("start %d %s \n", t, out_file_name);
            if(write_to_file) {
                output_file = fopen(out_file_name, "a");
                fprintf(output_file, "%e %e \n", i, execution_time);
                fclose(output_file);
            }
            printf("start %d \n", t);
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

	double duration = ( std::clock() - start ) / static_cast<double>(CLOCKS_PER_SEC);
	printf("operations per sec %d \n", t/duration);

    // Atspausdinamas paskutinis taškas
    if(write_to_file) {
        output_file = fopen(out_file_name, "a");
        fprintf(output_file, "%e %e\n", i, execution_time);
        fclose(output_file);
    }
    if (callback_crunched != NULL)
        callback_crunched(ptr, execution_time);

    concatenate_vals( last_g, G, point_count);
    concatenate_vals( last_pr, P, point_count);
    concatenate_vals( last_o2, O2, point_count);
    concatenate_vals( last_1ox, Ox1, point_count);
    concatenate_vals( last_2ox, Ox2, point_count);
    concatenate_vals( last_1red, Red1, point_count);
    concatenate_vals( last_2red, Red2, point_count);

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
