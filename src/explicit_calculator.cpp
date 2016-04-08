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
                          std::vector<double> * O2, \
                          std::vector<double> * tim, \
                          std::vector<double> * Ct_g,\
                          std::vector<double> * Ct_p, \
                          std::vector<double> * Ct_o2
                         ) {
    int a;

    // Srovės tankis
    // double s_g_t, s_pr_t, s_o2_t;


    // Žingsnių pagal erdvę masyvas
    double *space_steps;
    double *space_points;

    // Tinklo taškų skaičius per visus biojutiklio sluoksnius
    int point_count;

    // Iteracija pagal laiką
    int t = 0;

    // Simuliacijos laikas sekundėmis
    double execution_time;

    // Kintamasis nurodo ar jau pasiektas atsako laikas
    int response_time_reached;

    // Kinetikos dedamoji
    double mm_g = 0., mm_o2 = 0.;

    // Rezultatų saugojimui skirtas failas
    FILE *output_file;
    bool write_to_file = bio_info->write_to_file;

    // Sukuriami lokalūs kintamieji dėl optimizavimo
    double km1                    = bio_info->km1;
    double dt                    = bio_info->dt;
    int n                        = bio_info->n;
    enum resp_method resp_t_meth = bio_info->resp_t_meth;
    double min_t                 = bio_info->min_t;
    double resp_t                = bio_info->resp_t;
    char *out_file_name          = bio_info->out_file_name;
    double o2_0                  = bio_info->o2_0;
    double g_0                  = bio_info->g_0;
    double alpha                 = bio_info->alpha;
    double alpha_inverse = 1/alpha;
    int layer_count              = bio_info->layer_count;
    double oxi = static_cast<double>(bio_info->oxigen);

    double Dg, Dg0, Dg1;
    double Dpr, Dpr0, Dpr1;
    double Do2, Do20, Do21;

    double dr, dr0, dr1;
    double v_max1                  = bio_info->vmax1;
    //  double v_max2                  = bio_info->vmax2;
//	int16_t N_0 = 0, N_R0 = n, N_R1 = 2 * n, N_R = 3 * n;
    int16_t N_0 = 0, N_R0m = n, N_R0p = n + 1, N_R1 = 2 * n + 1, N_R = 3 * n + 1;

    // Sukuriamas rezultatų saugojimui skirtas failas
    if(write_to_file) {
        output_file = fopen(out_file_name, "w");
        fclose(output_file);
    }
    // Apskaičiuojamas tinklo taškų skaičius per visus biojutiklio sluoksnius
    point_count = layer_count * n + 2;
    printf("start ini %d \n", point_count);
    fflush(stdout);

    // Medžiagų koncentracijų masyvams išskiriama atmintis
    double *current_g = new double[point_count],
    *last_g = new double[point_count];
    double *current_pr = new double[point_count],
    *last_pr = new double[point_count];
    double *current_o2 = new double[point_count],
    *last_o2 = new double[point_count];

    //Debug mode
    FillArray(current_g, -1, N_0, N_R);
    FillArray(current_pr, -1, N_0, N_R);
    FillArray(current_o2, -1, N_0, N_R);

    // Priskiriamos pradinės ir kai kurios kraštinės sąlygos
    FillArray(last_pr, 0, N_0, N_R);
    FillArray(last_o2, o2_0, N_0, N_R);
    FillArray(last_o2, alpha*o2_0, N_0, N_R0m);
    FillArray(last_g, 0, N_0, N_R);
    FillArray(last_g, g_0, N_R0p, N_R);
    //PrintArray(last_g, point_count);
    //PrintArray(last_o2, point_count);

    // Kiekvienam sluoksniui apskaičiuojami žingsniai pagal erdvę
    space_steps = new double[layer_count];
    space_points = new double[point_count];
    space_points[N_0] = 0.;

    for (a = 0; a < layer_count; a++) {
        space_steps[a] = bio_info->layers[a].d / n;
    }


    int j;
    for(j = 1; j < N_R0p; j++) {
        space_points[j] = space_points[j-1] + space_steps[0];
        //printf("<-%d, %f\n", j, space_points[j]);
    }
    space_points[N_R0p] = space_points[N_R0m];
    for (a = 1; a < layer_count; a++) {
        int j;
        for(j = a*n+2; j < a*n+2+n; j++) {
            space_points[j] = space_points[j-1] + space_steps[a];
            //printf("<-%d, %f\n", j, space_points[j]);
        }
    }
    //PrintArray(space_points, point_count);

    double delta_inverse = 1/((pow(bio_info->layers[BAUDARY].d, 3) - pow(bio_info->layers[DIFFUSION].d, 3))/(3 * pow(bio_info->layers[DIFFUSION].d, 2)));

    std::clock_t start = std::clock();
    printf("start delta %f\n", delta_inverse);
    do {
        // Iteruojama per biojutiklio sluoksnius,
        // skaičiuojamos medžiagų koncentracijos

        // Surenkami pirmojo sluoksnio parametrai
        Dg  = bio_info->layers[MICROREACTOR].Dg;
        Dpr = bio_info->layers[MICROREACTOR].Dpr;
        Do2 = bio_info->layers[MICROREACTOR].Do2;
        dr = space_steps[MICROREACTOR];

        // Skaičiuojame MM taške 0
        mm_g = MM(last_g, 0, v_max1, km1);
        mm_o2 = MM(last_o2, 0, v_max1, km1);

        // Kraštinė substrato nepratekėjimo sąlyga centre r = 0.
        current_g[N_0]  = last_g[N_0]  +  dt * (Dg * LaplacePolar0(last_g, dr)   - 
((1 - oxi) * mm_g + oxi * mm_o2));
        current_pr[N_0] = last_pr[N_0] +  dt * (Dpr * LaplacePolar0(last_pr, dr) + 
((1 - oxi) * mm_g + oxi * mm_o2));
        current_o2[N_0] = last_o2[N_0] +  dt * (Do2 * LaplacePolar0(last_o2, dr) - 
((1 - oxi) * mm_g + oxi * mm_o2));

        //
        // Skaičiuojame sluoksnyje 0 < r < R_0
        for (a = N_0 + 1; a < N_R0m; a++) {
            mm_g = MM(last_g, a, v_max1, km1);
            mm_o2 = MM(last_o2, a, v_max1, km1);
            // Įskaičiuojama difuzijos įtaka
            current_g[a]  = last_g[a]  + dt * (Dg * LaplacePolar(last_g, a, dr, space_points[a])   -
((1 - oxi) * mm_g + oxi * mm_o2));
            current_pr[a] = last_pr[a] + dt * (Dpr * LaplacePolar(last_pr, a, dr, space_points[a]) + 
((1 - oxi) * mm_g + oxi * mm_o2));
            current_o2[a] = last_o2[a] + dt * (Do2 * LaplacePolar(last_o2, a, dr, space_points[a]) -
((1 - oxi) * mm_g + oxi * mm_o2));
        }

        // Sluoksnių sandūroms pritaikomos derinimo sąlygos taške R_0
        // TODO add normal with alpha
        {

            Dg0 = bio_info->layers[MICROREACTOR].Dg;
            Dpr0 = bio_info->layers[MICROREACTOR].Dpr;
            Do20 = bio_info->layers[MICROREACTOR].Do2;
            dr0 = space_steps[MICROREACTOR];

            Dg1 = bio_info->layers[DIFFUSION].Dg;
            Dpr1 = bio_info->layers[DIFFUSION].Dpr;
            Do21 = bio_info->layers[DIFFUSION].Do2;
            dr1 = space_steps[DIFFUSION];

            current_g[N_R0m] = alpha * (Dg1 * dr0 * last_g[N_R0p + 1] + \
                                        Dg0 * dr1 * last_g[N_R0m - 1]) / \
                               (Dg1 * dr0 + alpha * Dg0 * dr1);
            current_pr[N_R0m] = alpha_inverse * (Dpr1 * dr0 * last_pr[N_R0p + 1] + \
                                                 Dpr0 * dr1 * last_pr[N_R0m - 1]) / \
                                (Dpr1 * dr0 + alpha_inverse * Dpr0 * dr1);
            current_o2[N_R0m] = alpha * (Do21 * dr0 * last_o2[N_R0p + 1] + \
                                         Do20 * dr1 * last_o2[N_R0m - 1]) / \
                                (Do21 * dr0 + alpha * Do20 * dr1);

            current_g[N_R0p] = (Dg1 * dr0 * last_g[N_R0p + 1] + \
                                Dg0 * dr1 * last_g[N_R0m - 1]) / \
                               (Dg1 * dr0 + alpha * Dg0 * dr1);
            current_pr[N_R0p] = (Dpr1 * dr0 * last_pr[N_R0p + 1] + \
                                 Dpr0 * dr1 * last_pr[N_R0m - 1]) / \
                                (Dpr1 * dr0 + alpha_inverse * Dpr0 * dr1);
            current_o2[N_R0p] = (Do21 * dr0 * last_o2[N_R0p + 1] + \
                                 Do20 * dr1 * last_o2[N_R0m - 1]) / \
                                (Do21 * dr0 + alpha * Do20 * dr1);
        }


        // Surenkami antrojo sluoksnio parametrai
        Dg  = bio_info->layers[DIFFUSION].Dg;
        Dpr = bio_info->layers[DIFFUSION].Dpr;
        Do2 = bio_info->layers[DIFFUSION].Do2;
        dr = space_steps[DIFFUSION];

        // Skaičiuojame sluoksnyje R_0 < r < R_1
        for (a = N_R0p + 1; a < N_R1; a++) {
            // Įskaičiuojama difuzijos įtaka
            current_g[a]  = last_g[a]  + dt * Dg * LaplacePolar(last_g, a, dr, space_points[a]);
            current_pr[a] = last_pr[a] + dt * Dpr * LaplacePolar(last_pr, a, dr, space_points[a]);
            current_o2[a] = last_o2[a] + dt * Do2 * LaplacePolar(last_o2, a, dr, space_points[a]);

        }

        // Sluoksnių sandūroms pritaikomos derinimo sąlygos taške R_1
        {

            Dg0 = bio_info->layers[DIFFUSION].Dg;
            Dpr0 = bio_info->layers[DIFFUSION].Dpr;
            Do20 = bio_info->layers[DIFFUSION].Do2;
            dr0 = space_steps[DIFFUSION];

            Dg1 = bio_info->layers[BAUDARY].Dg;
            Dpr1 = bio_info->layers[BAUDARY].Dpr;
            Do21 = bio_info->layers[BAUDARY].Do2;
            dr1 = space_steps[BAUDARY];

            current_g[N_R1] = (Dg1 * dr0 * last_g[N_R1 + 1] + \
                               Dg0 * dr1 * last_g[N_R1 - 1]) / \
                              (Dg1 * dr0 + Dg0 * dr1);
            current_pr[N_R1] = (Dpr1 * dr0 * last_pr[N_R1 + 1] + \
                                Dpr0 * dr1 * last_pr[N_R1 - 1]) / \
                               (Dpr1 * dr0 + Dpr0 * dr1);
            current_o2[N_R1] = (Do21 * dr0 * last_o2[N_R1 + 1] + \
                                Do20 * dr1 * last_o2[N_R1 - 1]) / \
                               (Do21 * dr0 + Do20 * dr1);
        }

        // Surenkami trečiojo sluoksnio parametrai
        Dg  = bio_info->layers[BAUDARY].Dg;
        Dpr = bio_info->layers[BAUDARY].Dpr;
        Do2 = bio_info->layers[BAUDARY].Do2;
        dr = space_steps[BAUDARY];

        // Skaičiuojame sluoksnyje R_1 < r <= R
        for (a = N_R1 + 1; a < N_R + 1; a++) {
            // Įskaičiuojama difuzijos įtaka
            current_g[a]  = last_g[a]  -  dt * Dg  * delta_inverse * (last_g[N_R1]  - last_g[N_R1-1])/dr;
            current_pr[a] = last_pr[a] -  dt * Dpr * delta_inverse * (last_pr[N_R1] - last_pr[N_R1-1])/dr;
            current_o2[a] = last_o2[a] -  dt * Do2 * delta_inverse * (last_o2[N_R1] - last_o2[N_R1-1])/dr;
        }



        // Apskaičiuojamas laikas
        t++;
        execution_time = t * dt;


        // Spausdinami rezultatai
        if ((t % PRINT_RATE) == 0) {
            //printf("start %d %s \n", t, out_file_name);

            double c_g = 0., c_p = 0., c_o2 = 0.;
            for(a = N_R0p; a < N_R1; a++) {
                c_g += space_steps[DIFFUSION] * (current_g[a + 1] - current_g[a]) *
                       pow((space_points[a + 1] - space_points[a]), 2);
                c_p += space_steps[DIFFUSION] * (current_pr[a + 1] - current_pr[a]) *
                       pow((space_points[a + 1] - space_points[a]), 2);
                c_o2 += space_steps[DIFFUSION] * (current_o2[a + 1] - current_o2[a]) *
                        pow((space_points[a + 1] - space_points[a]), 2);
            }
            c_g += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_g[N_R];
            c_p += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_pr[N_R];
            c_o2 += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_o2[N_R];

            c_g *= 3/(pow(N_R, 3) - pow(N_R0p, 3));
            c_p *= 3/(pow(N_R, 3) - pow(N_R0p, 3));
            c_o2 *= 3/(pow(N_R, 3) - pow(N_R0p, 3));

            if(write_to_file) {
                output_file = fopen(out_file_name, "a");
                fprintf(output_file, "%e, %e, %e, %e\n", execution_time, c_g, c_p, c_o2);
                fclose(output_file);
            }
            //printf("simulated: %f \n", execution_time);
            tim->push_back(execution_time);
            Ct_g->push_back(c_g);
            Ct_p->push_back(c_p);
            Ct_o2->push_back(c_o2);
            //if (callback_crunched != NULL)
            //    callback_crunched(ptr, t);
        }

        // Masyvai sukeičiami vietomis
        SwapArrays(&current_g, &last_g);
        SwapArrays(&current_pr, &last_pr);
        SwapArrays(&current_o2, &last_o2);

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
            break;
        case FIXED_TIME:
            response_time_reached = (execution_time >= resp_t);
            break;
        }
    } while (!response_time_reached);

    double duration = (std::clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
    printf("operations per sec %d \n", t/duration);
    printf("total operation: %d, simulated time: %f \n", t, execution_time);
    fflush(stdout);

    double c_g = 0., c_p = 0., c_o2 = 0.;
    for(a = N_R0p; a < N_R1; a++) {
        c_g += space_steps[DIFFUSION] * (current_g[a + 1] - current_g[a]) *
               pow((space_points[a + 1] - space_points[a]), 2);
        c_p += space_steps[DIFFUSION] * (current_pr[a + 1] - current_pr[a]) *
               pow((space_points[a + 1] - space_points[a]), 2);
        c_o2 += space_steps[DIFFUSION] * (current_o2[a + 1] - current_o2[a]) *
                pow((space_points[a + 1] - space_points[a]), 2);
    }
    c_g += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_g[N_R];
    c_p += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_pr[N_R];
    c_o2 += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_o2[N_R];

    c_g *= 3/(pow(N_R, 3) - pow(N_R0p, 3));
    c_p *= 3/(pow(N_R, 3) - pow(N_R0p, 3));
    c_o2 *= 3/(pow(N_R, 3) - pow(N_R0p, 3));

    // Atspausdinamas paskutinis taškas
    if(write_to_file) {
        output_file = fopen(out_file_name, "a");
        fprintf(output_file, "%e, %e, %e, %e\n", execution_time, c_g, c_p, c_o2);
        fclose(output_file);
    }
    tim->push_back(execution_time);
    Ct_g->push_back(c_g);
    Ct_p->push_back(c_p);
    Ct_o2->push_back(c_o2);

    if (callback_crunched != NULL)
        callback_crunched(ptr, t);

    // Gražiname norimą rezultatą
    concatenate_vals(last_g, G, point_count);
    concatenate_vals(last_pr, P, point_count);
    concatenate_vals(last_o2, O2, point_count);

    // Atlaisvinama atmintis
    free(current_g);
    free(last_g);
    free(current_pr);
    free(last_pr);
    free(current_o2);
    free(last_o2);

    free(space_steps);
    free(space_points);
}
}  //  namespace BiserLikeModel
