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
                          std::vector<double> * points, \
                          std::vector<double> * P,\
                          std::vector<double> * L,\
                          std::vector<double> * O2, \
                          std::vector<double> * tim, \
                          std::vector<double> * Ct_l,\
                          std::vector<double> * Ct_p, \
                          std::vector<double> * Ct_o2
                         ) {

    int a;

    // Srovės tankis
    // double s_g_t, s_pr_t, s_o2_t;
    bool dimensionless = bio_info->dimensionless;
    bool return_all = false;

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
    double mm_l = 0., mm_o2 = 0., mm = 0., mm3 = 0.;

    // Rezultatų saugojimui skirtas failas
    FILE *output_file;
    bool write_to_file = bio_info->write_to_file;

    // Sukuriami lokalūs kintamieji dėl optimizavimo
    double km1                   = bio_info->km1;
    double km2                   = bio_info->km2;
    double dt                    = bio_info->dt;
    int n                        = bio_info->n;
    enum resp_method resp_t_meth = bio_info->resp_t_meth;
    double min_t                 = bio_info->min_t;
    double resp_t                = bio_info->resp_t;
    char *out_file_name          = bio_info->out_file_name;
    double o2_0                  = bio_info->o2_0;
    double l_0                   = bio_info->l_0;
    double rho                   = bio_info->rho;
    //double rho_inverse = 1/rho;
    int layer_count              = bio_info->layer_count;
    //double oxi = static_cast<double>(bio_info->oxigen);


    double Dl, Dl0, Dl1;
    double Dpr, Dpr0, Dpr1;
    double Do2, Do20, Do21;

    double dr, dr0, dr1;
    double v_max1                  = bio_info->vmax1;
    double v_max2                  = bio_info->vmax2;
    double beta = v_max1 / v_max2;
    double sigma1sq = 1., sigma2sq = 1.;
    int16_t N_0 = 0, N_R0m = n, N_R0p = n + 1, N_R1 = 2 * n + 1, N_R = 3 * n + 1;

    // Sukuriamas rezultatų saugojimui skirtas failas
    if(write_to_file) {
        output_file = fopen(out_file_name, "w");
        fclose(output_file);
    }
    // Apskaičiuojamas tinklo taškų skaičius per visus biojutiklio sluoksnius
    point_count = layer_count * n + 2;
    //printf("start ini %d \n", point_count);
    //fflush(stdout);

    // Medžiagų koncentracijų masyvams išskiriama atmintis
    double *current_l = new double[point_count],
    *last_l = new double[point_count];
    double *current_pr = new double[point_count],
    *last_pr = new double[point_count];
    double *current_o2 = new double[point_count],
    *last_o2 = new double[point_count];

    //Debug mode
    FillArray(current_l, -1, N_0, N_R);
    FillArray(current_pr, -1, N_0, N_R);
    FillArray(current_o2, -1, N_0, N_R);

    double delta_inverse = 0.;

    if (dimensionless) {
        sigma1sq = pow(bio_info->layers[MICROREACTOR].d, 2) * v_max1 / (bio_info->layers[MICROREACTOR].Dl * km1);
        sigma2sq = pow(bio_info->layers[MICROREACTOR].d, 2) * v_max2 / (bio_info->layers[MICROREACTOR].Dl * km2);
        dt = dt * bio_info->layers[MICROREACTOR].Dl / pow(bio_info->layers[MICROREACTOR].d, 2);
        delta_inverse = bio_info->layers[MICROREACTOR].d/((pow(bio_info->layers[MICROREACTOR].d + bio_info->layers[DIFFUSION].d + bio_info->layers[BAUDARY].d, 3) - pow(bio_info->layers[DIFFUSION].d + bio_info->layers[MICROREACTOR].d, 3))/(3 * pow(bio_info->layers[DIFFUSION].d + bio_info->layers[MICROREACTOR].d, 2)));
		for (a = layer_count - 1; a >= 0 ; a--) {
		    bio_info->layers[a].d = bio_info->layers[a].d / bio_info->layers[MICROREACTOR].d;
			bio_info->layers[a].Dpr = bio_info->layers[a].Dpr / bio_info->layers[MICROREACTOR].Dl;
			bio_info->layers[a].Do2 = bio_info->layers[a].Do2 / bio_info->layers[MICROREACTOR].Dl;
			bio_info->layers[a].Dl = bio_info->layers[a].Dl / bio_info->layers[MICROREACTOR].Dl;
		}
        o2_0 = o2_0 / km2;
        l_0 = l_0 / km1;

    } else {
        delta_inverse = 1/((pow(bio_info->layers[MICROREACTOR].d + bio_info->layers[DIFFUSION].d + bio_info->layers[BAUDARY].d, 3) - pow(bio_info->layers[DIFFUSION].d + bio_info->layers[MICROREACTOR].d, 3))/(3 * pow(bio_info->layers[DIFFUSION].d + bio_info->layers[MICROREACTOR].d, 2)));

    }

	// Priskiriamos pradinės ir kai kurios kraštinės sąlygos
	FillArray(last_pr, 0, N_0, N_R);
	FillArray(last_o2, o2_0, N_0, N_R);
	FillArray(last_o2, rho*o2_0, N_0, N_R0m);
	FillArray(last_l, 0, N_0, N_R);
	FillArray(last_l, l_0, N_R0p, N_R);
	//PrintArray(last_l, point_count);
	//PrintArray(last_o2, point_count);

    // Gražiname norimą rezultatą
    concatenate_vals(last_l, L, point_count);
    concatenate_vals(last_pr, P, point_count);
    concatenate_vals(last_o2, O2, point_count);

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
        //printfi"<-%d, %f\n", j, space_points[j]);
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
    concatenate_vals(space_points, points, point_count);



    std::clock_t start = std::clock();
	//double dt = std::pow(std::min(space_steps[MICROREACTOR], std::min(space_steps[DIFFUSION], space_steps[BAUDARY])), 2)/2;
	double c_l = averageConcentration(last_l, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R),
	c_p = averageConcentration(last_pr, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R),
	c_o2 = averageConcentration(last_o2, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R);
	//printf("simulated: %f \n", execution_time);
	tim->push_back(execution_time);
	Ct_l->push_back(c_l);
	Ct_p->push_back(c_p);
	Ct_o2->push_back(c_o2);

    //printf("start dt %.10f, stop on %f\n", dt, resp_t);
	const int t_total = (double)resp_t/dt;
	const int rate_div = std::max(100, t_total/DIVISION_RATE);
    //printf("total iterations T %d, rate %d\n", t_total, rate_div);

    printf("dt %f\n", dt);
    //printf("start delta %f\n", delta_inverse);
    printf("dimensionless %d\n", dimensionless);
    //printf("%f %f %f %f %f | %f %f %f %f\n", l_0, o2_0, sigma1sq, sigma2sq, beta, v_max1, v_max2, km1, km2);

    do {
        // Iteruojama per biojutiklio sluoksnius,
        // skaičiuojamos medžiagų koncentracijos

        // Surenkami pirmojo sluoksnio parametrai
        Dl  = bio_info->layers[MICROREACTOR].Dl;
        Dpr = bio_info->layers[MICROREACTOR].Dpr;
        Do2 = bio_info->layers[MICROREACTOR].Do2;
        dr = space_steps[MICROREACTOR];

        // Skaičiuojame MM taške 0
        mm_l  = MM(last_l[N_0], v_max1, km1);
        mm_o2 = MM(last_o2[N_0], v_max2, km2);
        mm3 = MM3(last_l[N_0], last_o2[N_0], beta);
        mm    = MM2(mm_l, mm_o2) * (1 - dimensionless) + dimensionless * mm3;

        // Kraštinė substrato nepratekėjimo sąlyga centre r = 0.
        current_l[N_0]  = last_l[N_0]  +  dt * (Dl  * LaplaceSperical0(last_l[N_0],  last_l[N_0 + 1],  dr) - 2 * sigma1sq * mm);
        current_pr[N_0] = last_pr[N_0] +  dt * (Dpr * LaplaceSperical0(last_pr[N_0], last_pr[N_0 + 1], dr) + 2 * sigma1sq * mm);
        current_o2[N_0] = last_o2[N_0] +  dt * (Do2 * LaplaceSperical0(last_o2[N_0], last_o2[N_0 + 1], dr) - sigma2sq * mm);

        //
        // Skaičiuojame sluoksnyje 0 < r < R_0
        for (a = N_0 + 1; a < N_R0m; a++) {
            mm_l  = MM(last_l[a], v_max1, km1);
            mm_o2 = MM(last_o2[a], v_max2, km2);
            mm3 = MM3(last_l[a], last_o2[a], beta);
            mm    = MM2(mm_l, mm_o2) * (1 - dimensionless) + dimensionless * mm3;
            // Įskaičiuojama difuzijos įtaka
            current_l[a]  = last_l[a]  + dt * (Dl  * LaplaceSperical(last_l[a - 1],  last_l[a],  last_l[a + 1],  dr, space_points[a]) - 2 * sigma1sq * mm);
            current_pr[a] = last_pr[a] + dt * (Dpr * LaplaceSperical(last_pr[a - 1], last_pr[a], last_pr[a + 1], dr, space_points[a]) + 2 * sigma1sq * mm);
            current_o2[a] = last_o2[a] + dt * (Do2 * LaplaceSperical(last_o2[a - 1], last_o2[a], last_o2[a + 1], dr, space_points[a]) - sigma2sq * mm);
		}

            // Sluoksnių sandūroms pritaikomos derinimo sąlygos taške R_0
            // TODO add normal with rho
            {

                Dl0 = bio_info->layers[MICROREACTOR].Dl;
                Dpr0 = bio_info->layers[MICROREACTOR].Dpr;
                Do20 = bio_info->layers[MICROREACTOR].Do2;
                dr0 = space_steps[MICROREACTOR];

                Dl1 = bio_info->layers[DIFFUSION].Dl;
                Dpr1 = bio_info->layers[DIFFUSION].Dpr;
                Do21 = bio_info->layers[DIFFUSION].Do2;
                dr1 = space_steps[DIFFUSION];

                current_l[N_R0m] = rho * (Dl1 * dr0 * last_l[N_R0p + 1] + \
                                          Dl0 * dr1 * last_l[N_R0m - 1]) / \
                                   (Dl1 * dr0 + rho * Dl0 * dr1);
                /*current_pr[N_R0m] = rho_inverse * (Dpr1 * dr0 * last_pr[N_R0p + 1] + \
                                 Dpr0 * dr1 * last_pr[N_R0m - 1]) / \
                                (Dpr1 * dr0 + rho_inverse * Dpr0 * dr1);*/
                current_pr[N_R0m] = rho * (Dpr1 * dr0 * last_pr[N_R0p + 1] + \
                                           Dpr0 * dr1 * last_pr[N_R0m - 1]) / \
                                    (Dpr1 * dr0 + rho * Dpr0 * dr1);
                current_o2[N_R0m] = rho * (Do21 * dr0 * last_o2[N_R0p + 1] + \
                                           Do20 * dr1 * last_o2[N_R0m - 1]) / \
                                    (Do21 * dr0 + rho * Do20 * dr1);

                current_l[N_R0p] = (Dl1 * dr0 * last_l[N_R0p + 1] + \
                                    Dl0 * dr1 * last_l[N_R0m - 1]) / \
                                   (Dl1 * dr0 + rho * Dl0 * dr1);
                /*current_pr[N_R0p] = (Dpr1 * dr0 * last_pr[N_R0p + 1] + \
                                 Dpr0 * dr1 * last_pr[N_R0m - 1]) / \
                                (Dpr1 * dr0 + rho_inverse * Dpr0 * dr1);*/
                current_pr[N_R0p] = (Dpr1 * dr0 * last_pr[N_R0p + 1] + \
                                     Dpr0 * dr1 * last_pr[N_R0m - 1]) / \
                                    (Dpr1 * dr0 + rho * Dpr0 * dr1);
                current_o2[N_R0p] = (Do21 * dr0 * last_o2[N_R0p + 1] + \
                                     Do20 * dr1 * last_o2[N_R0m - 1]) / \
                                    (Do21 * dr0 + rho * Do20 * dr1);
            }


            // Surenkami antrojo sluoksnio parametrai
            Dl  = bio_info->layers[DIFFUSION].Dl;
            Dpr = bio_info->layers[DIFFUSION].Dpr;
            Do2 = bio_info->layers[DIFFUSION].Do2;
            dr = space_steps[DIFFUSION];

            // Skaičiuojame sluoksnyje R_0 < r < R_1
            for (a = N_R0p + 1; a < N_R1; a++) {
                // Įskaičiuojama difuzijos įtaka
                current_l[a]  = last_l[a]  + dt * Dl  * LaplaceSperical(last_l[a - 1],  last_l[a],  last_l[a + 1],  dr, space_points[a]);
                current_pr[a] = last_pr[a] + dt * Dpr * LaplaceSperical(last_pr[a - 1], last_pr[a], last_pr[a + 1], dr, space_points[a]);
                current_o2[a] = last_o2[a] + dt * Do2 * LaplaceSperical(last_o2[a - 1], last_o2[a], last_o2[a + 1], dr, space_points[a]);

            }

            // Sluoksnių sandūroms pritaikomos derinimo sąlygos taške R_1
            /*{

                Dl0 = bio_info->layers[DIFFUSION].Dl;
                Dpr0 = bio_info->layers[DIFFUSION].Dpr;
                Do20 = bio_info->layers[DIFFUSION].Do2;
                dr0 = space_steps[DIFFUSION];

                Dl1 = bio_info->layers[BAUDARY].Dl;
                Dpr1 = bio_info->layers[BAUDARY].Dpr;
                Do21 = bio_info->layers[BAUDARY].Do2;
                dr1 = space_steps[BAUDARY];

                current_l[N_R1] = (Dl1 * dr0 * last_l[N_R1 + 1] + \
                                   Dl0 * dr1 * last_l[N_R1 - 1]) / \
                                  (Dl1 * dr0 + Dl0 * dr1);
                current_pr[N_R1] = (Dpr1 * dr0 * last_pr[N_R1 + 1] + \
                                    Dpr0 * dr1 * last_pr[N_R1 - 1]) / \
                                   (Dpr1 * dr0 + Dpr0 * dr1);
                current_o2[N_R1] = (Do21 * dr0 * last_o2[N_R1 + 1] + \
                                    Do20 * dr1 * last_o2[N_R1 - 1]) / \
                                   (Do21 * dr0 + Do20 * dr1);
            }*/

            // Surenkami trečiojo sluoksnio parametrai
            Dl  = bio_info->layers[BAUDARY].Dl;
            Dpr = bio_info->layers[BAUDARY].Dpr;
            Do2 = bio_info->layers[BAUDARY].Do2;
            dr = space_steps[BAUDARY];

            // Skaičiuojame sluoksnyje R_1 <= r <= R
            for (a = N_R1; a < N_R + 1; a++) {
                // Įskaičiuojama difuzijos įtaka
                current_l[a]  = last_l[a]  -  dt * Dl  * delta_inverse * (last_l[N_R1]  - last_l[N_R1-1])/dr;
                current_pr[a] = last_pr[a] -  dt * Dpr * delta_inverse * (last_pr[N_R1] - last_pr[N_R1-1])/dr;
                current_o2[a] = last_o2[a] -  dt * Do2 * delta_inverse * (last_o2[N_R1] - last_o2[N_R1-1])/dr;
            }



            // Apskaičiuojamas laikas
            t++;
            execution_time = t * dt;


            // Spausdinami rezultatai
            if ((t % rate_div) == 0) {
                //printf("start %d %s \n", t, out_file_name);
        		//double c_l = 0., c_p = 0., c_o2 = 0.;
				c_l = averageConcentration(current_l, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R);
				c_p = averageConcentration(current_pr, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R);
				c_o2 = averageConcentration(current_o2, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R);
                /*for(a = N_R0p; a < N_R1; a++) {
                    c_l += space_steps[DIFFUSION] * (current_l[a + 1] - current_l[a]) *
                           pow((space_points[a + 1] - space_points[a]), 2);
                    c_p += space_steps[DIFFUSION] * (current_pr[a + 1] - current_pr[a]) *
                           pow((space_points[a + 1] - space_points[a]), 2);
                    c_o2 += space_steps[DIFFUSION] * (current_o2[a + 1] - current_o2[a]) *
                            pow((space_points[a + 1] - space_points[a]), 2);
                }
                c_l += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_l[N_R1];
                c_p += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_pr[N_R1];
                c_o2 += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_o2[N_R1];

                c_l *= 3/(pow(N_R, 3) - pow(N_R0p, 3));
                c_p *= 3/(pow(N_R, 3) - pow(N_R0p, 3));
                c_o2 *= 3/(pow(N_R, 3) - pow(N_R0p, 3));*/

                if(write_to_file) {
                    output_file = fopen(out_file_name, "a");
                    fprintf(output_file, "%e, %e, %e, %e\n", execution_time, c_l, c_p, c_o2);
                    fclose(output_file);
                }
                //printf("simulated: %f \n", execution_time);
                if (return_all){
                    tim->push_back(execution_time);
                    Ct_l->push_back(c_l);
                    Ct_p->push_back(c_p);
                    Ct_o2->push_back(c_o2);
                }
				// Gražiname norimą rezultatą
				concatenate_vals(last_l, L, point_count);
				concatenate_vals(last_pr, P, point_count);
				concatenate_vals(last_o2, O2, point_count);

                //if (callback_crunched != NULL)
                //    callback_crunched(ptr, t);
            }

            // Masyvai sukeičiami vietomis
            SwapArrays(&current_l, &last_l);
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
            case HALF_TIME:
                if (c_o2 <= (o2_0 / 2)) {
                    response_time_reached = 1;
                    break;
                }
                break;
            }
        }
        while (!response_time_reached);

        double duration = (std::clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
        //printf("operations per sec %d \n", t/duration);
        //printf("total operation: %d, simulated time: %f \n", t, execution_time);
        //fflush(stdout);

		c_l = averageConcentration(current_l, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R);
		c_p = averageConcentration(current_pr, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R);
		c_o2 = averageConcentration(current_o2, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R);

        /*double c_l = 0., c_p = 0., c_o2 = 0.;*/
        /*for(a = N_R0p; a < N_R1; a++) {
            c_g += space_steps[DIFFUSION] * (current_l[a + 1] - current_l[a]) *
                   pow((space_points[a + 1] - space_points[a]), 2);
            c_p += space_steps[DIFFUSION] * (current_pr[a + 1] - current_pr[a]) *
                   pow((space_points[a + 1] - space_points[a]), 2);
            c_o2 += space_steps[DIFFUSION] * (current_o2[a + 1] - current_o2[a]) *
                    pow((space_points[a + 1] - space_points[a]), 2);
        }
        c_g += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_l[N_R];
        c_p += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_pr[N_R];
        c_o2 += (pow(N_R, 3) - pow(N_R1, 3))/3 * current_o2[N_R];

        c_g *= 3/(pow(N_R, 3) - pow(N_R0p, 3));
        c_p *= 3/(pow(N_R, 3) - pow(N_R0p, 3));
        c_o2 *= 3/(pow(N_R, 3) - pow(N_R0p, 3));*/

        // Atspausdinamas paskutinis taškas
        if(write_to_file) {
            output_file = fopen(out_file_name, "a");
            fprintf(output_file, "%e, %e, %e, %e\n", execution_time, c_l, c_p, c_o2);
            fclose(output_file);
        }
        tim->push_back(execution_time);
        Ct_l->push_back(c_l);
        Ct_p->push_back(c_p);
        Ct_o2->push_back(c_o2);

        if (callback_crunched != NULL)
            callback_crunched(ptr, t);

        // Gražiname norimą rezultatą
        concatenate_vals(last_l, L, point_count);
        concatenate_vals(last_pr, P, point_count);
        concatenate_vals(last_o2, O2, point_count);

        // Atlaisvinama atmintis
        free(current_l);
        free(last_l);
        free(current_pr);
        free(last_pr);
        free(current_o2);
        free(last_o2);

        free(space_steps);
        free(space_points);
    }



void two_layer_model(struct bio_params *bio_info, void *ptr, \
                     void (*callback_crunched)(void *, int),
                     std::vector<double> * points, \
                     std::vector<double> * P,\
                     std::vector<double> * L,\
                     std::vector<double> * tim, \
                     std::vector<double> * Ct_l,\
                     std::vector<double> * Ct_p,\
                     std::vector<double> * Ot_p, \
                     std::vector<double> * characteristics
                    ) {

    int a;

    // Srovės tankis
    // double s_g_t, s_pr_t, s_o2_t;
    bool return_all = false;

    // Žingsnių pagal erdvę masyvas
    double *space_steps;
    double *space_points;

    // Tinklo taškų skaičius per visus biojutiklio sluoksnius
    int point_count;

    // Iteracija pagal laiką
    int t = 0;

    // Simuliacijos laikas sekundėmis
    double execution_time = 0.;

    // Kintamasis nurodo ar jau pasiektas atsako laikas
    int response_time_reached;

    // Kinetikos dedamoji
    double mm = 0.;

    // Rezultatų saugojimui skirtas failas
    FILE *output_file;
    bool write_to_file = bio_info->write_to_file;

    // Sukuriami lokalūs kintamieji dėl optimizavimo
    double km1                   = bio_info->km1;
    double dt                    = bio_info->dt;
    int n                        = bio_info->n;
    enum resp_method resp_t_meth = bio_info->resp_t_meth;
    double min_t                 = bio_info->min_t;
    double resp_t                = bio_info->resp_t;
    char *out_file_name          = bio_info->out_file_name;
    double l_0                   = bio_info->l_0;
    double rho                   = bio_info->rho;
    //double rho_inverse = 1/rho;
    int layer_count              = bio_info->layer_count;
    //double oxi = static_cast<double>(bio_info->oxigen);
    double op = 0.;

    double Dl, Dl0, Dl1;
    double Dpr, Dpr0, Dpr1;


    double dr, dr0, dr1;
    double v_max1                  = bio_info->vmax1;

    int16_t N_0 = 0, N_R0m = n, N_R0p = n + 1, N_R1 = 2 * n + 1;

    // Sukuriamas rezultatų saugojimui skirtas failas
    if(write_to_file) {
        output_file = fopen(out_file_name, "w");
        fclose(output_file);
    }
    // Apskaičiuojamas tinklo taškų skaičius per visus biojutiklio sluoksnius
    point_count = layer_count * n + 2;
    //printf("start ini %d \n", point_count);
    //fflush(stdout);

    // Medžiagų koncentracijų masyvams išskiriama atmintis
    double *current_l = new double[point_count],
    *last_l = new double[point_count];
    double *current_pr = new double[point_count],
    *last_pr = new double[point_count];


    //Debug mode
    //FillArray(current_l, -1, N_0, N_R);
    //FillArray(current_pr, -1, N_0, N_R);
    double mse = 1.;

    double delta_inverse = 0.;

    delta_inverse = 1/((pow(bio_info->layers[MICROREACTOR].d + bio_info->layers[DIFFUSION].d + bio_info->layers[BAUDARY].d, 3) - pow(bio_info->layers[DIFFUSION].d + bio_info->layers[MICROREACTOR].d, 3))/(3 * pow(bio_info->layers[DIFFUSION].d + bio_info->layers[MICROREACTOR].d, 2)));


    // Priskiriamos pradinės ir kai kurios kraštinės sąlygos
    /*FillArray(last_pr, 0, N_0, N_R1);
    FillArray(last_l, 0, N_0, N_R1);
    FillArray(last_l, l_0, N_R1, N_R1);

    FillArray(current_pr, 0, N_0, N_R1);
    FillArray(current_l, 0, N_0, N_R1);
    FillArray(current_l, l_0, N_R1, N_R1);*/
    FillArray(last_pr, 0, N_0, N_R1);

    FillArray(last_l, 0, N_0, N_R0m);
    FillArray(last_l, 0, N_R0p, N_R1);


    // Gražiname norimą rezultatą
    concatenate_vals(last_l, L, point_count);
    concatenate_vals(last_pr, P, point_count);


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
        //printfi"<-%d, %f\n", j, space_points[j]);
    }
    space_points[N_R0p] = space_points[N_R0m];
    for (a = 1; a < layer_count; a++) {
        int j;
        for(j = a*n+2; j < a*n+2+n; j++) {
            space_points[j] = space_points[j-1] + space_steps[a];
            //printf("<-%d, %f\n", j, space_points[j]);
        }
    }



    //double dt = std::pow(std::min(space_steps[MICROREACTOR], std::min(space_steps[DIFFUSION], space_steps[BAUDARY])), 2)/2;
    double c_l = averageConcentration(last_l, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R1),
           c_p = averageConcentration(last_pr, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R1);

    if(write_to_file) {
        output_file = fopen(out_file_name, "a");
        fprintf(output_file, "%e, %e, %e\n", execution_time, c_l, c_p);
        fclose(output_file);
    }

    //PrintArray(space_points, point_count);
    concatenate_vals(space_points, points, point_count);
    //printf("simulated: %f \n", execution_time);
    tim->push_back(execution_time);
    Ct_l->push_back(c_l);
    Ct_p->push_back(c_p);
    Ot_p->push_back(bio_info->layers[DIFFUSION].Dpr * (current_pr[N_R0p + 1] - current_pr[N_R0p]) / space_steps[DIFFUSION]);

    //printf("start dt %.10f, stop on %f\n", dt, resp_t);
    const int t_total = (double)resp_t/dt;
    const int rate_div = std::max(100, t_total/DIVISION_RATE);
    printf("total iterations T %d, rate %d %d\n", t_total, rate_div, DIVISION_RATE);

    printf("dt %f\n", dt);
    //printf("start delta %f\n", delta_inverse);

    //std::clock_t start = std::clock();
    do {
        // Iteruojama per biojutiklio sluoksnius,
        // skaičiuojamos medžiagų koncentracijos

        // Surenkami pirmojo sluoksnio parametrai
        Dl  = bio_info->layers[MICROREACTOR].Dl;
        Dpr = bio_info->layers[MICROREACTOR].Dpr;

        dr = space_steps[MICROREACTOR];

        // Skaičiuojame MM taške 0
        mm  = MM(last_l[N_0], v_max1, km1);


        // Kraštinė substrato nepratekėjimo sąlyga centre r = 0.
        current_l[N_0]  = last_l[N_0]  +  dt * (Dl  * LaplaceSperical0(last_l[N_0],  last_l[N_0 + 1],  dr) -  mm);
        current_pr[N_0] = last_pr[N_0] +  dt * (Dpr * LaplaceSperical0(last_pr[N_0], last_pr[N_0 + 1], dr) +  mm);


        //
        // Skaičiuojame sluoksnyje 0 < r < R_0
        for (a = N_0 + 1; a < N_R0m; a++) {
            mm  = MM(last_l[a], v_max1, km1);
            // Įskaičiuojama difuzijos įtaka
            current_l[a]  = last_l[a]  + dt * (Dl  * LaplaceSperical(last_l[a - 1],  last_l[a],  last_l[a + 1],  dr, space_points[a]) -  mm);
            current_pr[a] = last_pr[a] + dt * (Dpr * LaplaceSperical(last_pr[a - 1], last_pr[a], last_pr[a + 1], dr, space_points[a]) +  mm);

        }




        // Surenkami antrojo sluoksnio parametrai
        Dl  = bio_info->layers[DIFFUSION].Dl;
        Dpr = bio_info->layers[DIFFUSION].Dpr;

        dr = space_steps[DIFFUSION];

        // Skaičiuojame sluoksnyje R_0 < r < R_1
        for (a = N_R0p + 1; a < N_R1; a++) {
            // Įskaičiuojama difuzijos įtaka
            current_l[a]  = last_l[a]  + dt * Dl  * LaplaceSperical(last_l[a - 1],  last_l[a],  last_l[a + 1],  dr, space_points[a]);
            current_pr[a] = last_pr[a] + dt * Dpr * LaplaceSperical(last_pr[a - 1], last_pr[a], last_pr[a + 1], dr, space_points[a]);


        }


        // Sluoksnių sandūroms pritaikomos derinimo sąlygos taške R_0
        {

            Dl0 = bio_info->layers[MICROREACTOR].Dl;
            Dpr0 = bio_info->layers[MICROREACTOR].Dpr;

            dr0 = space_steps[MICROREACTOR];

            Dl1 = bio_info->layers[DIFFUSION].Dl;
            Dpr1 = bio_info->layers[DIFFUSION].Dpr;

            dr1 = space_steps[DIFFUSION];

            current_l[N_R0m] = rho * (Dl1 * dr0 * current_l[N_R0p + 1] + \
                                      Dl0 * dr1 * current_l[N_R0m - 1]) / \
                               (Dl1 * dr0 + rho * Dl0 * dr1);
            /*current_pr[N_R0m] = rho_inverse * (Dpr1 * dr0 * last_pr[N_R0p + 1] + \
                             Dpr0 * dr1 * last_pr[N_R0m - 1]) / \
                            (Dpr1 * dr0 + rho_inverse * Dpr0 * dr1);*/
            current_pr[N_R0m] = rho * (Dpr1 * dr0 * current_pr[N_R0p + 1] + \
                                       Dpr0 * dr1 * current_pr[N_R0m - 1]) / \
                                (Dpr1 * dr0 + rho * Dpr0 * dr1);


            current_l[N_R0p] = (Dl1 * dr0 * current_l[N_R0p + 1] + \
                                Dl0 * dr1 * current_l[N_R0m - 1]) / \
                               (Dl1 * dr0 + rho * Dl0 * dr1);
            /*current_pr[N_R0p] = (Dpr1 * dr0 * last_pr[N_R0p + 1] + \
                             Dpr0 * dr1 * last_pr[N_R0m - 1]) / \
                            (Dpr1 * dr0 + rho_inverse * Dpr0 * dr1);*/
            current_pr[N_R0p] = (Dpr1 * dr0 * current_pr[N_R0p + 1] + \
                                 Dpr0 * dr1 * current_pr[N_R0m - 1]) / \
                                (Dpr1 * dr0 + rho * Dpr0 * dr1);
            op = Dpr1 * (current_pr[N_R0p + 1] - current_pr[N_R0p]) / dr1;
        }



        current_l[N_R1]  = l_0;
        current_pr[N_R1] = 0;

        // Skaičiuojame sluoksnyje R_1 <= r <= R
        /*for (a = N_R1; a < N_R + 1; a++) {
            current_l[a]  = l_0;
            current_pr[a] = current_pr[N_R1-1];

        }*/



        // Apskaičiuojamas laikas
        t++;
        execution_time = t * dt;


        // Spausdinami rezultatai
        if ((t % rate_div) == 0) {
            //printf("start %d %s \n", t, out_file_name);
            //double c_l = 0., c_p = 0., c_o2 = 0.;
            mse = calc_L2(current_pr, last_pr, point_count) + calc_L2(current_l, last_l, point_count);
            c_l = averageConcentration(current_l, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R1);
            c_p = averageConcentration(current_pr, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R1);
            //printf("error: %lf, %lf ? %d\n", mse, mse*1e28, mse < EPSILON);
            /*if ((t % (rate_div*10)) == 0) {
                double res = calc_L2(current_l, last_l, point_count);
                printf("error: %lf ? %d\n", res, res < EPSILON);
                fflush(stdout);
            }*/


            if(write_to_file) {
                output_file = fopen(out_file_name, "a");
                fprintf(output_file, "%e, %e, %e\n", execution_time, c_l, c_p);
                fclose(output_file);
            }
            //printf("simulated: %f \n", execution_time);
            if (return_all) {
                tim->push_back(execution_time);
                Ct_l->push_back(c_l);
                Ct_p->push_back(c_p);
                Ot_p->push_back(op);
                // Gražiname norimą rezultatą
                concatenate_vals(last_l, L, point_count);
                concatenate_vals(last_pr, P, point_count);


                double product_rate = -4. * M_PI * bio_info->layers[DIFFUSION].Dpr * \
                        space_points[N_R1] * space_points[N_R1] * \
                        (current_pr[N_R1] - current_pr[N_R1-1])/space_steps[DIFFUSION];
                double substrate_rate = 4/3 * M_PI * v_max1 * space_points[N_R1] * space_points[N_R1] * space_points[N_R1] * l_0 /(km1 + l_0);//-4. * M_PI * space_points[N_R1] *\
                        space_points[N_R1] * space_points[N_R1] * c_p / 3;
                        
                double efectiviness = 3 * averageRate(current_l, space_points, space_steps[DIFFUSION], N_R0p, N_R1, v_max1, km1)/ (space_points[N_R1]*space_points[N_R1]*space_points[N_R1] * MM(l_0, v_max1, km1));


                characteristics->push_back(product_rate);
                characteristics->push_back(substrate_rate);
                characteristics->push_back(product_rate/substrate_rate);
                characteristics->push_back(efectiviness);

            }


            //if (callback_crunched != NULL)
            //    callback_crunched(ptr, t);
        }


        // Masyvai sukeičiami vietomis
        SwapArrays(&current_l, &last_l);
        SwapArrays(&current_pr, &last_pr);


        // Nustatoma ar tęsti simuliaciją
        switch (resp_t_meth) {
        case MIN_TIME:
            if (execution_time < min_t) {
                response_time_reached = 0;
                break;
            }
            // Jeigu jau pasiekė minimalų laiką,


        case FIXED_TIME:
            response_time_reached = (execution_time >= resp_t);
            break;
            // tuomet tikrinama pagal DEFAULT_TIME sąlygas
        case DEFAULT_TIME:
           if (mse < EPSILON) {
                response_time_reached = 1;
                break;
            }

        }
    }
    while (!response_time_reached);

    /*double duration = (std::clock() - start);
    printf("Reference time: %f \n", duration / static_cast<double>(CLOCKS_PER_SEC));*/
    //printf("total operation: %d, simulated time: %f \n", t, execution_time);
    //fflush(stdout);

    c_l = averageConcentration(current_l, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R1);
    c_p = averageConcentration(current_pr, space_points, space_steps[DIFFUSION], N_R0p, N_R1, N_R1);


    // Atspausdinamas paskutinis taškas
    if(write_to_file) {
        output_file = fopen(out_file_name, "a");
        fprintf(output_file, "%e, %e, %e\n", execution_time, c_l, c_p);
        fclose(output_file);
    }


    tim->push_back(execution_time);
    Ct_l->push_back(c_l);
    Ct_p->push_back(c_p);
    Ot_p->push_back(op);
    // Gražiname norimą rezultatą
    concatenate_vals(last_l, L, point_count);
    concatenate_vals(last_pr, P, point_count);

    double product_rate = -4. * M_PI * bio_info->layers[DIFFUSION].Dpr * \
            space_points[N_R1] * space_points[N_R1] * \
            (current_pr[N_R1] - current_pr[N_R1-1])/space_steps[DIFFUSION];
    double substrate_rate = -4. * M_PI * bio_info->layers[DIFFUSION].Dpr * \
            space_points[N_R1] * space_points[N_R1] * \
            (current_l[N_R1] - current_l[N_R1-1])/space_steps[DIFFUSION];

    double efectiviness = 3 * averageRate(current_l, space_points, space_steps[DIFFUSION], N_0, N_R0m, v_max1, km1)/ (space_points[N_R0m]*space_points[N_R0m]*space_points[N_R0m] * MM(l_0, v_max1, km1));

    characteristics->push_back(product_rate);
    characteristics->push_back(substrate_rate);
    characteristics->push_back(product_rate/substrate_rate);
    characteristics->push_back(efectiviness);
    printf("error: %lf, %lf ? %d\n", mse, mse*1e28, mse < EPSILON);


    if (callback_crunched != NULL)
        callback_crunched(ptr, execution_time);




    // Atlaisvinama atmintis
    free(current_l);
    free(last_l);
    free(current_pr);
    free(last_pr);


    free(space_steps);
    free(space_points);
}
}  //  namespace BiserLikeModel
