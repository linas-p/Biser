/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */

#include "BiserLikeModel/utils.h"

namespace BiserLikeModel {

double Laplace(double valm, double valc, double valp, double dh) {
    double val = (valp - 2 * valc + valm) / (dh * dh);
    return val;
}


double LaplaceSperical(double valm, double valc, double valp, double dr, double r) {
    double val = (valp - 2 * valc + valm) / (dr * dr)
                 + (2. / r) * (valp - valm)/(dr);
    return val;
}

double Laplace0(double valc, double valp, double dh) {
    return (2 * (valp - valc) / std::pow(dh, 2));
}


double LaplaceSperical0(double valc, double valp, double dr) {
    return (2 * (valp - valc) / std::pow(dr, 2));
}

double MM(double val, double vmax, double km) {
    return (vmax * val) / (km + val);
}

double MM3(double l1, double o2, double beta){
    return (beta * l1 * o2) / (beta*l1 + (2 + beta) * l1 * o2 + 2 * o2);
}

double MM2(double v1, double v2) {
    double val = 0;
    if(Likely(v1 || v2)) {
        val = (v1 * v2) /(v1 + 2 * v2);
    } else {
        val = 0;
    }
    return val;
}

double averageConcentration(double *array, double *points, double delta, int r_0p, int r_1, int r) {
    double conc = 0.;
    int a;
	//printf("<- %f %f %f ", points[r_0p], points[r_1], points[r]);

    for(a = r_0p; a < r_1; a++) {
        conc += delta * ((array[a + 1] + array[a]) / 2) * pow((points[a + 1] + points[a]) / 2, 2);
    }

    conc += (pow(points[r], 3) - pow(points[r_1], 3))/3 * array[r_1];

	//printf("<- %f %f %f ", 3/(pow(points[r], 3) - pow(points[r_0p], 3)));
    conc *= 3/(pow(points[r], 3) - pow(points[r_0p], 3));

	return conc;
}


double averageRate(double *array, double *points, double delta, int r_0, int r_1, double vmax, double km) {
    double rate = 0.;
    int a;
	//printf("<- %f %f %f ", points[r_0p], points[r_1], points[r]);

    for(a = r_0; a < r_1; a++) {
        rate += delta * MM((array[a + 1] + array[a]) / 2, vmax, km) * pow((points[a + 1] + points[a]) / 2, 2);
    }



	return rate;
}

void SwapArrays(double **array1, double **array2) {
    double *temp;

    temp = *array1;
    *array1 = *array2;
    *array2 = temp;
}

void FillArray(double *array, double value, int from, int to) {
    for (int a = from; a < to + 1; a++) {
        array[a] = value;
    }
}

void condition_assing(double *array1, double *array2, int length, \
                      double value) {
    for (int a = 0; a < length; a++)
        array1[a] = array2[a]*value;
}

double calc_L2(double *array1, double *array2, int length) {
    double result = 0.;
    for (int a = 0; a < length; a++)
         result += (array2[a] - array1[a])*(array2[a] - array1[a]);
    result /= (double)length;
    return result;
}

void PrintArray(double *array, int length) {
    int a;

    printf("--- \n");
    for (a = 0; a < length; a++)
        printf(" %.10f", array[a]);

    printf("\n--- \n");
    fflush(stdout);
}

void concatenate_vals(double *x, std::vector<double> * out, int length) {
    std::vector<double> tmp(x, (x + length));
    out->insert(out->end(), tmp.begin(), tmp.end());
}

void FillArray(std::vector<double> * out, double value) {
    out->push_back(value);
}

}  // namespace BiserLikeModel
