/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */

#include <BiserLikeModel/utils.h>

namespace BiserLikeModel {

double LaplacePolar(double *array, int k, double dr, double r) {
    double val = (array[k+1] - 2 * array[k] + array[k-1])/(dr*dr)
                 + (1/r)* (array[k+1] - array[k-1])/(dr);
    return val;

}

double LaplacePolar0(double *array, double dr) {
    double val = 2*(array[1] - array[0])/(dr*dr);
    return val;
}

double MM(double *array, int k, double vmax, double km) {
    double val = vmax*array[k]/(km+array[k]);
    return val;
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
