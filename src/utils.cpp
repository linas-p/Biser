/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */

#include <BiserLikeModel/utils.h>

namespace BiserLikeModel {

void swap_arrays(double **array1, double **array2) {
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

}  // namespace BiserLikeModel
