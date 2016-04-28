/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */

#include <BiserLikeModel/utils.h>

namespace BiserLikeModel {

double LaplacePolar(double *array, int k, double dr, double r) {
    double val = (array[k + 1] - 2 * array[k] + array[k - 1]) / (dr * dr)
                 + (1 / r) * (array[k + 1] - array[k - 1])/(dr);
    return val;

}

double LaplacePolar0(double * _array, double _dr) {
    return (2 * (_array[1] - _array[0]) / std::pow(_dr, 2));
}

double MM(double _val, double _vmax, double _km) {
    return (_vmax * _val) / (_km + _val);
}

double MM2(double _v1, double _v2) {
    double val = 0;
    if(Likely(_v1 || _v2)) {
        val = (_v1 * _v2) /(_v1 + 2 * _v2);
    } else {
        val = 0;
    }
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
