/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */
#ifndef INCLUDE_BISERLIKEMODEL_UTILS_H_
#define INCLUDE_BISERLIKEMODEL_UTILS_H_

#include <stdio.h>
#include <vector>

namespace BiserLikeModel {
double LaplacePolar(double *array, int k, double dr, double r);
double LaplacePolar0(double *array, double dr);
double MM(double *array, int k, double vmax, double km);

void SwapArrays(double **array1, double **array2);
void condition_assing(double *array1, double *array2, int length, double value);

void FillArray(double *array, double value, int from, int to);
void PrintArray(double *array, int length);
void concatenate_vals(double *, std::vector<double> *, int );

}  //  namespace BiserLikeModel
#endif  //  INCLUDE_BISERLIKEMODEL_UTILS_H_
