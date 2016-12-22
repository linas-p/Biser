/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */
#ifndef INCLUDE_BISERLIKEMODEL_UTILS_H_
#define INCLUDE_BISERLIKEMODEL_UTILS_H_
#define Unlikely(x) __builtin_expect((bool)(x), 0)
#define Likely(x)   __builtin_expect((bool)(x), 1)

#include <stdio.h>
#include <vector>
#include <cmath>

namespace BiserLikeModel {

double MM2(double _v1, double _v2);
double MM(double _val, double _vmax, double _km);
double MM3(double l1, double o2, double beta);

double LaplacePolar0(double valc, double valp, double dr);
double LaplacePolar(double valm, double valc, double valp, double dr, double r);
double averageConcentration(double *array, double *points, double delta, int r_0p, int r_1, int r);


void SwapArrays(double **array1, double **array2);
void condition_assing(double *array1, double *array2, int length, double value);

void FillArray(double *array, double value, int from, int to);
void PrintArray(double *array, int length);
void concatenate_vals(double *, std::vector<double> *, int );

}  //  namespace BiserLikeModel
#endif  //  INCLUDE_BISERLIKEMODEL_UTILS_H_
