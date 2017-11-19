/*
 *  Copyright (c) Linas Petkevicius 2017
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

double Laplace0(double valc, double valp, double dr);
double Laplace(double valm, double valc, double valp, double dr);
double LaplaceSperical0(double valc, double valp, double dr);
double LaplaceSperical(double valm, double valc, double valp, double dr,
                       double r);

double averageConcentration(double *array, double *points, double delta,
                            int r_0p, int r_1, int r, bool scale);
double averageRate(double *array, double *points, double delta, int r_0,
                   int r_1, double vmax, double km);

void SwapArrays(double **array1, double **array2);
void condition_assing(double *array1, double *array2, int length, double value);

void FillArray(double *array, double value, int from, int to);
void PrintArray(double *array, int length);
void concatenate_vals(double *, std::vector<double> *, int );

double calc_L2(double *array1, double *array2, int length);

}  //  namespace BiserLikeModel
#endif  //  INCLUDE_BISERLIKEMODEL_UTILS_H_
