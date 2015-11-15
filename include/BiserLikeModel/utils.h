/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */
#ifndef INCLUDE_BISERLIKEMODEL_UTILS_H_
#define INCLUDE_BISERLIKEMODEL_UTILS_H_

#include <stdio.h>

namespace BiserLikeModel {
void swap_arrays(double **array1, double **array2);
void condition_assing(double *array1, double *array2, int length, double value);

void fill_array(double *array, int length, double value, int from);
void print_array(double *array, int length);
}  //  namespace BiserLikeModel
#endif  //  INCLUDE_BISERLIKEMODEL_UTILS_H_
