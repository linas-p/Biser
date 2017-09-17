/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */

#ifndef INCLUDE_BISERLIKEMODEL_EXPLICIT_CALCULATOR_H_
#define INCLUDE_BISERLIKEMODEL_EXPLICIT_CALCULATOR_H_
#include <vector>

namespace BiserLikeModel {

void calculate_explicitly(struct bio_params *bio_info, void *ptr, \
                          void (*callback_crunched)(void *, int),
                          std::vector<double> * points, \
                          std::vector<double> * P, \
                          std::vector<double> * L, \
                          std::vector<double> * O2, \
                          std::vector<double> * t, \
                          std::vector<double> * Ct_g,\
                          std::vector<double> * Ct_p, \
                          std::vector<double> * Ct_o2
                          );

void two_layer_model(struct bio_params *bio_info, void *ptr, \
                          void (*callback_crunched)(void *, int),
                          std::vector<double> * points, \
                          std::vector<double> * P, \
                          std::vector<double> * L, \
                          std::vector<double> * t, \
                          std::vector<double> * Ct_g,\
                          std::vector<double> * Ct_p,
                          std::vector<double> * Ot_p,\
                          std::vector<double> * Chr
                          );

}  //  namespace BiserLikeModel
#endif  //  INCLUDE_BISERLIKEMODEL_EXPLICIT_CALCULATOR_H_
