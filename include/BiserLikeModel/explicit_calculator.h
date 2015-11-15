/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */

#ifndef INCLUDE_BISERLIKEMODEL_EXPLICIT_CALCULATOR_H_
#define INCLUDE_BISERLIKEMODEL_EXPLICIT_CALCULATOR_H_

namespace BiserLikeModel {

void calculate_explicitly(struct bio_params *bio_info, void *ptr, \
                          void (*callback_crunched)(void *, int));
}  //  namespace BiserLikeModel
#endif  //  INCLUDE_BISERLIKEMODEL_EXPLICIT_CALCULATOR_H_
