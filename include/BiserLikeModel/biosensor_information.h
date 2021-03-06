/*
 *  Copyright (c) Linas Petkevicius 2017
 *  Vilnius University
 *  GNU General Public license
 * */

#ifndef INCLUDE_BISERLIKEMODEL_BIOSENSOR_INFORMATION_H_
#define INCLUDE_BISERLIKEMODEL_BIOSENSOR_INFORMATION_H_
#include <BiserLikeModel/constants.h>

namespace BiserLikeModel {

enum resp_method {
    DEFAULT_TIME,  // 0 - iki pusiausvyros
    MIN_TIME,      // 1 - iki pusiausvyros su nurodytu minimaliu laiku
    FIXED_TIME,     // 2 - fiksuotas laikas
    HALF_TIME,     // 3 - iki duotos salygos laikas
    CML_TIME     // 4 - iki duotos komulatyvios salygos laikas
};

enum layer {
    MICROREACTOR,  // 0 - Mikro-reaktoriaus sluoksnis(Difuzija+MM)
    DIFFUSION,     // 1 - Difuzijos sluoksnis(Difuzija)
    BAUDARY        // 2 - Isorinis sluoksnis salytis su isore
};


struct layer_params {
    // Laukas nurodo ar tai fermento sluoksnis
    int enz_layer;

    // Difuzijos koeficientai (M/s)
    double Dl;
    double Dpr;
    double Do2;

    // Sluoksnio storis (m)
    double d;
};

struct bio_params {

    // Pusiausvyros konstantos (M)
    double km1, km2;
    double rho;

    //double kcat1, kcat2;
    // Greicio konstantos (M/s)
    double vmax1, vmax2;

    // Žingsnis pagal laiką (s)
    double dt;
    // Į kiek dalių dalinami sluoksniai
    int n;
    // Metodas, kuriuo bus nustatomas atsako laikas:
    enum resp_method resp_t_meth;
    // Minimalus atsako laikas (s)
    double min_t;
    // Fiksuotas atsako laikas (s)
    double resp_t;
    // Išvedimo failas
    char *out_file_name;
    bool write_to_file;

    bool dimensionless;
    bool spherical;

    // Elektronų, dalyvaujančių krūvio pernešime, skaičius
    int ne;

    // Pradinės koncentracijos tirpale (M)
    double pr_0;
    double l_0;
    double o2_0;

    // Biojutiklio sluoksnių skaičius
    int layer_count;
    //double oxigen;

    // Biojutiklio sluoksnių masyvas
    struct layer_params *layers;
};
}  //  namespace BiserLikeModel
#endif  //  INCLUDE_BISERLIKEMODEL_BIOSENSOR_INFORMATION_H_
