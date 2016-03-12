/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */

#ifndef INCLUDE_BISERLIKEMODEL_BIOSENSOR_INFORMATION_H_
#define INCLUDE_BISERLIKEMODEL_BIOSENSOR_INFORMATION_H_

namespace BiserLikeModel {

enum resp_method {
    DEFAULT_TIME,  // 0 - iki pusiausvyros
    MIN_TIME,      // 1 - iki pusiausvyros su nurodytu minimaliu laiku
    FIXED_TIME     // 2 - fiksuotas laikas
};

struct layer_params {
    // Laukas nurodo ar tai fermento sluoksnis
    int enz_layer;

    // Difuzijos koeficientai (cm^2/s)
    double Dg;
    double Dpr;
    double Do2;

    // Sluoksnio storis (cm)
    double d;
};

struct bio_params {

    // Reakcijos greičio konstanta k2 (s^-1)
    double k1, k2;
    // Pusiausvyros konstantos (mol/cm^3)
    double km1, km2;

	double kcat1, kcat2;
    // double vmax1, vmax2;

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
    // Elektronų, dalyvaujančių krūvio pernešime, skaičius
    int ne;

    // Pradinės koncentracijos tirpale (mol/cm^3)
    double pr_0;
    double g_0;
    double o2_0;
    double e1ox_0;
    double e2ox_0;
    double e1red_0;
    double e2red_0;
    // Biojutiklio sluoksnių skaičius
    int layer_count;

    // Biojutiklio sluoksnių masyvas
    struct layer_params *layers;
};
}  //  namespace BiserLikeModel
#endif  //  INCLUDE_BISERLIKEMODEL_BIOSENSOR_INFORMATION_H_
