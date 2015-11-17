/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */

#include <Rcpp.h>
#include <json.hpp>
#include <BiserLikeModel/biosensor_information.h>
#include <BiserLikeModel/explicit_calculator.h>
#include <fstream>

using namespace BiserLikeModel;
using json = nlohmann::json;

void callback_crunched(void *ptr, int time) {
    printf("%ds simulated\n", time);
}

void json_fill(struct bio_params *bio_info, \
               std::string configs = "../config/params.json") {
    json j;
    std::ifstream ifs (configs, std::ifstream::in);
    ifs >> j;
    //std::cout << j.dump(4) << std::endl;

    // [s^-1]
    bio_info->k1 = j["equations_params"]["k1"];
    bio_info->k2 = j["equations_params"]["k2"];

    // [mol/l] -> [mol/cm^3]
    bio_info->km1 = j["equations_params"]["km_1"];
    bio_info->km2 = j["equations_params"]["km_2"];

    // [mol/l] -> [mol/cm^3]
    bio_info->vmax1 = j["equations_params"]["vmax_1"];
    bio_info->vmax2 = j["equations_params"]["vmax_2"];

    // [s]
    bio_info->dt = j["dt"];
    bio_info->n = j["n"];
    bio_info->resp_t_meth = FIXED_TIME;
    // [s]
    bio_info->min_t = 100;
    // [s]
    bio_info->resp_t = 20;

    bio_info->out_file_name = "output.dat";
    bio_info->write_to_file = false;
    bio_info->ne = 1;

    // [mol/l] -> [mol/cm^3]
    bio_info->pr_0 = j["initial_conditions"]["product_0"];
    bio_info->g_0 = j["initial_conditions"]["glucose_0"];
    bio_info->o2_0 = j["initial_conditions"]["oxygen_0"];
    bio_info->e1ox_0 = j["initial_conditions"]["e1_oxi_0"];
    bio_info->e2ox_0 = j["initial_conditions"]["e2_oxi_0"];
    bio_info->e1red_0 = j["initial_conditions"]["e1_red_0"];
    bio_info->e2red_0 = j["initial_conditions"]["e2_red_0"];

    bio_info->layer_count = j["layers"].size();
    bio_info->layers = new layer_params[bio_info->layer_count];

    // Užpildoma sluoksnių informacija
    // 0
    bio_info->layers[0].enz_layer = j["layers"]["layer1"]["enzyme"];
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Dg = j["layers"]["layer1"]["Diff_glucose"];
    bio_info->layers[0].Do2 = j["layers"]["layer1"]["Diff_oxygen"];
    bio_info->layers[0].Dpr = j["layers"]["layer1"]["Diff_product"];
    // [um] -> [cm]
    bio_info->layers[0].d = j["layers"]["layer1"]["layer_length"];

    // 1
    bio_info->layers[1].enz_layer =  j["layers"]["layer2"]["enzyme"];
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Dg =  j["layers"]["layer2"]["Diff_glucose"];
    bio_info->layers[1].Do2 = j["layers"]["layer2"]["Diff_oxygen"];
    bio_info->layers[1].Dpr = j["layers"]["layer2"]["Diff_product"];

    // [um] -> [cm]
    bio_info->layers[1].d = j["layers"]["layer2"]["layer_length"];

}

RcppExport SEXP calculate(SEXP x) {
    Rcpp::NumericVector params(x);
    Rcpp::NumericVector nv;
    bio_params* bio_info = new bio_params();
    json_fill(bio_info);

    std::vector<double> P, G, O2, Ox1, Ox2, Red1, Red2;
    std::clock_t start = std::clock();
    calculate_explicitly(bio_info, NULL, &callback_crunched, &P, &G, &O2, \
                         &Ox1, &Ox2, &Red1, &Red2);
    double time = (std::clock()-start)/ static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "\n all time " << time << std::endl;
    std::cout<< "Done!" << std::endl;

    Rcpp::NumericVector PP(P.begin(), P.end());
    Rcpp::NumericVector GG(G.begin(), G.end());
    Rcpp::NumericVector OO2(O2.begin(), O2.end());
    Rcpp::NumericVector OOx1(Ox1.begin(), Ox1.end());
    Rcpp::NumericVector OOx2(Ox2.begin(), Ox2.end());
    Rcpp::NumericVector RRed1(Red1.begin(), Red1.end());
    Rcpp::NumericVector RRed2(Red2.begin(), Red2.end());

    free(bio_info->layers);
    free(bio_info);

    return( Rcpp::List::create(Rcpp::Named("P")=PP,
                               Rcpp::Named("G")=GG,
                               Rcpp::Named("O2")=OO2,
                               Rcpp::Named("Ox1")=OOx1,
                               Rcpp::Named("Ox2")=OOx2,
                               Rcpp::Named("Red1")=RRed1,
                               Rcpp::Named("Red2")=RRed2));
}
