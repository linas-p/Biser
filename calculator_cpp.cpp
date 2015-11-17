/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <BiserLikeModel/biosensor_information.h>
#include <BiserLikeModel/explicit_calculator.h>
#include <json.hpp>

using namespace BiserLikeModel;
using json = nlohmann::json;

void callback_crunched(void *ptr, int time) {
    printf("%ds simulated\n", time);
}

void json_fill(struct bio_params *bio_info, std::string configs = "../config/params.json") {
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
    bio_info->write_to_file = true;
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


void static_fill(struct bio_params *bio_info) {

    // [s^-1]
    bio_info->k2 = 2*1e-3;
    bio_info->k1 = 4*1e-3;

    // [mol/l] -> [mol/cm^3]
    bio_info->km1 = 8 * 1e-5;
    bio_info->km2 = 2 * 1e-5;
    // [mol/l] -> [mol/cm^3]
    bio_info->vmax1 = 4 * 1e-5;
    bio_info->vmax2 = 1e-5;

    // [s]
    bio_info->dt = 1e-2;
    bio_info->n = 10;
    bio_info->resp_t_meth = MIN_TIME;
    // [s]
    bio_info->min_t = 100;
    // [s]
    bio_info->resp_t = 0;
    bio_info->out_file_name = "output.dat";
    bio_info->ne = 1;

    // [mol/l] -> [mol/cm^3]
    bio_info->pr_0 = 0 * 1e-3;
    bio_info->g_0 = 1e-5;
    bio_info->o2_0 = 1.2 * 1e-6;
    bio_info->e1ox_0 = 1e-9;
    bio_info->e2ox_0 = 1e-9;
    bio_info->e1red_0 = 0*1e-9;
    bio_info->e2red_0 = 0*1e-9;

    bio_info->layer_count = 2;
    bio_info->layers = new layer_params[bio_info->layer_count];

    // Užpildoma sluoksnių informacija
    // 0
    bio_info->layers[0].enz_layer = 1;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Dg = 2.1 * 1e-6;
    bio_info->layers[0].Do2 = 0.666 * 1e-5;
    bio_info->layers[0].Dpr = 2.1 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[0].d = 0.1;

    // 1
    bio_info->layers[1].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Dg = 6.3 * 1e-6;
    bio_info->layers[1].Do2 = 2 * 1e-5;
    bio_info->layers[1].Dpr = 6.3 * 1e-6;

    // [um] -> [cm]
    bio_info->layers[1].d = 0.9;

}

int main() {
    struct bio_params *bio_info = new bio_params;
    std::vector<double> P, G, O2, Ox1, Ox2, Red1, Red2;
    // static_fill(bio_info);
    json_fill(bio_info);
    calculate_explicitly(bio_info, NULL, &callback_crunched, &P, &G, &O2, \
                         &Ox1, &Ox2, &Red1, &Red2);


    free(bio_info->layers);
    free(bio_info);

    return 0;
}
