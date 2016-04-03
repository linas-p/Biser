/*
 *  Copyright (c) Linas Petkevicius 2016
 *  Vilnius University
 *  GNU General Public license
 * */

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <BiserLikeModel/biosensor_information.h>
#include <BiserLikeModel/explicit_calculator.h>


using namespace BiserLikeModel;

#ifdef JSON
#include <json.hpp>
using json = nlohmann::json;

void json_fill(struct bio_params *bio_info, std::string configs = "../config/params.json") {
    json j;
    std::ifstream ifs (configs, std::ifstream::in);
    ifs >> j;

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
	bio_info->alpha = j["alpha"];


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

    // 2
    bio_info->layers[2].enz_layer =  j["layers"]["layer3"]["enzyme"];
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[2].Dg =  j["layers"]["layer3"]["Diff_glucose"];
    bio_info->layers[2].Do2 = j["layers"]["layer3"]["Diff_oxygen"];
    bio_info->layers[2].Dpr = j["layers"]["layer3"]["Diff_product"];
    // [um] -> [cm]
    bio_info->layers[2].d = j["layers"]["layer3"]["layer_length"];

}
#endif

void callback_crunched(void *ptr, int time) {
    printf("%ds simulated\n", time);
}

void static_fill(struct bio_params *bio_info) {


    // [mol/l] -> [mol/cm^3]
    bio_info->km1 = 6.8 * 1e-3;
    bio_info->km2 = 6.8 * 1e-3;

	// [mol/l] -> [mol/cm^3]
    bio_info->vmax1 = 4 * 1e-5;
    bio_info->vmax2 = 4 * 1e-5;


    // [s]
    bio_info->dt = 1e-2;
    bio_info->n = 4;
    bio_info->resp_t_meth = FIXED_TIME;

	// [s]
    bio_info->min_t = 100;

	// [s]
    bio_info->resp_t = 0;

    bio_info->out_file_name = "output.dat";
    bio_info->write_to_file = false;

    bio_info->ne = 1;

    // [mol/l] -> [mol/cm^3]
    bio_info->pr_0 = 0 * 1e-3;
    bio_info->g_0 = 1e-3;
    bio_info->o2_0 = 2.5 * 1e-4;

	bio_info->alpha = 0.5;
    bio_info->layer_count = 3;
    bio_info->layers = new layer_params[bio_info->layer_count];

    // Užpildoma sluoksnių informacija
    // 0
    bio_info->layers[0].enz_layer = 1;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Dg = 2.2 * 1e-6;
    bio_info->layers[0].Do2 = 0.8 * 1e-5;
    bio_info->layers[0].Dpr = 2.2 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[0].d = 0.10;

    // 1
    bio_info->layers[1].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Dg = 6.8 * 1e-6;
    bio_info->layers[1].Do2 = 2.4 * 1e-5;
    bio_info->layers[1].Dpr = 6.8 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[1].d = 0.02;

    // 2
    bio_info->layers[2].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[2].Dg = 6.8 * 1e-6;
    bio_info->layers[2].Do2 = 2.4 * 1e-5;
    bio_info->layers[2].Dpr = 6.8 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[2].d = 0.03;
}

int main() {
    struct bio_params *bio_info = new bio_params;
    std::vector<double> P, G, O2, t, Ct_g, Ct_p, Ct_o2;

#ifdef JSON
    json_fill(bio_info);
#else
    static_fill(bio_info);
#endif
    calculate_explicitly(bio_info, NULL, &callback_crunched, &P, &G, &O2, &t, &Ct_g, &Ct_p, &Ct_o2);


    free(bio_info->layers);
    free(bio_info);

    return 0;
}
