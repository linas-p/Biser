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
    bio_info->l_0 = j["initial_conditions"]["laccase_0"];
    bio_info->o2_0 = j["initial_conditions"]["oxygen_0"];
	bio_info->rho = j["rho"];


    bio_info->layer_count = j["layers"].size();
    bio_info->layers = new layer_params[bio_info->layer_count];

    // Užpildoma sluoksnių informacija
    // 0
    bio_info->layers[0].enz_layer = j["layers"]["layer1"]["enzyme"];
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Dl = j["layers"]["layer1"]["Diff_laccase"];
    bio_info->layers[0].Do2 = j["layers"]["layer1"]["Diff_oxygen"];
    bio_info->layers[0].Dpr = j["layers"]["layer1"]["Diff_product"];
    // [um] -> [cm]
    bio_info->layers[0].d = j["layers"]["layer1"]["layer_length"];

    // 1
    bio_info->layers[1].enz_layer =  j["layers"]["layer2"]["enzyme"];
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Dl =  j["layers"]["layer2"]["Diff_laccase"];
    bio_info->layers[1].Do2 = j["layers"]["layer2"]["Diff_oxygen"];
    bio_info->layers[1].Dpr = j["layers"]["layer2"]["Diff_product"];
    // [um] -> [cm]
    bio_info->layers[1].d = j["layers"]["layer2"]["layer_length"];

    // 2
    bio_info->layers[2].enz_layer =  j["layers"]["layer3"]["enzyme"];
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[2].Dl =  j["layers"]["layer3"]["Diff_laccase"];
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


    // M
    bio_info->km1 = 9.6 * 1e-3;
    bio_info->km2 = 5 * 1e-4;

	// M/s
    bio_info->vmax1 = 1.9 * 1e-4;
    bio_info->vmax2 = 3.9 * 1e-4;


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

    // M
    bio_info->pr_0 = 0 * 1e-3;
    bio_info->l_0 = 2e-3;
    bio_info->o2_0 = 2.5 * 1e-4;

	bio_info->rho = 0.56;
    bio_info->layer_count = 3;
    bio_info->layers = new layer_params[bio_info->layer_count];

    // Užpildoma sluoksnių informacija
    // 0
    bio_info->layers[0].enz_layer = 1;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[0].Dl = 2.2 * 1e-6;
    bio_info->layers[0].Do2 = 0.8 * 1e-5;
    bio_info->layers[0].Dpr = 2.2 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[0].d = 0.025;

    // 1
    bio_info->layers[1].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[1].Dl = 6.7 * 1e-6;
    bio_info->layers[1].Do2 = 2.4 * 1e-5;
    bio_info->layers[1].Dpr = 6.7 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[1].d = 0.005;

    // 2
    bio_info->layers[2].enz_layer = 0;
    // [um^2/s] -> [cm^2/s]
    bio_info->layers[2].Dl = 6.7 * 1e-6;
    bio_info->layers[2].Do2 = 2.4 * 1e-5;
    bio_info->layers[2].Dpr = 6.7 * 1e-6;
    // [um] -> [cm]
    bio_info->layers[2].d = 0.054;
}

int main() {
    struct bio_params *bio_info = new bio_params;
    std::vector<double> P, L, O2, t, CP, CL, CO2;

#ifdef JSON
    json_fill(bio_info);
#else
    static_fill(bio_info);
#endif
    calculate_explicitly(bio_info, NULL, &callback_crunched, &P, &L, &O2, &t, &CL, &CP, &CO2);


    free(bio_info->layers);
    free(bio_info);

    return 0;
}
