##
##  Copyright (c) Linas Petkevicius 2016
##  Vilnius University
##  GNU General Public license
##

library(Rcpp)
cppFunction('double LaplacePolar(double valm, double valc, double valp, double dr, double r) {
    double val = (valp - 2 * valc + valm) / (dr * dr)
                 + (1 / r) * (valp - valm)/(dr);
    return val;

}')

cppFunction('double LaplacePolar0(double valc, double valp, double dr) {
    return (2 * (valp - valc) / std::pow(dr, 2));
}')

cppFunction('double MM(double val, double vmax, double km) {
    return (vmax * val) / (km + val);
}')

cppFunction('double MM2(double v1, double v2) {
    double val = 0;
    if((v1 || v2)) {
        val = (v1 * v2) /(v1 + 2 * v2);
    } else {
        val = 0;
    }
    return val;
}')


