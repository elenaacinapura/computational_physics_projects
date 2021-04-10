#include "util.h"
#include <verlet.h>
#include <lennard_jones.h> /* contains polar r */

#include <stdio.h>
#include <math.h>

double Fx (double x, double z, void *p) {
    Params *param = (Params *) p;
    double theta = param->theta;
    double denom = pow(r_polar(x - sin(theta), 0.0, z - cos(theta)), 3);
    return (x - sin(theta)) / (denom) - x / pow(r_polar(x, 0.0, z), 3);
}

double Fz (double x, double z, void *p) {
    Params *param = (Params *) p;
    double theta = param->theta;
}