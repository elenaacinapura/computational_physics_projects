#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gnuplot_i.h>
#include <integrate_notypesafe.h>
#include <print_routines.h>

#include "util.h"

double integrand (double x, void *p) {
    Params *params = (Params *) p;
    double a = params->alpha;
    double xi = params->xi;

    return sqrt(2.0/M_PI) / a * exp(-2.0 * x*x /(a*a)) * (2.0*xi / (a*a) - 4.0*xi*x*x/(pow(a, 4)) - pow(cosh(x), -4));
}

double energy (double a, void *p) {
    Params *params = (Params *) p;
    params->alpha = a;
    double cutoff = 50.0;
    return integrate(integrand, -cutoff, cutoff, 1000, params);
}

