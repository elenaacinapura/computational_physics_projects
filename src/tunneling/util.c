#include <numerov.h>
#include <print_routines.h>
#include "util.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

double E;

double F_square (double x, void *p) {
    Param_F *param = (Param_F *) p;
    double xi = param->xi;

    if (fabs(x) < 0.5) {
        return xi * (1.0 - E);
    }
    else return -xi*E;
}

double F_gauss (double x, void *p) {
    Param_F *param = (Param_F *) p;
    double xi = param->xi;

    return xi * (-exp(-x*x/2) - E);
}

double F_asymm_L (double x, void *p) {
    Param_F *param = (Param_F *) p;
    double xi = param->xi;

    if (fabs(x) < 0.5) {
        if (x < 0.0) {
            return xi * (1.0 - E);
        } else {
            return xi * (0.5 - E);
        }
    } else return -xi*E;
}

double F_asymm_R (double x, void *p) {
    Param_F *param = (Param_F *) p;
    double xi = param->xi;

    if (fabs(x) < 0.5) {
        if (x > 0.0) {
            return xi * (1.0 - E);
        } else {
            return xi * (0.5 - E);
        }
    } else return -xi*E;
}

double F_cosh (double x, void *p) {
    Param_F *param = (Param_F *) p;
    double xi = param->xi;
    return xi*( pow(cosh(x),-4) - E);
    
}

void solve_numerov (double x[], complex double phi[], int dim, double dx, double F (double, void *p), void *p, bool printoutput, FILE *outfile) {
    /* Assuming x and phi have initial conditions in position 0 and 1 */

    Param_F *param = (Param_F *) p;
    double xi = param->xi;

    int i = 2;
    while (i < dim) {
        /* numerov */
        complex double new_phi = numerov_1D(x[i-1], phi[i-1], phi[i-2], dx, F, p);
        phi[i] = new_phi;
        x[i] = x[i-1] + dx;
        /* print values */
        if (printoutput) {
            fprint_double(outfile, -x[i]);
            fprint_double(outfile, creal(phi[i]));
            fprint_double(outfile, cimag(phi[i]));
            fprint_double(outfile, F(x[i], p)/xi + E);
            fprintf(outfile, "\n");
        }
        i++;
    }
}

complex double T_coeffcient(double x[], complex double phi[], int dim, double k) {
    int i1 = dim-1;
    int i2 = dim -10;
    double x1 = x[i1];
    double x2 = x[i2];
    complex double phi1 = phi[i1];
    complex double phi2 = phi[i2];

    return (cexp(I*k*(x2 - x1)) - cexp(I*k*(x1 - x2))) / (phi1 * cexp(I*k*x2) - phi2 * cexp(I*k*x1));
}
