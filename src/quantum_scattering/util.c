#include "util.h"

#include <assert.h>
#include <complex.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <numerov.h>
#include <print_routines.h>
#include <stdio.h>

double F_hard(double r, void *param) {
	Params *p = (Params *)param;
	double xi = p->xi;
    double E = p->E;
	int l = p->l;
	return (double)(l*(l-1)/(r*r)) - xi*E;
}

void execute_numerov (double x[], double u[], int dim, double dx, double F (double, void *), void *p) {
    int i = 2; 
    while (i < dim) {
        double u_curr = numerov_1D(x[i-1], u[i-1], u[i-2], dx, F, p);
        u[i] = u_curr;
        x[i] = x[i-1] + dx;
        i++;
    }
}

double find_delta (double r1, double r2, double u1, double u2, int l, double k) {
    double K = u2*r1/(u1*r2);
    return atan((gsl_sf_bessel_jl(l, k*r2) - K*gsl_sf_bessel_jl(l, k*r1)) / (gsl_sf_bessel_yl(l, k*r2) - K * gsl_sf_bessel_yl(l, k*r1)));
}