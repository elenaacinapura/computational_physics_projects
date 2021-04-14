#include "util.h"

#include <complex.h>
#include <math.h>
#include <numerov.h>
#include <print_routines.h>


double F_cosh (double x, void *param) {
    Params *p = (Params *) param;
    double a = p->a;
    double E = p->E;
    return -a * (pow(cosh(x), -4) + E);
}

void execute_numerov(double x[], double phi [], double dx, int dim, double F (double, void*), void *p) {
    int i = 2;
    do {
        phi[i] = numerov_1D(x[i-1], phi[i-1], phi[i-2], dx, F, p);
        x[i] = x[i-1] + dx;
        i++;
    } while (i < dim);
}

double calculate_delta (double x0, double phiF[], double phiB[], double dx, int dimF, int dimB, double F (double, void *), void *p) {
    return 1.0/dx * (phiF[dimF - 2] + phiB[dimB - 2] - phiF[dimB-1]*(2.0 + dx*F(x0, p))); 
}