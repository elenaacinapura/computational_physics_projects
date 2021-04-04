#include "util.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

Array phi, X;
double dx, E, L;

void fill_X () {
    int dim_e = (int)ceil(2.0*L/dx);
    assert(dim_e < MAXDIM);
    X.v[0] = -L;
    X.dim = 1;
    for (int i = 1; i < dim_e; i++) {
        X.v[i] = X.v[i-1] + dx;
        X.dim++;
    }
}

double V_square (double x, void *p) {
    Param_V_square *param = (Param_V_square *) p;
    double V0 = param->V0;
    double a = param->a;
    if (fabs(x) < a/2.0) {
        return V0;
    }
    return 0.0;
}

double F (double f (double, void *), double x, void *p) {
    return 2 * (f(x, p) - E);
}

void numerov_step (int idx, double f (double, void *), void *p) {
    assert(idx >= 1);
    phi.v[idx+1] = (phi.v[idx] * (2.0 + 5.0/6.0 * dx*dx * F(f, X.v[idx], p)) - phi.v[idx - 1] * (1.0 - dx*dx / 12.0 * F(f, X.v[idx - 1], p)))/(1 - dx*dx/12.0 * F(f, X.v[idx + 1], p));
    phi.dim++;
}

double calculate_T (int idx1, int idx2) {
    double k = 2.0*E;
    double x1 = X.v[idx1];
    double x2 = X.v[idx2];
    return (cexp(I*k*(x2-x1)) - cexp(-I*k*(x2-x1)))/(phi.v[idx1]*cexp(I*k*x2) - phi.v[idx2]*cexp(I*k*x1));
} 

/* PRINT AND CALCULATION ROUTINES */

double square_cabs (complex double z) {
    return cabs(z)*cabs(z);
}

void fprint_vec (FILE *file, double v [], int dim) {
    for (int i = 0; i < dim; i++) {
        fprintf(file, "%lf\n", v[i]);
    }
    fprintf(file, "\n");
}

void fprint_double (FILE *file, double d) {
    fprintf(file, "%lf\t", d);
}
