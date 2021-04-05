#include "util.h"
#include <numerov.h>

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

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



double F (double f (double, void *), double x, void *p) {
    return 2 * (f(x, p) - E);
}


void numerov_step (int idx, double f (double, void *), void *p) {
    assert(idx >= 1);
    phi.v[idx+1] = (phi.v[idx] * (2.0 + 5.0/6.0 * dx*dx * F(f, X.v[idx], p)) - phi.v[idx - 1] * (1.0 - dx*dx / 12.0 * F(f, X.v[idx - 1], p)))/(1 - dx*dx/12.0 * F(f, X.v[idx + 1], p));
    phi.dim++;
}

void execute_numerov (double F (double, void *), double V (double, void *p), void *p, FILE *file, bool output) {
    /* initialize X */
	fill_X();
    /* initialize phi */
	phi.dim = 2;
	double k = 2.0 * E;
	phi.v[0] = cexp(-I * k * L);
	phi.v[1] = cexp(-I * k * (L + dx));

    int idx = 1;
	
    while (idx < X.dim - 1) {
		complex double phi_next = numerov_1D(X.v[idx], phi.v[idx], phi.v[idx-1], dx, F_square, p);
		phi.v[idx+1] = phi_next;
		phi.dim++;
        if (output) {
            fprint_double(file, X.v[idx]);
            fprint_double(file, V_square(X.v[idx], &(Empty_struct){}));
            fprint_double(file, creal(phi.v[idx]));
            fprint_double(file, cimag(phi.v[idx]));
            fprintf(file, "\n");
        }
		idx++;
	}
}

complex double calculate_T (int idx1, int idx2) {
    double k = 2.0*E;
    complex double x1 = X.v[idx1];
    complex double x2 = X.v[idx2];
    return (cexp(I*k*(x2-x1)) - cexp(-I*k*(x2-x1)))/(phi.v[idx1]*cexp(I*k*x2) - phi.v[idx2]*cexp(I*k*x1));
} 

complex double calculate_R (int idx1, int idx2) {
    double k = 2.0*E;
    complex double x1 = X.v[idx1];
    complex double x2 = X.v[idx2];
    return (cexp(-I*k*(x2))*phi.v[idx1] - cexp(-I*k*x1)*phi.v[idx2])/(phi.v[idx1]*cexp(I*k*x2) - phi.v[idx2]*cexp(I*k*x1));
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

