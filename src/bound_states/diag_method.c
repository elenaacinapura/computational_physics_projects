#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

#include "util.h"

#define POT 0

int main() {
	/* Environment variables */
	double R = 10.0;
	const int N = 970;
	double dx = 2.0 * R / (N - 1);

	double xi = 0.01;
	double a = 2.0 / xi;
	Params_cosh p;
	p.a = a;

	double diag[N];
	double subdiag[N - 1];
	double res_eigval[N];
	double res_eigvec[N][N];

	/* Fill diag and subdiag */
	for (int i = 0; i < N; i++) {
		double x = -R + i * dx;
		if (POT == 0) {
			diag[i] = V_cosh(x, &p) + xi/(dx*dx);
		}
		if (i < N - 1) {
			subdiag[i] = -xi/(2.0*dx*dx);
		}
	}

	int info = diagonalize_tridiag_double(N, diag, subdiag, (double*)res_eigvec, res_eigval);
	assert(info == 0);

    printf("------------------------------------------------\nDIAGONALIZATION METHOD\nSearching bound energies for xi = %.3lf\n\n", xi);
    int cnt_bound = 0;
	for (int i = 0; i < N; i++) {
		if (res_eigval[i] > 0.0) {
			break;
		}
		printf("E_%d = %g\n", i, res_eigval[i]);
        cnt_bound++;
	}
    printf("\nNumber of bound states found: %d\n", cnt_bound);
	printf("------------------------------------------------\n");
}