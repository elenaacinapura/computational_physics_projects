#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

#include "util.h"

int main() {
	/* Environment variables */
	double R = 1.0;
	const int N = 400;
	double dx = R / (N - 1);
	double xi, a;

	xi = 2.0;
	a = 2.0 / xi;

	Params_periodic p;
	p.a = a;

	/* Matrix */
	complex double H[N][N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				double x = i * dx;
				H[i][j] = V_periodic(x, &p) + xi / (dx * dx);
			} else if (j == i + 1 || j == i - 1) {
				H[i][j] = -xi / (2.0 * dx * dx);
			} else {
				H[i][j] = 0.0;
			}
		}
	}
    double res_eigval[N];
    complex double res_eigvec[N][N];

	FILE *f;
	f = fopen("periodic.csv", "w");

	double K = -M_PI;
	double dK = 0.1;
	double K_end = M_PI;

	while (K < K_end) {
		H[0][N - 1] = -xi / (2.0 * dx * dx) * cexp(-I * K);
		H[N - 1][0] = -xi / (2.0 * dx * dx) * cexp(I * K);


		diagonalize_herm(N, (_Complex double *)H, (_Complex double *)res_eigvec, res_eigval);

		fprint_double(f, K);
		for (int i = 0; i < 3; i++) {
			fprint_double(f, res_eigval[i]);
		}
		fprintf(f, "\n");

		K += dK;
	}
    fclose (f);
}