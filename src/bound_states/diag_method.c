#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

#include "util.h"

int main() {
	/* Set potential type */
	int POT;
	printf("\nWhich potential do you want?\n \t- cosh -> type 0\n\t- Lennard Jones -> type 1\n");
	scanf("%d", &POT);
	assert(POT == 0 || POT == 1);

	/* Environment variables */
	double R = 20.0;
	const int N = 970;
	double dx;
	double xi, a;
	int l;
	Params_cosh p_cosh;
	Params_lj p_lj;

	if (POT == 0) {
		printf("------------------------------------------------\nDIAGONALIZATION METHOD - COSH POTENTIAL \nSearching bound energies for xi = %.3lf\n\n", xi);

		xi = 0.01;
		a = 2.0 / xi;
		p_cosh.a = a;

		dx = 2.0 * R / (N - 1);

		double diag[N];
		double subdiag[N - 1];
		double res_eigval[N];
		double res_eigvec[N][N];

		/* Fill diag and subdiag */
		for (int i = 0; i < N; i++) {
			double x = -R + i * dx;
			diag[i] = V_cosh(x, &p_cosh) + xi / (dx * dx);
			if (i < N - 1) {
				subdiag[i] = -xi / (2.0 * dx * dx);
			}
		}

		int info = diagonalize_tridiag_double(N, diag, subdiag, (double*)res_eigvec, res_eigval);
		assert(info == 0);

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

	} else if (POT == 1) {
		printf("------------------------------------------------\nDIAGONALIZATION METHOD - LENNARD JONES POTENTIAL\nSearching bound energies for l = %d\n\n", l);

		double kB = 1.38e-23;
		double eps = 35.7 * kB;
		double s = 2.79e-10;
		double hbar = 1.05e-34;
		double m = 10 * 1.66e-27;

		xi = pow(hbar, 2) / (m * s * s * eps);
		a = 2.0 / xi;

		double r0 = 0.4;
		dx = (R - r0) / ((double)(N - 1));

		p_lj.a = a;

		int cnt_bound = 0;

		for (int l = 0; l < 10; l++) {
			p_lj.l = l;

			double diag[N];
			double subdiag[N - 1];
			double res_eigval[N];
			double res_eigvec[N][N];

			/* Fill diag and subdiag */
			for (int i = 0; i < N; i++) {
				double x = r0 + i * dx;
				diag[i] = V_lj(x, &p_lj) + xi / (dx * dx);
				if (i < N - 1) {
					subdiag[i] = -xi / (2.0 * dx * dx);
				}
			}

			int info = diagonalize_tridiag_double(N, diag, subdiag, (double*)res_eigvec, res_eigval);
			assert(info == 0);

			printf("l = %d\n", l);
			for (int i = 0; i < N; i++) {
				if (res_eigval[i] > 0.0) {
					break;
				}
				printf("E_%d = %g\n", i, res_eigval[i]);
				cnt_bound++;
			}
		}
		printf("\nNumber of bound states found: %d\n", cnt_bound);
		printf("------------------------------------------------\n");
	}
}