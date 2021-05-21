#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

/*================ STRUCTURES ================*/
typedef struct Params_cosh {
	double a, E; /* a = 2 m s^2 V0 / hbar^2 */
} Params_cosh;

/*================ FUNCTION HEADERS ================*/
double V_cosh(double x, void *param);

int main() {
	double R = 20.0;
	const int N = 970;
	double dx;
	Params_cosh p_cosh;

	double xi = 0.01;
	double a = 2.0 / xi;
	p_cosh.a = a;

	dx = 2.0 * R / (N - 1);

	/*=================== WELCOME =================*/
	printf("=================================================\n");
	printf("BOUND STATES WITH DIAGONALIZATION\n");
	printf("=================================================\n");
	printf("Parameters:\n");
	printf("\tV(x) = -V0*cosh(x/a)^(-4)\n");
	printf("\txi = hbar^2 / (m V0 a^2) = %.4lf\n", xi);
	printf("\tN = %d\n", N);
	printf("\tR = %.1lf\n", R);
	printf("Starting calculating...\n\n");

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

	int info = diagonalize_tridiag_double(N, diag, subdiag, (double *)res_eigvec, res_eigval);
	assert(info == 0);

	printf("Starting calculating...\n\n");

	int cnt_bound = 0;
	for (int i = 0; i < N; i++) {
		if (res_eigval[i] > 0.0) {
			break;
		}
		printf("E_%d = %g\n", i, res_eigval[i]);
		cnt_bound++;
	}
	printf("\nNumber of bound states found: %d\n", cnt_bound);
	printf("=================================================\n");
}

/*================ FUNCTIONS ================*/
double V_cosh(double x, void *param) {
	return -1.0 / (pow(cosh(x), 4));
}