#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

/*================ STRUCTURES ================*/
typedef struct Params_periodic {
	double a, E;
} Params_periodic;

/*================ FUNCTION HEADERS ================*/
double V_periodic(double x, void *param);

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

	/*=================== WELCOME =================*/
	printf("===============================================================\n");
	printf("BOUND STATES OF PERIODIC POTENTIAL WITH DIAGONALIZATION\n");
	printf("===============================================================\n");
	printf("Parameters:\n");
	printf("\tN = %d\n", N);
	printf("\tR = %.1lf\n", R);
	printf("\txi = %.1lf\n", xi);
	printf("Starting calculating...\n\n");

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
	fprintf(f, "K\tE0\tE1\tE2\n");

	double K = -M_PI;
	double dK = 0.1;
	double K_end = M_PI;

	while (K < K_end) {
		H[0][N - 1] = -xi / (2.0 * dx * dx) * cexp(-I * K);
		H[N - 1][0] = -xi / (2.0 * dx * dx) * cexp(I * K);

		diagonalize_herm(N, (_Complex double *)H, (_Complex double *)res_eigvec, res_eigval);

		fprint_double(f, K);
		fprint_double(f, res_eigval[0]);
		fprint_double(f, res_eigval[1]);
		fprint_double_newline(f, res_eigval[2]);

		K += dK;
	}
	fclose(f);
	printf("Calculations ended successfully!\n");
	printf("===============================================================\n");
}

/*================ FUNCTIONS ================*/
double V_periodic(double x, void *param) {
	if (x < 0.3) {
		return -1.0;
	}
	return 0.0;
}