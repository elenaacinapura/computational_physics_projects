#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

/*================ CONSTANTS ================*/
const double kB = 1.38e-23;
const double eps = 35.7 * kB;
const double s = 2.79e-10;
const double hbar = 1.05e-34;
const double m = 10 * 1.66e-27;

/*================ STRUCTURES ================*/
typedef struct Params_lj {
	double a, E;
	int l;
} Params_lj;

/*================ FUNCTION HEADERS ================*/
double V_lj(double r, void *param);

int main() {
	double R = 5.0;
	const int N = 970;
	double r0 = 0.4;
	double dx = (R - r0) / ((double)(N - 1));
	double xi = pow(hbar, 2) / (m * s * s * eps);
	double a = 2.0 / xi;

	Params_lj p_lj;
	p_lj.a = a;

	/*=================== WELCOME =================*/
	printf("=================================================\n");
	printf("BOUND STATES WITH DIAGONALIZATION\n");
	printf("=================================================\n");
	printf("Parameters:\n");
	printf("\tV(x) = Lennard-Jones for Neon 10\n");
	printf("\tN = %d\n", N);
	printf("\tR = %.1lf\n", R);
	printf("\tr0 = %.1lf\n", r0);
	printf("Starting calculating...\n\n");

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

		int info = diagonalize_tridiag_double(N, diag, subdiag, (double *)res_eigvec, res_eigval);
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
	printf("=================================================\n");
}

/*================ FUNCTIONS ================*/
double V_lj(double r, void *param) {
	Params_lj *p = (Params_lj *)param;
	int l = p->l;
	double a = p->a;
	double xi = 2.0 / a;

	return 4.0 * (pow(r, -12) - pow(r, -6)) + xi * (double)(l * (l + 1)) / (r * r);
}