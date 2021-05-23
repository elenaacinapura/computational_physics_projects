#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

/*================ STRUCTURES ================*/
typedef struct Params {
	double a, E; /* a = 2 m s^2 V0 / hbar^2 */
} Params;

/*================ FUNCTION HEADERS ================*/
double V(double x, void *param);

int main() {
	double R = 3.0;
	const int N = 1024/2 - 1;
	double dx;
	Params p;

	double xi_true = 0.015;
    double xi = 2.0 * xi_true;
	double a = 2.0 / xi;
	p.a = a;

	dx = R / (N - 1);

	/*=================== WELCOME =================*/
	printf("=================================================\n");
	printf("BOUND STATES WITH DIAGONALIZATION\n");
	printf("=================================================\n");
	printf("Parameters:\n");
	printf("\txi = hbar^2 / (2 m V0 a^2) = %.4lf\n", xi_true);
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
		diag[i] = V(x, &p) + xi / (dx * dx);
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
    printf("Printing eigenfunction to file.\n");
    FILE *f;
    f = fopen("eigenfunction_diag.csv", "w");
    // fprintf(f, "x\tpsi\n");
    double norm = 0.0;
    for (int i = 0; i < N; i++) {
        norm += res_eigvec[0][i] * res_eigvec[0][i] * dx;
    }
    for (int i = 0; i < N; i++) {
        res_eigvec[0][i] /= sqrt(norm);
    }
    for (int i = 0; i < N; i++) {
        fprint_double(f, -R + i * dx);
        fprint_double_newline(f, res_eigvec[0][i]);
    }
	printf("=================================================\n");
}

/*================ FUNCTIONS ================*/
double V(double x, void *param) {
	return -x*x * (x + 1);
}