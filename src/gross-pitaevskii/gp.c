#include <assert.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>

/*============== PARAMETERS =============*/
const double G = 5.0;
const double R = 5.0;
const int N = 1024;
const double alpha = 0.1;
double dr;

/* ===================== FUNCTION HEADERS ====================== */
double u_0(double r);
void normalize(double u[N]);
void init_u(double u[N]);
double Vscf(double u, double r);
double Vtot(double u, double r);
void fill_H(double u[N], double diag[N], double subdiag[N - 1]); /* given u */
double dist_L2(double u1[], double u2[]);

/* ===================== FUNCTION HEADERS ====================== */
int main() {
	dr = R / (N - 1);

	double u[N];
	init_u(u);

	double diag[N];
	double subdiag[N - 1];

	double dist, mu;
	int cnt = 0;

	do {
		cnt++;
		fill_H(u, diag, subdiag);

		double res_eigval[N];
		double res_eigvec[N][N];
		int info = diagonalize_tridiag_double(N, diag, subdiag, (double *)res_eigvec, res_eigval);
		assert(info == 0);

		mu = res_eigval[0];
		double u_new[N];
		for (int i = 0; i < N; i++) {
			u_new[i] = res_eigvec[0][i];
		}
		normalize(u_new);

		for (int i = 0; i < N; i++) {
			u_new[i] = alpha * u_new[i] + (1.0 - alpha) * u[i];
		}
		normalize(u_new);

		dist = dist_L2(u, u_new);
		printf("\rcnt = %d\t dist = %lf", cnt, dist);
		fflush(stdout);

		vec_copy(N, u_new, u);
	} while (dist > 1e-5);

	FILE *f;
	f = fopen("solution.csv", "w");
	fprintf(f, "0\t0\n");
	for (int i = 0; i < N; i++) {
		fprint_double(f, (i + 1) * dr);
		double r = (i == 0 ? dr : (i + 1) * dr);
		fprint_double_newline(f, fabs(u[i] / (sqrt(4.0 * M_PI) * r)));
	}
	fclose(f);

	printf("\nMu = ");
	fprint_double_newline(stdout, mu);
}

/* ===================== FUNCTIONS ====================== */
double u_0(double r) {
	return exp(-r * r) * sqrt(4.0 * M_PI) * r;
}
void normalize(double u[N]) {
	double No = 0.0;
	for (int i = 0; i < N; i++) {
		No += u[i] * u[i] * dr;
	}
	for (int i = 0; i < N; i++) {
		u[i] /= sqrt(No);
	}
}
void init_u(double u[N]) {
	for (int i = 0; i < N; i++) {
		double r = (i + 1) * dr;
		u[i] = u_0(r);
	}
	normalize(u);
}
double Vscf(double u, double r) {
	if (r < dr) {
		return 0;
	}
	return 2.0 * G * pow(u / r, 2);
}
double Vtot(double u, double r) {
	return r * r + Vscf(u, r);
}
void fill_H(double u[N], double diag[N], double subdiag[N - 1]) {
	for (int i = 0; i < N; i++) {
		double r = (i + 1) * dr;
		diag[i] = Vtot(u[i], r) + 2.0 / (dr * dr);
		if (i < N - 1) {
			subdiag[i] = -1.0 / (dr * dr);
		}
	}
}
double dist_L2(double u1[], double u2[]) {
	double res = 0.0;
	for (int i = 0; i < N; i++) {
		res += fabs(u1[i] - u2[i]) * dr;
	}
	return res;
}