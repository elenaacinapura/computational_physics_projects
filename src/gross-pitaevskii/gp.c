#include <assert.h>
#include <linear_algebra/blas_wrappers.h>
#include <linear_algebra/lapack_wrappers.h>
#include <math.h>
#include <print_routines.h>

/*============== PARAMETERS =============*/
const double G_target = 5.0;
const double R = 5.0;
const int N = 1024;
const double alpha = 0.1;
double dr, G;

/* ===================== FUNCTION HEADERS ====================== */
double u_0(double r);
void normalize(double u[N]);
void init_u(double u[N]);
double Vscf(double u, double r);
double Vtot(double u, double r);
void fill_H(double u[N], double diag[N], double subdiag[N - 1]); /* given u */
double dist_L2(double u1[], double u2[]);
double calculate_ex_energy(double u[]);
double calculate_V(double u[]);

/* ===================== FUNCTION HEADERS ====================== */
int main() {
	dr = R / (N - 1);

	/* Welcome */
	printf("=====================================================\n");
	printf("GROSS-PYTAEVSKII'S METHOD\n");
	printf("=====================================================\n");
	printf("Parameters:\n");
	printf("\tN = %d\n", N);
	printf("\tR = %.1lf\n", R);
	printf("\talpha (aka \"gentleness\") = %.3lf\n", alpha);
	printf("\tG_target = a * N / lambda_H = %.1lf\n", G_target);
	printf("=====================================================\n");
	printf("Calculating...\n");

	double u[N];
	init_u(u);

	double diag[N];
	double subdiag[N - 1];

	double dist, mu;
	int cnt = 0;

	if (G_target > 20.0) {
		double G_start = 1.0;
		double dG = 1.0;
		G = G_start;
		while (G <= G_target + 0.1 * dG) {
			fflush(stdout);
			printf("\rCurrent G = %.1lf", G);
			fflush(stdout);
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

				vec_copy(N, u_new, u);
			} while (dist > 1e-5);
			G += dG;
		}
	} else {
		G = G_target;
		do {
			fflush(stdout);
			printf("\rDistance from previous approximation = %lf", dist);
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

			vec_copy(N, u_new, u);
		} while (dist > 1e-5);
	}

	double V_ext = calculate_ex_energy(u);
	double V_tot = calculate_V(u);
	double K = mu - V_tot;

	FILE *f;
	f = fopen("solution.csv", "w");
	fprintf(f, "0\t0\n");
	for (int i = 0; i < N; i++) {
		fprint_double(f, (i + 1) * dr);
		double r = (i == 0 ? dr : (i + 1) * dr);
		fprint_double_newline(f, fabs(u[i] / (sqrt(4.0 * M_PI) * r)));
	}
	fclose(f);

	printf("\n=====================================================\n");
	printf("Results:\n");
	printf("\tMu = %lf\n", mu);
	printf("\t<V_tot> = %lf\n", V_tot);
	printf("\t<V_ext> = %lf\n", V_ext);
	printf("\t<K> = %lf\n", K);
	printf("=====================================================\n");
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
double calculate_ex_energy(double u[]) {
	double res = 0.0;
	for (int i = 0; i < N; i++) {
		double r = (i + 1) * dr;
		double phi = u[i] / (sqrt(4 * M_PI) * r);
		res += 4 * M_PI * phi * phi * pow(r, 4) * dr;
	}
	return res;
}
double calculate_V(double u[]) {
	double ex_V = calculate_ex_energy(u);
	double int_V = 0.0;
	for (int i = 0; i < N; i++) {
		double r = (i + 1) * dr;
		double phi = u[i] / (sqrt(4 * M_PI) * r);
		double V = 8 * M_PI * G * phi * phi;
		int_V += 4 * M_PI * phi * phi * V * r * r * dr;
	}
	return ex_V + int_V;
}