#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <math.h>
#include <numerical_methods/fft_radial.h>
#include <physics/lennard_jones.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>

/*============== CONSTANTS ==============*/
const int N = 1024;
const double R = 30.0;
const double rho = 0.7;
const double T = 0.5;
const double alpha = 0.1;
const int N_REPS = 10;

/*============== FUNCTIONS HEADERS ==============*/
double h0(double r);
double init(double r[], double h[], double c[]);
double norm(double h_new[], double h_old[]);

/*============== MAIN ==============*/
int main() {
	double r[N];
	double h[N];
	double c[N];
	double F[N];

	init(r, h, c);

	for (int cnt = 0; cnt < N_REPS; cnt++) {
		double c_t[N], h_t[N], c_new[N], h_new[N]; /* transforms and new vectors */
		vec_copy(N, c, c_t);
		vec_copy(N, h, h_t);

		fft_radial_forward(h_t, N, R);
		fft_radial_forward(c_t, N, R);

		for (int i = 0; i < N; i++) {
			F[i] = h_t[i] * c_t[i] * rho * 2.0 * M_PI * R / N;
			if (isnan(F[i])) {
				printf("Nan value found in rep %d\n", cnt);
				exit(0);
			}
		}

		fft_radial_inverse(F, N, R);

		for (int i = 0; i < N; i++) {
			c_new[i] = (F[i] + 1.0) * h0(r[i]);
			h_new[i] = c_new[i] + F[i];

			c_new[i] = alpha * c_new[i] + (1.0 - alpha) * c[i];
			h_new[i] = alpha * h_new[i] + (1.0 - alpha) * h[i];
		}
		double nor = norm(h_new, h);
		// fprint_double_newline(stdout, nor);

		vec_copy(N, h_new, h);
		vec_copy(N, c_new, c);
	}

	FILE* file;
	file = fopen("percus.csv", "w");
	fprintf(file, "r\tg\n");

	for (int i = 0; i < N; i++) {
		fprint_double(file, r[i]);
		fprint_double_newline(file, h[i] + 1.0);
	}
	fclose(file);
}

/*============== FUNCTIONS ==============*/
double h0(double r) {
	if (r < 1e-5) {
		return -1.0;
	}
	double V = 4.0 * (pow(r, -12) - pow(r, -6));
	return exp(-V / T) - 1.0;
}
double init(double r[], double h[], double c[]) {
	for (int n = 0; n < N; n++) {
		r[n] = n * R / N;
		h[n] = h0(r[n]);
		c[n] = h[n];
	}
}
double norm(double h_new[], double h_old[]) {
	double res = 0.0;
	double dr = R / N;
	for (int i = 0; i < N; i++) {
		res += fabs(h_new[i] - h_old[i]) * dr;
	}
	return res;
}
