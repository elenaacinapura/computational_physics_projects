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
const double R = 6.0;
const double T = 0.5;
const double alpha = 0.01;
const double rho_target = 0.5;

/*============== FUNCTIONS HEADERS ==============*/
/**
 * Repulsive part of Lennard-Jones interaction
 */
double h0(double r);

/**
 * Initialize the vectors.
 */
double init(double r[], double h[], double c[]);

/**
 * Calculate the distance between two functions.
 */
double Dist(double h_new[], double h_old[]);

/*============== MAIN ==============*/
int main() {
	/* Welcome */
	printf("============================================================\n");
	printf("PERCUS - YEVICK SIMULATION\n");
	printf("============================================================\n");
	printf("Parameters of the simulation:\n");
	printf("\tN = %d\n\tR = %.1lf\n\tT = %.1lf\n\trho = %.1lf\n\talpha = %.2lf\n", N, R, T, rho_target, alpha);
	printf("\nCalculating...\n\n");
	/* My vectors */
	double r[N];
	double h[N];
	double c[N];
	double F[N];

	double rho_start = 0.1;
	double d_rho = 0.001;

	double norm;
	
	init(r, h, c);
	double rho = rho_start;
	while (rho < rho_target) {
		do {
			double c_t[N], h_t[N], c_new[N], h_new[N]; /* transforms and new vectors */
			vec_copy(N, c, c_t);
			vec_copy(N, h, h_t);

			fft_radial_forward(h_t, N, R);
			fft_radial_forward(c_t, N, R);

			for (int i = 0; i < N; i++) {
				F[i] = h_t[i] * c_t[i] * rho * 2.0 * M_PI * R / N;
				assert(!isnan(F[i]));
			}

			fft_radial_inverse(F, N, R);

			for (int i = 0; i < N; i++) {
				c_new[i] = (F[i] + 1.0) * h0(r[i]);
				h_new[i] = c_new[i] + F[i];

				c_new[i] = alpha * c_new[i] + (1.0 - alpha) * c[i];
				h_new[i] = alpha * h_new[i] + (1.0 - alpha) * h[i];
			}

			norm = Dist(h_new, h);
			// fprint_double_newline(stdout, norm);

			vec_copy(N, h_new, h);
			vec_copy(N, c_new, c);
		} while (norm > 1e-6);
		rho += d_rho;
	}
	printf("Calculations ended successfully!\n");
	printf("============================================================\n");

	/*======================= PRINT ========================*/
	FILE* file;
	file = fopen("percus.csv", "w");
	fprintf(file, "r\tg\n");

	for (int i = 0; i < N; i++) {
		fprint_double(file, r[i]);
		fprint_double_newline(file, h[i] + 1.0);
	}
	fclose(file);

	FILE *f_par;
	f_par = fopen("parameters.csv", "w");
	fprintf(f_par, "rho\tT\n");
	fprint_double(f_par, rho_target);
	fprint_double_newline(f_par, T);
	fclose(f_par);
}

/*============== FUNCTIONS ==============*/
double h0(double r) {
	if (r < 1e-5) {
		return -1.0;
	}
	double V = 4.0 * pow(r, -12);
	return exp(-V / T) - 1.0;
}
double init(double r[], double h[], double c[]) {
	for (int n = 0; n < N; n++) {
		r[n] = n * R / N;
		h[n] = h0(r[n]);
		c[n] = h[n];
	}
}
double Dist(double h_new[], double h_old[]) {
	double res = 0.0;
	double dr = R / N;
	for (int i = 0; i < N; i++) {
		res += fabs(h_new[i] - h_old[i]) * dr;
	}
	return res;
}
