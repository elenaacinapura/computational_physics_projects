#include <assert.h>
#include <linear_algebra/array_routines.h>
#include <linear_algebra/blas_wrappers.h>
#include <math.h>
#include <numerical_methods/fft_radial.h>
#include <physics/lennard_jones.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>

/*======================= CONSTANTS ========================*/
const int N = 4096;
const double R = 8.0;
const double rho_target = 0.2;
const double T = 1.35;
const double alpha = 0.01;			 /* adjusting parameter */
const double alpha_potential = 20.0; /* parameter in the potential */

/*======================= FUNCTION HEADERS ========================*/
double potential(double r);
double deriv_potential(double r);
double h0(double r);
double init(double r[], double h[], double c[]);
double Dist(double h_new[], double h_old[]);
double pressure(double r[], double g[]);
double int_energy(double r[], double g[]);

/*======================= MAIN ========================*/
int main() {
	/*======================= WELCOME ========================*/
	printf("============================================================\n");
	printf("THE HYPERNETTED CHAIN APPROXIMATION\n");
	printf("============================================================\n");
	printf("Parameters of the simulation:\n");
	printf("\tN = %d\n\tR = %.1lf\n\tT = %.2lf\n\tRho target = %.1lf\n\tAlpha = %.2lf (aka \"gentleness\")\n", N, R, T, rho_target, alpha);
	printf("============================================================\n");

	printf("Calculating...\n");
	/*======================= PERCUS-YEVICK METHOD ========================*/
	/* My vectors */
	double r[N];
	double h[N];
	double c[N];
	double F[N];

	double rho_start = 0.05;
	double d_rho = 0.005;

	double norm;
	int stuck_in_loop = 0;

	init(r, h, c);
	double rho = rho_start;
	printf("Current value of rho = %.3lf", rho);
	while (rho <= rho_target + 0.1 * d_rho) { /* some times rho=rho_target fails just for binary representation*/
		printf("\rCurrent value of rho = %.3lf", rho);
		fflush(stdout);
		do {
			double c_t[N], h_t[N], c_new[N], h_new[N]; /* transforms and new vectors */
			/*========= Convolution =========*/
			vec_copy(N, c, c_t);
			vec_copy(N, h, h_t);

			fft_radial_forward(h_t, N, R);
			fft_radial_forward(c_t, N, R);

			for (int i = 0; i < N; i++) {
				F[i] = h_t[i] * c_t[i] * rho * 2.0 * M_PI * R / N;
				assert(!isnan(F[i]));
			}

			fft_radial_inverse(F, N, R);

			/*========= Update h and c =========*/
			for (int i = 0; i < N; i++) {
				c_new[i] = h0(r[i]) * exp(F[i]) + exp(F[i]) - 1.0 - F[i];
				h_new[i] = c_new[i] + F[i];

				c_new[i] = alpha * c_new[i] + (1.0 - alpha) * c[i];
				h_new[i] = alpha * h_new[i] + (1.0 - alpha) * h[i];
			}
			/*========= Calculate norm =========*/
			double norm_old = norm;
			norm = Dist(h_new, h);
			if (norm_old == norm) {
				stuck_in_loop++;
			}
			if (stuck_in_loop > 5) {
				printf("Calculation failed, stuck in a loop. Decrease the value of alpha.\n");
				printf("============================================================\n");
				exit(0);
			}
			vec_copy(N, h_new, h);
			vec_copy(N, c_new, c);
		} while (norm > 1e-6);
		rho += d_rho;
	}
	printf("\n============================================================\n");
	/*======================= PRESSURE AND INTERNAL ENERGY ========================*/
	double g[N];
	double P, U;
	for (int i = 0; i < N; i++) {
		g[i] = h[i] + 1.0;
	}
	P = pressure(r, g);
	U = int_energy(r, g);
	printf("Results of the calculations:\n");
	printf("\tP = %.3lf\n\tU per particles = %.3lf\n", P, U);
	printf("============================================================\n");

	/*======================= PRINT ========================*/
	FILE* file;
	file = fopen("percus.csv", "w");
	fprintf(file, "r\tg\n");

	for (int i = 0; i < N; i++) {
		fprint_double(file, r[i]);
		fprint_double_newline(file, g[i]);
	}
	fclose(file);

	FILE* f_par;
	f_par = fopen("parameters.csv", "w");
	fprintf(f_par, "rho\tT\n");
	fprint_double(f_par, rho_target);
	fprint_double_newline(f_par, T);
	fclose(f_par);
}

/*============== FUNCTIONS ==============*/
/**
 * Potential.
 */
double potential(double r) {
	if (r <= 0.5) {
		r = 0.5;
	}
	return 1.0 / (1.0 - 6.0 / alpha_potential) * (6.0 / alpha_potential * exp(alpha_potential - alpha_potential * r) - pow(r, -6));
}
/**
 * Potential.
 */
double deriv_potential(double r) {
	if (r <= 0.5) {
		return 0;
	}
	return 1.0 / (1.0 - 6.0 / alpha_potential) * (-6.0 * exp(alpha_potential - alpha_potential * r) + 6.0 * pow(r, -7));
}
/**
 * Initial value for h at low densities
 */
double h0(double r) {
	double V = potential(r);
	return exp(-V / T) - 1.0;
}
/**
 * Initialize the vectors.
 */
double init(double r[], double h[], double c[]) {
	for (int n = 0; n < N; n++) {
		r[n] = n * R / N;
		h[n] = h0(r[n]);
		c[n] = h[n];
	}
}
/**
 * Calculate the distance between two functions (L^2)
 */
double Dist(double h_new[], double h_old[]) {
	double res = 0.0;
	double dr = R / N;
	for (int i = 0; i < N; i++) {
		res += fabs(h_new[i] - h_old[i]) * dr;
	}
	return res;
}
/**
 * Calculate the pressure of the fluid.
 */
double pressure(double r[], double g[]) {
	double res = 0.0;
	double dr = R / N;
	for (int n = 1; n < N; n++) {
		res += pow(r[n], 3) * g[n] * deriv_potential(r[n]) * dr;
	}
	res *= -2.0 * M_PI * pow(rho_target, 2) / 3.0;
	return res;
}
/**
 * Calculate the average internal energy per particle of the fluid.
 */
double int_energy(double r[], double g[]) {
	double res = 0.0;
	double dr = R / N;
	for (int i = 1; i < N; i++) {
		res += potential(r[i]) * g[i] * pow(r[i], 2) * dr;
	}
	res *= rho_target * 2.0 * M_PI;
	return res;
}
