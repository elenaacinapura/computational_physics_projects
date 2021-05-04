#include <assert.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>
#include <math.h>
#include <print_routines.h>

#define MOVING_FRAME 0

const int N = 256;

double omega(double k, double h);

complex double f0(double x, double L);

double group_velocity(double k_max, double h);

/***********************************************
 * MAIN
 ***********************************************/

int main() {
	/* Parameters */
	double T = 1000.0;
	double dt = 0.5;
	double h = 10.0;
	double L = 50.0;
	double dx = L / N;
	double dk = 2.0 * M_PI / L;

	/* Welcome */
	printf("\n***************************************************\n");
	printf("A wave package in a gravitational field\n\nCalculating...\n\n");

	FILE *f_par;
	f_par = fopen("parameters.csv", "w");
	assert(f_par != NULL);
	fprintf(f_par, "N\tT\tdt\tL\n%d\t%lf\t%lf\t%lf\n", N, T, dt, L);
	fclose(f_par);

	/* Structures */
	double f[2 * N];
	double x[N];

	FILE *file;
	file = fopen("f.csv", "w");
	assert(f != NULL);
	fprintf(file, "x\tf\n");

	for (int i = 0; i < N; i++) {
		x[i] = i * dx;
		f[2 * i] = creal(f0(x[i], L));
		f[2 * i + 1] = cimag(f0(x[i], L));
	}

	gsl_fft_complex_radix2_forward(f, 1, N);

	/* Find k_max */
	double k_max;
	double f_max = cabs(f[0] + I * f[1]);
	for (int i = 0; i < N; i++) {
		double k = dk * (i <= N / 2 ? i : i - N);
		complex double f_curr = f[2 * i] + I * f[2 * i + 1];
		if (cabs(f_curr) > f_max) {
			f_max = cabs(f_curr);
			k_max = k;
		}
	}
	double v = group_velocity(k_max, h);

	double t = 0.0;
	while (t < T) {
		/* Evolution */
		for (int i = 0; i < N; i++) {
			int k = dk * (i <= N / 2 ? i : i - N);
			double om = omega(k, h);
		
			double fre = f[2 * i];
			double fim = f[2 * i + 1];
			complex double f_c = fre + I * fim;
			f_c *= cexp(-I * om * dt);
			if (MOVING_FRAME) {
				f_c *= cexp(I * k * v*dt);
			}

			f[2 * i] = creal(f_c);
			f[2 * i + 1] = cimag(f_c);
		}

		gsl_fft_complex_radix2_inverse(f, 1, N);

		for (int n = 0; n < N; n++) {
			fprint_double(file, x[n]);
			fprint_double_newline(file, f[2 * n]);
		}

		gsl_fft_complex_radix2_forward(f, 1, N);

		t += dt;
	}

	fclose(file);
	printf("Simulation ended successfully.\n\n");
	printf("Results:\nk_max = %lf\nGroup velocity = %lf\n", k_max, v);
	printf("\n***************************************************\n");
}

/***********************************************
 * END OF MAIN
 ***********************************************/

double omega(double k, double h) {
	return sqrt(k * tanh(k * h));
}

complex double f0(double x, double L) {
	double x0 = L * 0.3;
	return exp(-(x - x0) * (x - x0) / 8.0) * cexp(I * 2 * M_PI * (x - x0));
}

double group_velocity(double k_max, double h) {
	double eps = 1e-4;
	return 0.5 / eps * (omega(k_max + eps, h) - omega(k_max - eps, h));
}