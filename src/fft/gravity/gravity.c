#include <assert.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>
#include <math.h>
#include <print_routines.h>
#define N 256

const double g = 9.81;

double omega(double k, double h);

complex double f0(double x, double L);

double group_velocity(double k_max, double h);

/***********************************************
 * MAIN
 ***********************************************/

int main() {
	/* Parameters */
	double T = 10;
	double dt = 0.02;
	double h = 10.0;
	double L = 30.0;
	double dx = L / N;
	double dk = 2.0 * M_PI / L;

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

	double t = 0.0;
	int start = 1;
	double f_max, v;
	int k_max;

	while (t < T) {
		gsl_fft_complex_radix2_forward(f, 1, N);

		if (start) {
			start = 0;
			f_max = cabs(f[0] + I * f[1]);
			for (int i = 0; i < N; i++) {
				int k = dk * (i <= N / 2 ? i : i - N);
				complex double f_curr = f[2 * i] + I * f[2 * i + 1];
				if (cabs(f_curr) > f_max) {
					f_max = cabs(f_curr);
					k_max = k;
				}
			}
			v = group_velocity(k_max * dk, h);
		}

		for (int i = 0; i < N; i++) {
			int k = dk * (i <= N / 2 ? i : i - N);
			double om = omega(k, h);
			/* Real and imag parts of e^iwt */
			double c = cos(om * t);
			double s = sin(-om * t);

			double fre = f[2 * i];
			double fim = f[2 * i + 1];

			/* Complex multiplication */
			f[2 * i] = fre * c - fim * s;
			f[2 * i + 1] = fre * s + fim * c;
		}

		gsl_fft_complex_radix2_inverse(f, 1, N);

		for (int n = 0; n < N; n++) {
			fprint_double(file, x[n] -v*t);
			fprint_double_newline(file, f[2 * n]);
		}

		t += dt;
	}

	fclose(file);
	printf("k_max_index = %d\nk_max = %lf\ngroup velocity = %lf\n", k_max, k_max * dk, v);
}

/***********************************************
 * END OF MAIN
 ***********************************************/

double omega(double k, double h) {
	return sqrt(k * tanh(k * h));
}

complex double f0(double x, double L) {
	double x0 = L * 0.5;
	return exp(-(x - x0) * (x - x0) / 8.0) * cexp(I * 2 * M_PI * (x - x0));
}

double group_velocity(double k_max, double h) {
	double eps = 1e-4;
	return 0.5 / eps * (omega(k_max + eps, h) - omega(k_max - eps, h));
}