#include <assert.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>
#include <math.h>
#include <print_routines.h>

const int N = 256;
const double L = 10.0;
const int T = 1000;

double V(double x) {
	if (fabs(x - L / 2) < 1.0) {
		return 1.0;
	}
	return 0.0;
}

int main() {
	double dx = L / N;
	double dk = 2.0 * M_PI / L;
	double dt = 5e-2;

	double alpha = 0.025;
	double E = 1.0;
	double k_E = sqrt(E / alpha);

	double x[N];
	double complex rho[N];
	double complex xi[N];
	double phi[2 * N];

	/* Welcome */
	printf("\n***************************************************\n");
	printf("Quantum particle on a 1D barrier\n\nCalculating...\n\n");

	/* Print parameters */
	FILE *f_par;
	f_par = fopen("parameters.csv", "w");
	assert(f_par != NULL);
	fprintf(f_par, "N\tT\tdt\tL\n%d\t%d\t%lf\t%lf\n", N, T, dt, L);
	fclose(f_par);

	for (int n = 0; n < N; n++) {
		x[n] = dx * n;

		rho[n] = cexp(-I * V(x[n]) * dt);

		double k = dk * (n <= N / 2 ? n : n - N);
		complex double Delta_k_2 = -k * k;
		xi[n] = cexp(I * alpha * dt * Delta_k_2);

		/* Free particle */
		phi[2 * n] = cos(k_E * x[n]);
		phi[2 * n + 1] = sin(I * k_E * x[n]);
	}

	FILE *file;
	file = fopen("schrodinger.csv", "w");
	assert(file != NULL);
	fprintf(file, "x\tf\n");

	for (int t = 0; t < T; t++) {
		/* Print */
		for (int n = 0; n < N; n++) {
			fprint_double(file, x[n]);
			fprint_double_newline(file, phi[2 * n]);
		}
		/* Multiply by rho[n] */
		for (int n = 0; n < N; n++) {
			double complex phi_c = phi[2 * n] + I * phi[2 * n + 1];
			phi_c *= rho[n];
			phi[2 * n] = creal(phi_c);
			phi[2 * n + 1] = cimag(phi_c);
		}
		gsl_fft_complex_radix2_forward(phi, 1, N);
		for (int n = 0; n < N; n++) {
			double complex phi_c = phi[2 * n] + I * phi[2 * n + 1];
			phi_c *= xi[n];
			phi[2 * n] = creal(phi_c);
			phi[2 * n + 1] = cimag(phi_c);
		}
		gsl_fft_complex_radix2_inverse(phi, 1, N);
	}
	fclose(file);

	printf("Simulation ended successfully.\n");
	printf("\n***************************************************\n");
}