#include <assert.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>
#include <math.h>
#include <print_routines.h>

/*======================= CONSTANTS ========================*/
const int N = 4096;
const double L = 400.0;
const int num_T = 500;
const double E = 1.0;
double alpha = 0.025;

/*======================= FUNCTIONS HEADERS ========================*/
complex double f0(double x, double k);
double V(double x);

/*======================= MAIN ========================*/
int main() {
	/*======================= Welcome ========================*/
	printf("\n=======================================================\n");
	printf("Quantum particle on a 1D barrier\n");
	printf("=======================================================\n");
	printf("Parameters:\n");
	printf("\t N = %d \n\t L = %.1lf \n\t Num_timesteps = %d \n\t E = %.1lf \n\t alpha = h_bar^2 / (2mV0a^2) = %lf\n", N, L, num_T, E, alpha);
	
	/*=========================================================*/
	double dx = L / N;
	double dk = 2.0 * M_PI / L;

	double x[N];
	double complex rho[N];
	double complex xi[N];
	double phi[2 * N];

	double dt = 1.0 / sqrt(E);
	double k_E = sqrt(E / alpha);

	double norm_start, norm_R, norm_T;
	norm_start = 0.0;
	norm_R = 0.0;
	norm_T = 0.0;
	printf("\t dt = %lf\n", dt);
	printf("\nCalculating...\n\n");
	/*======================= Initialize arrays ========================*/
	for (int n = 0; n < N; n++) {
		/* x */
		x[n] = dx * n;
		/* rho */
		rho[n] = cexp(-I * V(x[n]) * dt);
		/* xi */
		double k = dk * (double)(n <= N / 2 ? n : n - N);
		complex double Delta_k_2 = -k * k;
		xi[n] = cexp(I * alpha * dt * Delta_k_2);
		/* wave packet: gaussian moving to the right */
		phi[2 * n] = creal(f0(x[n], k_E));
		phi[2 * n + 1] = cimag(f0(x[n], k_E));
		/* norm of phi */
		complex double phi_c = phi[2 * n] + I * phi[2 * n + 1];
		norm_start += phi_c * conj(phi_c);
	}

	/* ============================= Evolution ===============================*/
	FILE *file;
	file = fopen("schrodinger.csv", "w");
	assert(file != NULL);
	fprintf(file, "x\tf\n");

	int cnt = -1;
	double t = 0.0 - dt;
	while (cnt < num_T) {
		printf("\rTimestep = %d", cnt);
		t += dt;
		cnt++;
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
		/* Transform */
		gsl_fft_complex_radix2_forward(phi, 1, N);
		/* Multiply by xi[n] */
		for (int n = 0; n < N; n++) {
			double complex phi_c = phi[2 * n] + I * phi[2 * n + 1];
			phi_c *= xi[n];
			phi[2 * n] = creal(phi_c);
			phi[2 * n + 1] = cimag(phi_c);
		}
		/* Antitransform */
		gsl_fft_complex_radix2_inverse(phi, 1, N);
		if (fabs(phi[2*N - 2]) > 1e-2) {
			break;
		}
	}
	fclose(file);
	/* ============================= End of evolution ===============================*/
	/* Calculate norm_R and norm_T */
	for (int n = 0; n < N; n++) {
		complex double phi_c = phi[2 * n] + I * phi[2 * n + 1];
		if (x[n] < L * 0.5 - 1.0) {
			norm_R += phi_c * conj(phi_c);
		} else if (x[n] > L * 0.5 + 1) {
			norm_T += phi_c * conj(phi_c);
		}
	}
	double T = norm_T / norm_start;
	printf("Simulation ended successfully!\n");
	printf("\nResults of the simulation:\n");
	printf("\t T = %lf\n", T);
	printf("=======================================================\n");


	/* Print parameters */
	FILE *f_par;
	f_par = fopen("parameters.csv", "w");
	assert(f_par != NULL);
	fprintf(f_par, "N\tT\tL\tE\n%d\t%d\t%lf\t%lf\n", N, num_T, L, E);
	fclose(f_par);
}

/*======================= FUNCTIONS ========================*/
complex double f0(double x, double k) {
	double sigma = 40.0 / k;
	double x0 = L / 2.0 - 5.0 * sigma;
	assert(x0 > 0.0);
	return exp(-(x - x0) * (x - x0) / pow(sigma, 2)) * cexp(I * k * (x - x0));
}
double V(double x) {
	if (fabs(x - L / 2) < 1.0) {
		return 1.0;
	}
	return 0.0;
}