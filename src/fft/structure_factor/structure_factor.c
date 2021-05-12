#include <assert.h>
#include <math.h>
#include <numerical_methods/fft_radial.h>
#include <print_routines.h>

/*============== CONSTANTS ==============*/
const int N = 1024;
const double R = 50.0;
const double T = 0.8;
const double rho = 0.5;

/*============== FUNCTIONS ==============*/
double h(double r) {
	if (fabs(r) < 1e-5) {
		return -1.0;
	}
	double V = 4.0 * (pow(r, -12) - pow(r, -6));
	return exp(-V/T) - 1.0;
}

/*============== MAIN ==============*/
int main() {
	double dr = R / N;
	double dk = M_PI / R;
	double r[N];
	double K[N];
	double k;
	double S[N];
	double S_integral[N];


	for (int n = 0; n < N; n++) {
		r[n] = (double)n * dr;
		K[n] = n * dk;
		S[n] = h(r[n]);

	}
	/* Calculate S(K) via transform */
	fft_radial_forward(S, N, R);

	/* Calculate S(K) via integral */
	for (int k = 0; k < N; k++) {
		double sum = 0.0;
		for (int n = 0; n < N; n++) {
			sum += 2.0 * h(r[n]) * r[n] / K[k] * sin(K[k] * r[n]);
		}
		S_integral[k] = sum;
	}

	FILE *file;
	file = fopen("structure.csv", "w");
	assert(file != NULL);

	fprintf(file, "r\tS\tS_int\n");

	for (int n = 0; n < N; n++) {
		S[n] = S[n] * rho + 1.0;
		S_integral[n] = S_integral[n] * rho + 1.0;
		fprint_double(file, K[n]);
		fprint_double(file, S[n]);
		fprint_double_newline(file, S_integral[n]);
	}
	fclose(file);
}