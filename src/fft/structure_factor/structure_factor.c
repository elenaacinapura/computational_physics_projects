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
	double r[N];
	double S[N];

	for (int n = 0; n < N; n++) {
		r[n] = (double)n * dr;
		S[n] = h(r[n]);

	}
	fft_radial_forward(S, N, R);

	FILE *file;
	file = fopen("structure.csv", "w");
	assert(file != NULL);

	fprintf(file, "r\tS\n");

	for (int n = 0; n < N; n++) {
		double k = M_PI / R * n;
		S[n] = S[n] * rho + 1.0;

		fprint_double(file, k);
		fprint_double_newline(file, S[n]);
	}
	fclose(file);
}