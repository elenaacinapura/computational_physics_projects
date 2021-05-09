#include <assert.h>
#include <math.h>
#include <numerical_methods/fft_radial.h>
#include <print_routines.h>

const int N = 256;
const double R = 50.0;
const double rho = 0.1;

double h(double r) {
	if (fabs(r) < 1e-5) {
		return -1.0;
	}
	struct Empty {
	} s;
	double V = 4.0 * (pow(r, -12) - pow(r, -6));
	return exp(-V) - 1.0;
}

int main() {
	double dr = R / N;
	double r[N];
	double S[2 * N];

	for (int n = 0; n < N; n++) {
		r[n] = (double)n * dr;

		if (isnan(h(r[n]))) {
			printf("nan value found in h(r) for n = %d ", n);
		}
		S[2 * n] = h(r[n]);
		S[2 * n + 1] = 0.0;
	}
	fft_radial_forward(S, N, R);

	FILE *file;
	file = fopen("structure.csv", "w");
	assert(file != NULL);

	fprintf(file, "r\tS\n");

	for (int n = 0; n < N; n++) {
		double k = M_PI / R * n;
		S[2 * n] = S[2 * n] * rho + 1.0;

		fprint_double(file, k);
		fprint_double_newline(file, S[2 * n]);
	}
	fclose(file);
}