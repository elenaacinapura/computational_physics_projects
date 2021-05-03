#include <assert.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>
#include <math.h>
#include <print_routines.h>

double func(double x) {
	return 2.0*x;
}

int main() {
	int N = 512;
	double L = 10.0;
	double dx = L / N;
	double dk = 2.0 * M_PI / L;

    FILE *par;
    par = fopen("parameters.csv", "w");
    fprintf(par, "N\tL\n");
    fprintf(par, "%d\t%lf\n", N, L);
    fclose(par);

	double x[N], f[2 * N];
    FILE *file;
    file = fopen("func.csv", "w");
    fprintf(file, "x\tf\n");

	for (int n = 0; n < N; n++) {
        x[n] = dx * n;
        f[2*n] = func(x[n]);
        f[2*n + 1] = 0.0;
        fprint_double(file, x[n]);
        fprint_double_newline(file, f[2*n]);
	}
    fclose(file);

    gsl_fft_complex_radix2_forward(f, 1, N);

    for (int k = 0; k < N; k++) {
        double K = dk * (k <= N / 2 ? k : k - N);
        complex double f_complex = f[2*k] + I * f[2*k + 1];
        f_complex *= I*K;
        f[2*k] = creal(f_complex);
        f[2*k + 1] = cimag(f_complex);
    }

    gsl_fft_complex_radix2_inverse(f, 1, N);
    
    file = fopen("derivative.csv", "w");
    fprintf(file, "x\tf\n");
    for (int n = 0; n < N; n++) {
        fprint_double(file, x[n]);
        fprint_double_newline(file, f[2*n]);
    }
    fclose(file);
}
