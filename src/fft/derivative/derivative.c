#include <assert.h>
#include <complex.h>
#include <gsl/gsl_fft_complex.h>
#include <math.h>
#include <print_routines.h>

double func(double x) {
	return sin(x);
}

int main() {
    /* Parameters */
	int N = 512;
	double L = 50.0;
	double dx = L / N;
	double dk = 2.0 * M_PI / L;

    /* Print parameters to file */
    FILE *par;
    par = fopen("parameters.csv", "w");
    fprintf(par, "N\tL\n");
    fprintf(par, "%d\t%lf\n", N, L);
    fclose(par);

	double x[N], f[2 * N];

    FILE *file;
    file = fopen("func.csv", "w");
    fprintf(file, "x\tf\n");

    /* Fill f(x) and print it to file */
	for (int n = 0; n < N; n++) {
        x[n] = dx * n;
        f[2*n] = func(x[n]);
        f[2*n + 1] = 0.0;
        fprint_double(file, x[n]);
        fprint_double_newline(file, f[2*n]);
	}
    fclose(file);

    /* Transform */
    gsl_fft_complex_radix2_forward(f, 1, N);

    /* Derivative */
    for (int k = 0; k < N; k++) {
        complex double K =  I * dk * (k <= N / 2 ? k : k - N);
        complex double f_complex = f[2*k] + I * f[2*k + 1];
        f_complex *= K;
        f[2*k] = creal(f_complex);
        f[2*k + 1] = cimag(f_complex);
    }

    /* Antitransform */
    gsl_fft_complex_radix2_inverse(f, 1, N);
    
    /* Print derivative to file */
    file = fopen("derivative.csv", "w");
    fprintf(file, "x\tf\n");
    for (int n = 0; n < N; n++) {
        fprint_double(file, x[n]);
        fprint_double_newline(file, f[2*n]);
    }
    fclose(file);
}
