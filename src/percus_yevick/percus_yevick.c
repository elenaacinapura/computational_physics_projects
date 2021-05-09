#include <assert.h>
#include <math.h>
#include <numerical_methods/fft_radial.h>
#include <physics/lennard_jones.h>
#include <print_routines.h>
#include <stdio.h>

/*============== CONSTANTS ==============*/
const int N = 512;
const double R = 50.0;
const double rho = 0.7;
const double alpha = 0.5;

/*============== FUNCTIONS HEADERS ==============*/
double h0(double r);
double init(double h[], double c[]);

/*============== MAIN ==============*/
int main() {
	double r[N];
	double h[N];
	double c[N];
	double F[N];

	init(h, c);

	for (int cnt = 0; cnt < 100; cnt++) {
		fft_radial_forward(h, N, R);
		fft_radial_forward(c, N, R);

		for (int i = 0; i < N; i++) {
			F[i] = h[i] * c[i] * rho * 2.0 * M_PI * R / N;
		}

        fft_radial_inverse(F, N, R);
        fft_radial_inverse(h, N, R);
        fft_radial_inverse(c, N, R);

        for (int i = 0; i < N; i++) {
            double r = i * R/N;
            double c_old = c[i]; 
            double h_old = h[i];

            c[i] = (F[i] + 1.0) * h0(r);
            h[i] = c[i] + F[i];

            c[i] = alpha * c[i] + (1.0 - alpha) * c_old;
            h[i] = alpha * h[i] + (1.0 - alpha) * h_old;
        }
	}

    FILE * file;
    file = fopen("percus.csv", "w");
    fprintf(file, "r\tg\n");

    for (int i = 0; i < N; i++) {
        fprint_double(file, i * R/N);
        fprint_double_newline(file, h[i] + 1.0);
    }
    fclose (file);
}

/*============== FUNCTIONS ==============*/
double h0(double r) {
	if (r < 1e-5) {
		return -1.0;
	}
	struct Empty {
	} s;
	double V = lj_potential(r, &s);
	return exp(-V) - 1.0;
}
double init(double h[], double c[]) {
	for (int n = 0; n < N; n++) {
		double r = n * R / N;
		h[n] = h0(r);
		c[n] = h[n];
	}
}