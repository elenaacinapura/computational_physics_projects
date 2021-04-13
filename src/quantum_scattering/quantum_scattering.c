#include <assert.h>
#include <complex.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

#include "util.h"

int main() {
	double L = 50.0;
	double dx = 0.001;
	int dim = (int)(L / dx);
	double xi = 1.0; /* xi = (2 m a^2 V0) / h^2 */
    double E = 1.0;
    double k = sqrt(xi*E);

	double x[dim];
	double u[dim];

	double sigma = 0.0;

	for (int l = 0; l < 10; l++) {
		/* Set initial conditions for hard sphere*/
		double x0 = 1.0;
		double u0 = 1.0;
		x[0] = x0;
		x[1] = x0 + dx;
		u[0] = 0.0;
		u[0] = u0;

		Params p = {xi, 1, E};

		execute_numerov(x, u, dim, dx, F_hard, &p);

		int idx1 = dim -1;
        int idx2 = dim - 10;
        double r1 = x[idx1];
        double r2 = x[idx2];
        double u1 = u[idx1];
        double u2 = u[idx2];

        double delta = find_delta(r1, r2, u1, u2, l, k);

        sigma += 4*M_PI /(k*k) * (double)(2*l+1) * sin(delta)*sin(delta);
	}
    printf("Sigma = %lf\n", sigma);
}