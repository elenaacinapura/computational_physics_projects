#include <assert.h>
#include <complex.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

int main() {
	/* 	Potential type: 
		0 for hard sphere
		1 for lennard-jones
	*/
	int potential_type = 0;
	/* Parameters */
	double xi = 1.0; /* xi = (2 m a^2 V0) / h^2 */
	double E;
	double E_start = 0.0;
	double E_end = 10.0;
	double dE = 0.5;
	/* How far I will go */
	double L = 50.0;
	double dx = 0.001;
	int dim = (int)(L / dx);

	double x[dim];
	double u[dim];

	E = E_start;
	double sigma;
	FILE *file;
	file = fopen("sigma_tot.csv", "w");
	assert(file != NULL);

	while (E <= E_end) {
		E += dE;
		double k = sqrt(xi * E);
		sigma = 0.0;

		for (int l = 0; l < 20; l++) {
			/* Set initial conditions*/
			double x0, u0, u1;
			if (potential_type == 0) { /* Hard sphere */
				x0 = 1.0;
				u0 = 0.0;
				u1 = 1.0;
			} else if (potential_type == 1) { /* Lennard-Jones */
				x0 = 0.5;
				u0 = exp(-sqrt(4.0 * xi / 25.0) * pow(x0, -5));
				u1 = exp(-sqrt(4.0 * xi / 25.0) * pow(x0 + dx, -5));
			}
			x[0] = x0;
			u[0] = u0;
			x[1] = x0 + dx;
			u[1] = u1;

			Params p = {xi, E, l};

			execute_numerov(x, u, dim, dx, F_hard, &p);

			int idx1 = dim - 1;
			int idx2 = dim - 20;
			double r1 = x[idx1];
			double r2 = x[idx2];
			double u_r1 = u[idx1];
			double u_r2 = u[idx2];

			double delta_l = find_delta(r1, r2, u_r1, u_r2, l, k);

			sigma += 4 * M_PI / (k * k) * (double)(2 * l + 1) * sin(delta_l) * sin(delta_l);
		}
		fprint_double(file, E);
		fprint_double(file, sigma);
		fprintf(file, "\n");
	}
	fclose(file);

	system("gnuplot sigma_tot.gp -p");
}