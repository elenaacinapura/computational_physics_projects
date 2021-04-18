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
	int potential_type = 1;

	/* Parameters */
	double kB = 1.38e-23;
	double eps = 68.5 * kB;
	double s = 3.18e-10;
	double hbar = 1.05e-34;
	double m = 1.66e-27;

	double xi, E_start, E_end, dE, E;
	if (potential_type == 0) {
		xi = 1.0;
		dE = 0.01;
		E_start = 0.001;
		E_end = 10.0;
	} else if (potential_type = 1){
		xi = (2.0*m*s*s*eps/(pow(hbar, 2))); /* xi = (2 m a^2 V0) / h^2 */
		E_start = 0.05*1.6e-22 / eps;
		E_end = 5.0*1.6e-22 / eps;
		dE = 0.05*1.6e-22 / eps;

	}

	/* How far I will go */
	double L = 50.0;
	double dx = 0.005;
	int dim = (int)(L / dx);
	/* Arrays for position ang wave function */
	double x[dim];
	double u[dim];

	/* Preparing for the loop */
	E = E_start;
	double sigma;
	FILE *file;
	file = fopen("sigma_tot.csv", "w");
	assert(file != NULL);

	while (E <= E_end) {	/* Loop on energies */
		double k = sqrt(xi * E);
		sigma = 0.0;

		for (int l = 0; l < 20; l++) {	/* Loop on l */
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

			/* Execute numerov */
			Params p = {xi, E, l};
			if (potential_type == 0) {
				execute_numerov(x, u, dim, dx, F_hard, &p);
			} else if (potential_type == 1) {
				execute_numerov(x, u, dim, dx, F_lj, &p);
			}

			/* Calculate delta_l and sigma_tot */
			int idx1 = dim - 1;
			int idx2 = dim - 20;
			double r1 = x[idx1];
			double r2 = x[idx2];
			double u_r1 = u[idx1];
			double u_r2 = u[idx2];

			double delta_l = find_delta(r1, r2, u_r1, u_r2, l, k);

			sigma += 4 * M_PI / (k * k) * (double)(2 * l + 1) * sin(delta_l) * sin(delta_l);
		}
		if (potential_type == 1) {
			fprint_double(file, E / 1.6e-22 * eps);
		}
		if (potential_type == 0) {
			fprint_double(file, E);
		}
		fprint_double(file, sigma * s*s / (1e-20));
		fprintf(file, "\n");

		E += dE;

	}
	fclose(file);

	system("gnuplot sigma_tot.gp -p");
}