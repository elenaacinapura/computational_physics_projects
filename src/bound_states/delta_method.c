#include <assert.h>
#include <complex.h>
#include <differential_eq/numerov.h>
#include <math.h>
#include <numerical_methods/zero_bisection.h>
#include <print_routines.h>

#include "util.h"

int main() {
	double dx = 1e-3;
	double a, xi;
	double L, x0, A, B;

	double E_start = -1.0;
	double dE = 0.0005;
	double E_end = 0.0;

	FILE *file;
	/******************************************/
	/* Cosh potential */
	/******************************************/
	xi = 0.05;
	a = 2.0 / xi;

	L = 5.0;
	x0 = 0.3;
	A = 1.0;
	B = 1.0;

	Params_delta p_cosh;
	p_cosh.L = L;
	p_cosh.dx = dx;
	p_cosh.x0 = x0;
	p_cosh.a = a;
	p_cosh.A = A;
	p_cosh.B = B;

	file = fopen("delta_cosh.csv", "w");
	assert(file != NULL);

	printf("------------------------------------------------\nDELTA METHOD - COSH POTENTIAL\nSearching bound energies for xi = %.3lf\n\n", xi);

	double E = E_start + dE;
	double delta, delta_old;
	int cnt_bound = 0;
	int cnt = 0;

	while (E < E_end) {
		double delta_new = Delta_E_cosh(E, &p_cosh);

		/* Print */
		fprint_double(file, E);
		fprint_double(file, delta_new);
		fprintf(file, "\n");

		if (cnt > 1) {
			if (delta_new * delta <= 0.0) {
				if ((delta_new > delta && delta > delta_old) || (delta_new < delta && delta < delta_old)) {
					double E_bound;
					E_bound = zero_bisection(Delta_E_cosh, E - dE, E, &p_cosh);
					printf("E%d = %lf\n", cnt_bound, E_bound);
					cnt_bound++;
				}
			}
		}
		delta_old = delta;
		delta = delta_new;
		E += dE;
		cnt++;
	}

	printf("\nNumber of bound states found: %d\n", cnt_bound);
	printf("------------------------------------------------\n");

	fclose(file);
	/******************************************/
	/* Lennard-Jones potential */
	/******************************************/
	/* Parameters */
	double kB = 1.38e-23;
	double eps = 35.7 * kB;
	double s = 2.79e-10;
	double hbar = 1.05e-34;
	double m = 10 * 1.66e-27;

	a = 2.0 * m * s * s * eps / pow(hbar, 2);
	xi = 2.0 / a;

	Params_delta p_lj;

	dx = 1e-3;
	L = 3.0;
	x0 = 1.1;
	A = -3.0;
	B = 2.0;

	p_lj.L = L;
	p_lj.dx = dx;
	p_lj.x0 = x0;
	p_lj.a = a;
	p_lj.A = A;
	p_lj.B = B;

	file = fopen("delta_lj.csv", "w");
	assert(file != NULL);
	printf("------------------------------------------------\nDELTA METHOD - LENNARD JONES POTENTIAL\n\n");
	cnt_bound = 0;

	for (int l = 0; l < 10; l++) {
		printf("l = %d\n", l);
		p_lj.l = l;

		E = E_start + dE;
		cnt = 0;

		while (E < E_end) {
			double delta_new = Delta_E_lj(E, &p_lj);

			/* Print */
			fprint_double(file, E);
			fprint_double(file, delta_new);
			fprintf(file, "\n");

			if (cnt > 1) {
				if (delta_new * delta <= 0.0) {
					if ((delta_new > delta && delta > delta_old) || (delta_new < delta && delta < delta_old)) {
						double E_bound = zero_bisection(Delta_E_lj, E - dE, E, &p_lj);
						printf("E%d = %lf\n", cnt_bound, E_bound);
						cnt_bound++;
					}
				}
			}
			delta_old = delta;
			delta = delta_new;
			E += dE;
			cnt++;
		}
	}
	printf("\nNumber of bound states found: %d\n", cnt_bound);
	printf("------------------------------------------------\n");
	fclose(file);
}