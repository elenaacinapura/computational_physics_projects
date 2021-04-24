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
	Params_delta p;

	double E_start = -1.0;
	double dE = 0.001;
	double E_end = 0.0;

	FILE *file;
	/******************************************/
	/* Cosh potential */
	/******************************************/
	xi = 0.01;
	a = 2.0 / xi;

	L = 5.0;
	x0 = 0.3;
	A = 1.0;
	B = 1.0;

	p.L = L;
	p.dx = dx;
	p.x0 = x0;
	p.a = a;
	p.A = A;
	p.B = B;

	file = fopen("delta_cosh.csv", "w");
	assert(file != NULL);

	printf("------------------------------------------------\nDELTA METHOD - COSH POTENTIAL\nSearching bound energies for xi = %.3lf\n\n", xi);

	double E = E_start + dE;
	double delta, delta_old;
	int cnt_bound = 0;
	int cnt = 0;

	while (E < E_end) {
		double delta_new = Delta_E_cosh(E, &p);

		/* Print */
		fprint_double(file, E);
		fprint_double(file, delta_new);
		fprintf(file, "\n");

		if (cnt > 1) {
			if (delta_new * delta <= 0.0) {
				if ((delta_new > delta && delta > delta_old) || (delta_new < delta && delta < delta_old)) {
					double E_bound;
					E_bound = zero_bisection(Delta_E_cosh, E - dE, E, &p);
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

	a = 2.0 * m * s * s * eps / (pow(hbar, 2));
	xi = 2.0 / a;
	int l = 1;

	L = 3.0;
	x0 = 1.1;
	A = 1.0;
	B = 1.0;

	p.L = L;
	p.dx = dx;
	p.x0 = x0;
	p.a = a;
	p.A = A;
	p.B = B;
	p.l = l;

	file = fopen("delta_lj.csv", "w");
	assert(file != NULL);
	printf("------------------------------------------------\nDELTA METHOD - LENNARD JONES POTENTIAL\nSearching bound energies for l = %d\n\n", l);

	E = E_start + dE;
	cnt_bound = 0;
	cnt = 0;

	while (E < E_end) {
		double delta_new = Delta_E_lj(E, &p);

		/* Print */
		fprint_double(file, E);
		fprint_double(file, delta_new);
		fprintf(file, "\n");

		if (cnt > 1) {
			if (delta_new * delta <= 0.0) {
				if ((delta_new > delta && delta > delta_old) || (delta_new < delta && delta < delta_old)) {
						double E_bound = zero_bisection(Delta_E_lj, E - dE, E, &p);
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
}