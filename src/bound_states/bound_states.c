#include <assert.h>
#include <complex.h>
#include <math.h>
#include <differential_eq/numerov.h>
#include <print_routines.h>
#include <numerical_methods/zero_bisection.h>

#include "util.h"

int main() {
	/* Environment variables */
	double L = 5.0;
	double dx = 1e-3;
	double x0 = 0.15;
	double A = -1.0;
	double B = -2.0;
	double xi = 0.005;
	double a = 2.0 / xi;
	Params_delta p;
	p.L = L; p.dx = dx; p.x0 = x0; p.a = a; p.A = A; p.B = B;

	double E_start = -1.0;
	double dE = 0.0001;
	double E_end = 0.0;

	FILE *file;
	file = fopen("delta.csv", "w");
	assert(file != NULL);

	double E = E_start+dE;
	double delta;
	int first_time = 1;
	int cnt_bound = 0;
	while (E < E_end) {
		double delta_new = Delta_E_cosh(E, &p);
		/* PRINT */
		fprint_double(file, E);
		fprint_double(file, delta_new);
		fprintf(file, "\n");
		if (!first_time){
			if (delta_new * delta <= 0.0) {
				double E_bound = zero_bisection(Delta_E_cosh, E - dE, E, &p);
				printf("New bound state found with E = %lf\n", E_bound);
				cnt_bound++;
			}
		}
		delta = delta_new;
		E += dE;
		first_time = 0;
	}
	printf("\nNumber of bound states found: %d\n", cnt_bound);
	fclose(file);

	system("gnuplot delta.gp");
}