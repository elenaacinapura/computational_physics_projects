#include <assert.h>
#include <complex.h>
#include <differential_eq/numerov.h>
#include <math.h>
#include <numerical_methods/zero_bisection.h>
#include <print_routines.h>

#include "../routines.h"

/*================ PARAMETERS ===============*/
const double xi = 0.05;
const double L = 5.0;
const double x0 = 0.3;
const double dx = 0.001;
/*================ STRUCTURES ===============*/
typedef struct Params_cosh {
	double a, E;
} Params_cosh; /* a = 2 m s^2 V0 / hbar^2 */
typedef struct Params_delta {
	double L, dx, x0, a, A, B;
} Params_delta;
/*================ FUNCTION HEADERS ==============*/
double F_cosh(double x, void *param);
double Delta_E_cosh(double E, void *param);

/*================ MAIN ===============*/
int main() {
    double a = 2.0 / xi;
	double A = 1.0;
	double B = 1.0;

    double E_start = -1.0;
	double dE = 0.0005;
	double E_end = 0.0;

	Params_delta p_cosh;
	p_cosh.L = L;
	p_cosh.dx = dx;
	p_cosh.x0 = x0;
	p_cosh.a = a;
	p_cosh.A = A;
	p_cosh.B = B;

	/*============ Welcome ============*/
    printf("=================================================\n");
	printf("BOUND STATES WITH NUMEROV ALGORITHM\n");
	printf("=================================================\n");
    printf("Parameters:\n");
	printf("\tV(x) = -V0*cosh(x/a)^(-4)\n");
	printf("\txi = hbar^2 / (m V0 a^2) = %.4lf\n", xi);
    printf("\tL = %.1lf\n", L);
    printf("\tdx = %.4lf\n", dx);
    printf("\tx0 = %.1lf\n", x0);
	printf("Starting calculating...\n\n");
	printf("Finding bound states:\n\n");

	/*============ Cycle on E ============*/
    FILE *file;
	file = fopen("cosh.csv", "w");
	assert(file != NULL);

	double E = E_start + dE;
	double delta, delta_old;
	int cnt_bound = 0;
	int cnt = 0;

	while (E < E_end) {
		/* Do numerov and calculate difference in derivatives */
		double delta_new = Delta_E_cosh(E, &p_cosh);

		/* Print to file */
		fprint_double(file, E);
		fprint_double(file, delta_new);
		fprintf(file, "\n");

		/* Check if delta changes sign and search exact energy */
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
		/* Advance */
		delta_old = delta;
		delta = delta_new;
		E += dE;
		cnt++;
	}
	fclose(file);

	printf("\nNumber of bound states found: %d\n", cnt_bound);
	printf("=================================================\n");
}
/*================ FUNCTIONS ==============*/
double F_cosh(double x, void *param) {
	Params_cosh *p = (Params_cosh *)param;
	double a = p->a;
	double E = p->E;
	return -a * (pow(cosh(x), -4) + E);
}
double Delta_E_cosh(double E, void *param) {
	/* Extract parameters */
	Params_delta *p = (Params_delta *)param;
	double L, dx, x0, a, A, B;
	L = p->L;
	dx = p->dx;
	x0 = p->x0;
	a = p->a;
	A = p->A;
	B = p->B;
	/* Create arrays */
	int dimF = (int)ceil((L + x0) / dx) + 1;
	int dimB = (int)ceil((L - x0) / dx) + 1;
	double xF[dimF], xB[dimB];
	double phiF[dimF], phiB[dimB];

	Params_cosh param_numerov;
	param_numerov.a = a;
	param_numerov.E = E;

	double k = sqrt(a * (-E));
	/* phiF */
	xF[0] = -L;
	xF[1] = -L + dx;
	phiF[0] = A * exp(k * (-L));
	phiF[1] = A * exp(k * (-L + dx));
	execute_numerov(xF, phiF, dx, dimF, F_cosh, &param_numerov);
	/* phiB */
	xB[0] = L;
	xB[1] = L - dx;
	phiB[0] = B * exp(-k * L);
	phiB[1] = B * exp(-k * (L - dx));
	execute_numerov(xB, phiB, -dx, dimB, F_cosh, &param_numerov);

	double R = phiF[dimF - 1] / phiB[dimB - 1];
	for (int i = 0; i < dimB; i++) {
		phiB[i] *= R;
	}

	double delta = calculate_delta(x0, phiF, phiB, dx, dimF, dimB, F_cosh, &param_numerov);

	return delta;
}