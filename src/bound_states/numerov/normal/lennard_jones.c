#include <assert.h>
#include <complex.h>
#include <differential_eq/numerov.h>
#include <math.h>
#include <numerical_methods/zero_bisection.h>
#include <print_routines.h>

#include "../routines.h"

/*================ CONSTANTS ===============*/
const double kB = 1.38e-23;
const double eps = 35.7 * kB;
const double s = 2.79e-10;
const double hbar = 1.05e-34;
const double m = 10 * 1.66e-27;

/*================ STRUCTURES ===============*/
typedef struct Params_lj {
    double a, E;
    int l;
} Params_lj;
typedef struct Params_delta {
    double L, dx, x0, r_start, a, A, B;
    int l;
} Params_delta;

/*================ FUNCTION HEADERS ==============*/
double F_lj(double r, void *param);
double Delta_E_lj(double E, void *param);

/*================ MAIN ===============*/
int main() {
	/*==================== PARAMETERS ==================*/
	double L = 3.0;
	double dx = 1e-3;
	double r0 = 1.1;
	double r_start = 0.7;
	double A = -3.0;
	double B = 2.0;

	double a = 2.0 * m * s * s * eps / pow(hbar, 2);

	Params_delta p_lj;
	p_lj.L = L;
	p_lj.dx = dx;
	p_lj.x0 = r0;
	p_lj.r_start = r_start;
	p_lj.a = a;
	p_lj.A = A;
	p_lj.B = B;

	/*============ Welcome ============*/
	printf("=================================================\n");
	printf("BOUND STATES WITH NUMEROV ALGORITHM\n");
	printf("=================================================\n");
	printf("Parameters:\n");
	printf("\tV(x) = Lennard-Jones for Neon 10\n");
	printf("\tL = %.1lf\n", L);
	printf("\tdx = %.4lf\n", dx);
	printf("\tx0 = %.1lf\n", r0);
	printf("\tr_start = %.1lf\n", r_start);
	printf("Starting calculating...\n\n");
	printf("Finding bound states:\n\n");

	/*============ Find boun states ============*/
	FILE *file;
	file = fopen("lj.csv", "w");
	assert(file != NULL);
	int cnt_bound = 0;

	/*============ Cycle on the angular momentum ============*/
	for (int l = 0; l < 10; l++) {
		printf("l = %d\n", l);
		p_lj.l = l;

		/*============ Cycle on E ============*/
		double E_start = -1.0;
		double dE = 0.0005;
		double E_end = 0.0;
		double E = E_start + dE;
		int cnt = 0;
		double delta, delta_old;
		while (E < E_end) {
			/* Do numerov and calculate difference in derivatives */
			double delta_new = Delta_E_lj(E, &p_lj);

			/* Print to file */
			fprint_double(file, E);
			fprint_double(file, delta_new);
			fprintf(file, "\n");

			/* Check if delta changes sign and search exact energy */
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
	printf("=================================================\n");
	fclose(file);
}
/*================ FUNCTIONS ==============*/
double F_lj(double r, void *param) {
    Params_lj *p = (Params_lj *)param;
	double a = p->a;
    double E = p->E;
    int l = p->l;
    assert(r != 0.0);

    return a*4.0*(pow(r, -12) - pow(r, -6)) + (double)(l*(l+1))/(r*r) - a*E;
}
double Delta_E_lj(double E, void *param) {
	Params_delta *p = (Params_delta *)param;
	double L, dx, x0, a, A, B, r_start;
	L = p->L;
	dx = p->dx;
	x0 = p->x0;
	r_start = p->r_start;
	a = p->a;
	A = p->A;
	B = p->B;
    int l = p->l;

	/* Create arrays */
	int dimF = (int)ceil((x0 - r_start)/ dx) + 1;	/* x0 is the meeting point */
	int dimB = (int)ceil((L - x0) / dx) + 1;
	double xF[dimF], xB[dimB];
	double phiF[dimF], phiB[dimB];

	double k = sqrt(a * (-E));
	Params_lj param_numerov;
    param_numerov.a = a;
    param_numerov.E = E;
    param_numerov.l = l;

	/* phiF */
	xF[0] = r_start;
	xF[1] = r_start + dx;
	phiF[0] = A * exp(-sqrt(4.0*a/25.0) * pow(r_start, -5));
	phiF[1] = A * exp(-sqrt(4.0*a/25.0) * pow(r_start + dx, -5));
	execute_numerov(xF, phiF, dx, dimF, F_lj, &param_numerov);

	/* phiB */
	xB[0] = L;
	xB[1] = L - dx;
	phiB[0] = B * exp(-k * L);
	phiB[1] = B * exp(-k * (L - dx));
	execute_numerov(xB, phiB, -dx, dimB, F_lj, &param_numerov);

	double R = phiF[dimF - 1] / phiB[dimB - 1];
	for (int i = 0; i < dimB; i++) {
		phiB[i] *= R;
	}

    double delta  = calculate_delta(x0, phiF, phiB, dx, dimF, dimB, F_lj, &param_numerov);

    return delta;
}
