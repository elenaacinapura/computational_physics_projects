#include <assert.h>
#include <complex.h>
#include <differential_eq/numerov.h>
#include <math.h>
#include <numerical_methods/zero_bisection.h>
#include <print_routines.h>

#include "../routines.h"

/*================ PARAMETERS ===============*/
const double xi = 0.015;
const double LEFT_LIM = -4.1;
const double RIGHT_LIM = 0.0;
const double dx = 0.001;
const double x0 = -0.6; /* meeting point */
const double LEFT_COEFF = 2.0;
const double RIGHT_COEFF = 0.0;
/*================ STRUCTURES ===============*/
typedef struct Params {
	double a, E;
} Params;
/*================ FUNCTION HEADERS ==============*/
double Potential(double x);
double F(double x, void *param);
double Delta_E(double E, void *param);
void print_eigenfunction(double E, void *param);
/*================ MAIN ===============*/
int main() {
    int dimF = (int)ceil((x0 - LEFT_LIM) / dx) + 1;
	int dimB = (int)ceil((RIGHT_LIM - x0) / dx) + 1;
    int N = dimF + dimB -1;
	/*============ Welcome ============*/
	printf("=================================================\n");
	printf("FIND THE GROUND STATE WITH NUMEROV ALGORITHM\n");
	printf("=================================================\n");
	printf("Parameters:\n");
	printf("\txi = hbar^2 / (m V0 a^2) = %.3lf\n", xi);
	printf("\tInterval: [%.1lf, %.1lf]\n", LEFT_LIM, RIGHT_LIM);
	printf("\tdx = %.4lf\n", dx);
	printf("\tMeeting point x0 = %.1lf\n", x0);
    printf("\tNumber of points N = %d\n", N);
	printf("\nStarting calculating...\n\n");

	double a = 1.0 / xi;
	Params p;
	p.a = a;

	double E_start = -0.2;
	double E_end = 0.0;
	double dE = 1e-4;

	double delta, delta_old, E_ground;
	E_ground = -1.0;

	double E = E_start;
	int cnt = 0;
	while (E < E_end) {
		p.E = E;

		/* Do numerov and calculate difference in derivatives */
		double delta_new = Delta_E(E, &p);

		/* Check if delta changes sign and search exact energy */
		if (cnt > 1) {
			if (delta_new * delta <= 0.0) {
				if ((delta_new > delta && delta > delta_old) || (delta_new < delta && delta < delta_old)) {
					E_ground = zero_bisection(Delta_E, E - dE, E, &p);
					printf("E0 = %lf\n", E_ground);
					break;
				}
			}
		}
		/* Advance */
		delta_old = delta;
		delta = delta_new;
		E += dE;
		cnt++;

		E += dE;
	}
	if (E_ground > -0.9) {
		printf("\nPrinting eigenfunction to file.\n");
		print_eigenfunction(E, &p);
	}
}
/*================ FUNCTIONS ==============*/
double Potential(double x) {
	return -x * x * (x + 1);
}
double F(double x, void *param) {
	Params *p = (Params *)param;
	double a = p->a;
	double E = p->E;
	double V = Potential(x);
	return a * (V - E);
}
double Delta_E(double E, void *param) {
	Params *p = (Params *)param;
	p->E = E;

	/* Create arrays: forward and backward*/
	int dimF = (int)ceil((x0 - LEFT_LIM) / dx) + 1;
	int dimB = (int)ceil((RIGHT_LIM - x0) / dx) + 1;
	double xF[dimF], xB[dimB];
	double phiF[dimF], phiB[dimB];

	/* phiF */
	xF[0] = LEFT_LIM;
	xF[1] = LEFT_LIM + dx;
	phiF[0] = 0.0;
	phiF[1] = -1e-3;
	// fprint_double_newline(stdout, phiF[1]);
	execute_numerov(xF, phiF, dx, dimF, F, p);
	/* phiB */
	xB[0] = RIGHT_LIM;
	xB[1] = RIGHT_LIM - dx;
	phiB[0] = 0.0;
	phiB[1] = -1.0;
	execute_numerov(xB, phiB, -dx, dimB, F, p);

	double R = phiF[dimF - 1] / phiB[dimB - 1];
	assert(isfinite(R));
	for (int i = 0; i < dimB; i++) {
		phiB[i] *= R;
	}

	double delta = calculate_delta(x0, phiF, phiB, dx, dimF, dimB, F, p);

	return delta;
}

void print_eigenfunction(double E, void *param) {
	Params *p = (Params *)param;
	p->E = E;

	/* Create arrays: forward and backward*/
	int dimF = (int)ceil((RIGHT_LIM - LEFT_LIM + x0) / dx) + 1;
	int dimB = (int)ceil((RIGHT_LIM - LEFT_LIM - x0) / dx) + 1;
	double xF[dimF], xB[dimB];
	double phiF[dimF], phiB[dimB];

	/* phiF */
	xF[0] = LEFT_LIM;
	xF[1] = LEFT_LIM + dx;
	phiF[0] = 0.0;
	phiF[1] = 2.0;
	execute_numerov(xF, phiF, dx, dimF, F, p);
	/* phiB */
	xB[0] = RIGHT_LIM;
	xB[1] = RIGHT_LIM - dx;
	phiB[0] = 0.0;
	phiB[1] = -1.0;
	execute_numerov(xB, phiB, -dx, dimB, F, p);

	double R = phiF[dimF - 1] / phiB[dimB - 1];
	for (int i = 0; i < dimB; i++) {
		phiB[i] *= R;
	}

	/* normalize  */
	double N = 0.0;
	for (int i = 0; i < dimF; i++) {
		N += phiF[i];
	}
	for (int i = 1; i < dimB; i++) {
		N += phiB[i];
	}
	for (int i = 0; i < dimF; i++) {
		phiF[i] /= N;
	}
	for (int i = 1; i < dimB; i++) {
		phiB[i] /= N;
	}

	/* print to file */
	FILE *file;
	file = fopen("eigenfunction.csv", "w");
	fprintf(file, "x\tpsi\n");
	for (int i = 0; i < dimF; i++) {
		fprint_double(file, xF[i]);
		fprint_double_newline(file, phiF[i]);
	}
	for (int i = dimB - 1; i > 0; i--) {
		fprint_double(file, xB[i]);
		fprint_double_newline(file, phiB[i]);
	}
	fclose(file);
}