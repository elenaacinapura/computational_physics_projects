#include <assert.h>
#include <complex.h>
#include <differential_eq/numerov.h>
#include <math.h>
#include <print_routines.h>

#include "../routines.h"

/*================= PARAMETERS ================*/
int N = 1000;
double L = 1.0;
double xi = 1.0;

/*================= STRUCTURES ================*/
typedef struct Params {
	double alfa, E;
} Params;

/*================= FUNCTION HEADERS ================*/
double F(double x, void *param);

/*================= MAIN ================*/
int main() {
	/*================= WELCOME ================*/
	printf("======================================================\n");
	printf("BOUND STATES OF A PERIODIC POTENTIAL USING NUMEROV\n");
	printf("======================================================\n");
	printf("Parameters:\n");
	printf("\tL = %.1lf\n\tN = %d\n\txi = hbar^2 / (2 m L^2 V0) = %.1lf\n", L, N, xi);
	printf("Calculating...\n");
	printf("Status:\t%.0lf%%", 0.0);
	fflush(stdout);
	

	double dx = L / (N - 1);
	double alfa = 1.5;
	Params p;
	p.alfa = alfa;

	double x[N];
	double phi1[N];
	double phi2[N];
	double complex psi[N];

	double K_start = -M_PI;
	double K_end = M_PI;
	double dK = 0.5;
	int num_K = ceil(2.0*K_end/dK);
	double percentage_step = 100.0/num_K;

	double dE = 0.001;
	double E_start = -1.0;

	FILE *file;
	file = fopen("periodic.csv", "w");
	fprintf(file, "K\tE0\tE1\tE2\n");

	/* Cycle on K */
	double K = K_start;
	int cnt_K = 0;
	while (K <= K_end) {
		fprint_double(file, K);

		/* Cycle on E */
		int cnt_bound = 0;
		double E = E_start;
		while (cnt_bound <= 2) {
			p.E = E;

			/* Initialize x and phis */
			x[0] = 0.0;
			x[1] = dx;
			phi1[0] = 0.0;
			phi1[1] = 1.0 / L * dx;
			phi2[0] = 1.0;
			phi2[1] = 1.0;

			execute_numerov(x, phi1, dx, N, F, &p);
			execute_numerov(x, phi2, dx, N, F, &p);

			double complex beta = phi1[N - 1] / (cexp(I * K * L) - phi2[N - 1]);

			for (int i = 0; i < N; i++) {
				psi[i] = phi1[i] + beta * phi2[i];
			}

			double Delta = cabs((cexp(I * K * L) * psi[1] + psi[N - 2] - psi[N - 1] * (2.0 + dx * dx * F(x[N - 1], &p))) / dx);


			if (fabs(Delta) < 1e-3) {
				cnt_bound++;
				if (cnt_bound <= 2) {
					fprint_double(file, E);
				} else {
					fprint_double_newline(file, E);
				}
				E += 1.0;
				continue;
			}
			E += dE;
		}
		K += dK;
		
		cnt_K++;
		printf("\rStatus:\t%.0lf%%", percentage_step*cnt_K);
		fflush(stdout);
	}
	fclose(file);
	printf("\nCalculations ended successfully!\n");
	printf("======================================================\n");
}

/*================= FUNCTIONS ================*/
double F(double x, void *param) {
	Params *p = (Params *)param;
	double alfa = p->alfa;
	double E = p->E;
	if (x < 0.3) {
		return alfa * (1.0 - E);
	}
	return -alfa * E;
}