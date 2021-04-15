#include <assert.h>
#include <complex.h>
#include <math.h>
#include <numerov.h>
#include <print_routines.h>

#include "util.h"

int main() {
	/* Environment variables */
	double L = 5.0;
	double dx = 1e-3;
	double x0 = 0.3;
	double A = 1.0;
	double B = 2.0;

	double xi = 0.05;
	double a = 2.0 / xi;

	int dimF = (int)ceil((L + x0) / dx) + 1;
	int dimB = (int)ceil((L - x0) / dx) + 1;
	double xF[dimF], xB[dimB];
	double phiF[dimF], phiB[dimB];

	double E_start = -1.0;
	double dE = 0.001;
	double E_end = 0.0;

	FILE *file;
	file = fopen("delta.csv", "w");
	assert(file != NULL);

	double E = E_start+dE;
	while (E < E_end) {
		double k = sqrt(a * (-E));
		Params p = {a, E};

		/* phiF */
		xF[0] = -L;
		xF[1] = -L + dx;
		phiF[0] = A * exp(k * (-L));
		phiF[1] = A * exp(k * (-L + dx));
		execute_numerov(xF, phiF, dx, dimF, F_cosh, &p);

		/* phiB */
		xB[0] = L;
		xB[1] = L - dx;
		phiB[0] = B * exp(-k * L);
		phiB[1] = B * exp(-k * (L - dx));
		execute_numerov(xB, phiB, -dx, dimB, F_cosh, &p);

		// printf("%lf\t%lf\n", phiF[dimF-1], phiB[dimB-1]);

		double R = phiF[dimF - 1] / phiB[dimB - 1];


		for (int i = 0; i < dimB; i++) {
			phiB[i] *= R;
		}

		double delta = calculate_delta(x0, phiF, phiB, dx, dimF, dimB, F_cosh, &p);
		fprint_double(file, E);
		fprint_double(file, delta);
		fprintf(file, "\n");

		E += dE;
	}
	fclose(file);

	system("gnuplot delta.gp -p");
}