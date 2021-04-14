#include <complex.h>
#include <math.h>
#include <numerov.h>
#include <print_routines.h>

#include "util.h"

int main() {
	double L = 5.0;
	double dx = 1e-3;
	double x0 = 0.3;
	double E = 1.0;
	double xi = 1.0;
	double a = 2.0 / xi;
	double k = sqrt(a * E);
	double A = 1.0;
    double B = 2.0;

	Params p = {a, E};

	int dimF = (int)ceil((L + x0) / dx)+1;
	int dimB = (int)ceil((L - x0) / dx)+1;
	double xF[dimF], xB[dimB];
	double phiF[dimF], phiB[dimB];

	/* phiF */
	xF[0] = -L;
	xF[1] = -L + dx;
	phiF[0] = A * exp(k * (-L));
	phiF[1] = A * exp(k * (-L + dx));
    execute_numerov(xF, phiF, dx, dimF, F_cosh, &p);

    printf("%lf\n", xF[dimF-1]);

    /* phiB */
	xB[0] = L;
	xB[1] = L - dx;
	phiB[0] = B * exp(- k * L);
	phiB[1] = B * exp(-k * (L - dx));
    execute_numerov(xB, phiB, -dx, dimB, F_cosh, &p);

    double R = phiF[dimF-1] / phiB[dimB-1];

    for (int i = 0; i < dimB; i++) {
        phiB[i] *= R;
        printf("%lf\t%lf\n", phiF[i], phiB[i]);
    }

    double delta = calculate_delta(x0, phiF, phiB, dx, dimF, dimB, F_cosh, &p);

    printf("%lf\n", xB[dimB-1]);

    fprint_double(stdout, delta);
}