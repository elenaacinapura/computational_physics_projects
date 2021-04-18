#include "util.h"

#include <complex.h>
#include <differential_eq/numerov.h>
#include <math.h>
#include <print_routines.h>
#include <assert.h>

double F_cosh(double x, void *param) {
	Params_cosh *p = (Params_cosh *)param;
	double a = p->a;
	double E = p->E;
	return -a * (pow(cosh(x), -4) + E);
}

double F_lj(double r, void *param) {
    Params_lj *p = (Params_lj *)param;
	double a = p->a;
    double E = p->E;
    int l = p->l;
    assert(r != 0.0);
    return a*4.0*(pow(r, -12) - pow(r, -6)) + (double)(l*(l+1))/(r*r) - a*E;
}

void execute_numerov(double x[], double phi[], double dx, int dim, double F(double, void *), void *p) {
	int i = 2;
	do {
		phi[i] = numerov_1D(x[i - 1], phi[i - 1], phi[i - 2], dx, F, p);
		x[i] = x[i - 1] + dx;
		i++;
	} while (i < dim);
}

double calculate_delta(double x0, double phiF[], double phiB[], double dx, int dimF, int dimB, double F(double, void *), void *p) {
	return 1.0 / dx * (phiF[dimF - 2] + phiB[dimB - 2] - phiF[dimB - 1] * (2.0 + dx * F(x0, p)));
}

double Delta_E_cosh(double E, void *param) {
	Params_delta *p = (Params_delta *)param;
	double L, dx, x0, a, A, B;
	L = p->L;
	dx = p->dx;
	x0 = p->x0;
	a = p->a;
	A = p->a;
	B = p->B;
	/* Create arrays */
	int dimF = (int)ceil((L + x0) / dx) + 1;
	int dimB = (int)ceil((L - x0) / dx) + 1;
	double xF[dimF], xB[dimB];
	double phiF[dimF], phiB[dimB];

	double k = sqrt(a * (-E));
	Params_cosh param_numerov = {a, E};

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

    double delta  = calculate_delta(x0, phiF, phiB, dx, dimF, dimB, F_cosh, &p);

    return delta;
}

double Delta_E_lj(double E, void *param) {
	Params_delta *p = (Params_delta *)param;
	double L, dx, x0, a, A, B;
	L = p->L;
	dx = p->dx;
	x0 = p->x0;
	a = p->a;
	A = p->a;
	B = p->B;
    int l = p->l;
	/* Create arrays */
	int dimF = (int)ceil(x0 / dx) + 1;
	int dimB = (int)ceil((L - x0) / dx) + 1;
	double xF[dimF], xB[dimB];
	double phiF[dimF], phiB[dimB];

	double k = sqrt(a * (-E));
	Params_lj param_numerov;
    param_numerov.a = a;
    param_numerov.E = E;
    param_numerov.l = l;

	/* phiF */
    double r0 = 1.0;
	xF[0] = r0;
	xF[1] = r0 + dx;
	phiF[0] = A * exp(-sqrt(4.0*a/25.0) * pow(r0, -5));
	phiF[1] = A * exp(-sqrt(4.0*a/25.0) * pow(r0 + dx, -5));
	execute_numerov(xF, phiF, dx, dimF, F_lj, &param_numerov);
    // for (int i = 0; i < dimF; i++) {
    //     printf("%lf\t%lf\n", xF[i], phiF[i]);
    // }
	/* phiB */
	xB[0] = L;
	xB[1] = L - dx;
	phiB[0] = B * exp(-k * L);
	phiB[1] = B * exp(-k * (L - dx));
	execute_numerov(xB, phiB, -dx, dimB, F_lj, &param_numerov);
    // for (int i = 0; i < dimF; i++) {
    //     printf("%lf\t%lf\n", xB[i], phiB[i]);
    // }

	double R = phiF[dimF - 1] / phiB[dimB - 1];
	for (int i = 0; i < dimB; i++) {
		phiB[i] *= R;
	}

    double delta  = calculate_delta(x0, phiF, phiB, dx, dimF, dimB, F_lj, &p);

    return delta;
}