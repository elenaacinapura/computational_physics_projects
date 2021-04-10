#include <assert.h>
#include <complex.h>
#include <gnuplot_i.h>
#include <math.h>
#include <numerov.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

int main() {
	double L = 20.0;
	double dx = 1e-3;
	double xi_true;		/* xi_true = h*h/2*m*|V_0|*a*a */
	Param_F param;

	int dim = (int)(2.0 * L / dx);
	double x[dim];
	complex double phi[dim];

	E = 0.0;
	double dE = 0.01;
	double k;
	complex double T;

    FILE *file;
    file = fopen("T.csv", "w");

	while (E <= 5.0) {

		fprint_double(file,E);
		
		/* primo valore */
		xi_true = 0.2;
		param.xi = 1/xi_true;

        k = sqrt(param.xi) * sqrt(E);
		x[0] = -L;
		x[1] = -L + dx;
		phi[0] = cexp(-I * k * (-L));
		phi[1] = cexp(-I * k * (-L + dx));

        solve_numerov(x, phi, dim, dx, F_cosh, &param, 0, file);

        T = T_coeffcient(x, phi, dim, k);

        fprint_double(file, T*conj(T));
		fprint_double(file, 1-T*conj(T));

		/* secondo valore */
		xi_true = 0.05;
		param.xi = 1/xi_true;

        k = sqrt(param.xi) * sqrt(E);
		x[0] = -L;
		x[1] = -L + dx;
		phi[0] = cexp(-I * k * (-L));
		phi[1] = cexp(-I * k * (-L + dx));

        solve_numerov(x, phi, dim, dx, F_cosh, &param, 0, file);

        T = T_coeffcient(x, phi, dim, k);

        fprint_double(file, T*conj(T));
		fprint_double(file, 1-T*conj(T));

		/* terzo valore */
		xi_true = 0.01;
		param.xi = 1/xi_true;

        k = sqrt(param.xi) * sqrt(E);
		x[0] = -L;
		x[1] = -L + dx;
		phi[0] = cexp(-I * k * (-L));
		phi[1] = cexp(-I * k * (-L + dx));

        solve_numerov(x, phi, dim, dx, F_cosh, &param, 0, file);

        T = T_coeffcient(x, phi, dim, k);

        fprint_double(file, T*conj(T));
		fprint_double(file, 1-T*conj(T));

        fprintf(file, "\n");
		E += dE;
	}

    fclose(file);


    /* Plot */
    gnuplot_ctrl *h;
    h = gnuplot_init();
	gnuplot_cmd(h, "load \'plot_TR.gp\'");
	system("eog T.png");
	system("eog R.png");

}
