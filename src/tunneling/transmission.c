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
	Param_F_square param = {1.0};

	int dim = (int)(2.0 * L / dx);
	double x[dim];
	complex double phi[dim];

	E = 0.0;
	double dE = 0.01;
	double k;
	complex double T;

    FILE *file;
    file = fopen("T.csv", "w");

	while (E <= 4.0) {

        k = param.sqrt_xi * sqrt(E);
		x[0] = -L;
		x[1] = -L + dx;
		phi[0] = cexp(-I * k * (-L));
		phi[1] = cexp(-I * k * (-L + dx));

        solve_numerov(x, phi, dim, dx, F_square, &param, 0, file);

        T = T_coeffcient(x, phi, dim, k);

        fprint_double(file, E);
        fprint_double(file, T*conj(T));
        fprintf(file, "\n");

		E += dE;
	}

    fclose(file);

    /* Plot */
    gnuplot_ctrl *h;
    h = gnuplot_init();
    gnuplot_set_term_png(h, "T.png");
	gnuplot_set_xlabel(h, "E");
    gnuplot_set_ylabel(h, "|T|^2");
	gnuplot_cmd(h, "set title \'T(E)\'");
    gnuplot_cmd(h, "unset key");
	gnuplot_cmd(h, "plot \'T.csv\' using 1:2 w lines");
	system("eog T.png"); /* open png from terminal */
}
