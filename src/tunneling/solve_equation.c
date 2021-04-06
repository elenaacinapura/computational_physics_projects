#include <gnuplot_i.h>
#include <numerov.h>
#include <print_routines.h>
#include "util.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main () {
    double L = 20.0;
    double dx = 1e-3;

    /* Parameters that I choose */
    E = 0.5;
    Param_F_square param = {1.0};

    int dim = (int)(2.0*L/dx);
    double k = param.sqrt_xi*sqrt(E);

    double x [dim];
    complex double phi [dim];

    x[0] = -L;
    x[1] = -L + dx;
    phi[0] = cexp(-I*k*(-L));
    phi[1] = cexp(-I*k*(-L + dx));

    FILE *file;
    file = fopen("solution.csv", "w");
    solve_numerov(x, phi, dim, dx, F_square, &param, 1, file);
    fclose(file);

    complex double T = T_coeffcient(x, phi, dim, k);
    printf("|T|^2 = %lf\n", T*conj(T));

    /* Plot */
    gnuplot_ctrl *h;
    h = gnuplot_init();
    gnuplot_set_term_png(h, "solution.png");
	gnuplot_cmd(h, "set xlabel \'x\'");
	char cmd[1000];
	sprintf(cmd, "set title \'E = %lf, xi = %lf\'", E, param.sqrt_xi*param.sqrt_xi);
	gnuplot_cmd(h, cmd);
	gnuplot_cmd(h, "plot \'solution.csv\' using 1:2 title \'Re(phi(x))\' w lines");
	system("eog solution.png"); /* open png from terminal */
}