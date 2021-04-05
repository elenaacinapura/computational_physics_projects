/********************************************************/
/* SINGLE SIMULATION TO PLOT THE SOLUTION */
/********************************************************/
#include <assert.h>
#include <complex.h>
#include <gnuplot_i.h>
#include <math.h>
#include <numerov.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

double F_square(double x, void *p) {
    /* Extract parameters */
    Param_V_square *param = (Param_V_square *) p;
    double V0 = param->V0;
    double a = param->a;

    double xi = 2.0*a*a*V0;

	x = (double)x;  /* just to be sure */

	if (fabs(x) <  0.5) {
		return xi * 2.0 * (1.0 - E);
	}
	return -xi *2.0 * E;
}

double V_square(double x, void *p) {
	if (fabs(x) < 1.0 / 2.0) {
		return 1.0;
	}
	return 0.0;
}

int main() {
	/* Boundary conditions */
	L = 30.0;
	E = 10.0;
	dx = 1e-3;
	double a = 1.0;
	double V0 = 1.0;

	FILE *outputfile;
	outputfile = fopen("tunneling_phi.csv", "w");

	execute_numerov(F_square, V_square, &(Param_V_square){V0, a}, outputfile, 1);

	fclose(outputfile);

	complex double T = calculate_T(phi.dim - 10, phi.dim - 30);
	complex double R = calculate_R(phi.dim - 10, phi.dim - 30);
	printf("|T|^2 = %lf\n", square_cabs(T));
	printf("|R|^2 = %lf\n", square_cabs(R));

	/* Plot */
	gnuplot_ctrl *h;
	h = gnuplot_init();
	gnuplot_set_term_png(h, "phi.png");
	gnuplot_cmd(h, "set xlabel \'x\'");
	char cmd[1000];
	sprintf(cmd, "set title \'E = %lf\'", E);
	gnuplot_cmd(h, cmd);
	gnuplot_cmd(h, "plot \'tunneling_phi.csv\' using 1:2 title \'V(x)\' w lines, \'tunneling_phi.csv\' using 1:3 title \'Re(phi(x))\' w lines");
	system("eog phi.png"); /* open png from terminal */
}