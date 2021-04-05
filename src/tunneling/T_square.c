/********************************************************/
/* CALCULATION OF T(E) */
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
	Param_V_square *param = (Param_V_square *)p;
	double V0 = param->V0;
	double a = param->a;

	double xi = 2.0 * a * a * V0;

	x = (double)x; /* just to be sure */

	if (fabs(x) < 0.5) {
		return xi * 2.0 * (1.0 - E);
	}
	return -xi * 2.0 * E;
}

double V_square(double x, void *p) {
	if (fabs(x) < 1.0 / 2.0) {
		return 1.0;
	}
	return 0.0;
}

int main(int argc, char **argv) {
	L = 50.0;
	dx = 1e-3;
	double a = 1.0;
	double V0 = 1.0;

	FILE *outputfile;
	outputfile = fopen("tunneling_T.csv", "w");
	/*
	assert(argc > 1);
	outputfile = fopen(argv[1], "w");
	*/
	E = 0.0;
	double dE = 0.01;
	while (E <= 3.0) {
		
		execute_numerov(F_square, V_square, &(Param_V_square){V0, a}, outputfile, 0);

		double T = calculate_T(phi.dim - 2, phi.dim - 10);
		fprint_double(outputfile, E);
		fprint_double(outputfile, square_cabs(T));
		fprintf(outputfile, "\n");

		E += dE;
	}
	fclose(outputfile);

	/* Plot */
	gnuplot_ctrl *h;
	h = gnuplot_init();
	gnuplot_set_term_png(h, "T.png");
	gnuplot_cmd(h, "set xlabel \'E\'");
	gnuplot_set_ylabel(h, "|T|^2");
	gnuplot_cmd(h, "set title \'T(E)\'");
	gnuplot_cmd(h, "plot \'tunneling_T.csv\' using 1:2 w lines");

	system("eog T.png");
}