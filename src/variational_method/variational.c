#include <assert.h>
#include <gnuplot_i.h>
#include <integrate_notypesafe.h>
#include <math.h>
#include <minimum.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

int main() {
	double E;
	Params p;
	double xi[] = {0.025, 0.005, 0.0025};
	p.xi = 0.01;
	double integral_cutoff = 50.0;
	double dalpha = 0.01;

	FILE *file;
	file = fopen("E.csv", "w");

	double alpha = 0.08;
	while (alpha <= 1.2) {
		p.alpha = alpha;
		fprint_double(file, alpha);

		for (int i = 0; i < 3; i++) {
			p.xi = xi[i];
			E = integrate(integrand, -integral_cutoff, integral_cutoff, 1000, &p);

			fprint_double(file, E);
		}
		fprintf(file, "\n");

		alpha += dalpha;
	}
	fclose(file);

    double alpha_low = 0.1; 
    double alpha_high = 2.0;

	for (int i = 0; i < 3; i++) {
        p.xi = xi[i];
        double a0 = minimum(energy, alpha_low, alpha_high, &p);
        p.alpha = a0;
        printf("\nxi = %lf\t a0 = %lf\tE0 = %lf\n", xi[i], a0, energy(a0, &p));
	}
    printf("\n");

	system("gnuplot plot_E.gp -p");

	// gnuplot_ctrl *h;
	// h = gnuplot_init();
	// gnuplot_cmd(h, "load \'plot_E.gp\'");
}