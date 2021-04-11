#include <gnuplot_i.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>
#include <verlet.h>

#include "util.h"

int main() {
	double dt = 1e-3;
	double x0 = -10.0;
	double z0 = 15.0;
	double cutoff = z0;
	double dz = 0.1;	/* width of the bins of pi(z) */
	double dz_start = 0.005;	/* z step for shooting particles */
	int dim = 2 * (int)rint(cutoff / dz);

	double pi[6][dim];
	double pos[2], v[2], a[2];
	double E[] = {0.1, 0.3, 0.5, 1.0, 2.0, 5.0};
	int E_dim = 6;

	int cnt;

	/*------------------------------------------------*/
	/* theta = 0 */
	/*------------------------------------------------*/

	double theta = 0.0;

	for (int e = 0; e < E_dim; e++) { /* cycle on E */
		/* set E */
		Params par = {E[e], theta};

		/* empty pi(z) */
		for (int i = 0; i < dim; i++) {
			pi[e][i] = 0.0;
		}
		/* loop on z from z0 to z0 */
		double z = -z0;
		cnt = 0;
		while (z <= z0) {
			pos[0] = x0; /* initial x */
			pos[1] = z;	 /* initial z */
			v[0] = sqrt(2.0 * E[e]);	/* initial v_x */
			v[1] = 0.0;					/* initial v_z */

			execute_verlet(pos, v, a, dt, &par);	/* self-explanatory*/

			int idx = (int)rint(pos[1] / dz) + dim / 2;		/* position of final z in p(z) */
			if (idx >= 0 && idx < dim) {					/* restict to z in the range (-z0, z0) */
				pi[e][idx] += 1.0;
			}

			z += dz_start;
			cnt++;
		}
	}

	FILE *file;
	file = fopen("pi.csv", "w");

	for (int i = 0; i < dim; i++) {
		fprint_double(file, -cutoff + dz * i);
		for (int e = 0; e < E_dim; e++) {
			pi[e][i] /= cnt;
			fprint_double(file, pi[e][i]);
		}
		fprintf(file, "\n");
	}
	fclose(file);

	/*------------------------------------------------*/
	/* theta = 0 */
	/*------------------------------------------------*/

	theta = -M_PI_4;

	for (int e = 0; e < E_dim; e++) { /* cycle on E */
		/* set E */
		Params par = {E[e], theta};

		/* empty pi(z) */
		for (int i = 0; i < dim; i++) {
			pi[e][i] = 0.0;
		}
		/* loop on z from z0 to z0 */
		double z = -z0;
		cnt = 0;
		while (z <= z0) {
			pos[0] = x0; /* initial x */
			pos[1] = z;	 /* initial z */
			v[0] = sqrt(2.0 * E[e]);	/* initial v_x */
			v[1] = 0.0;					/* initial v_z */

			execute_verlet(pos, v, a, dt, &par);	/* self-explanatory*/

			int idx = (int)rint(pos[1] / dz) + dim / 2;		/* position of final z in p(z) */
			if (idx >= 0 && idx < dim) {					/* restict to z in the range (-z0, z0) */
				pi[e][idx] += 1.0;
			}

			z += dz_start;
			cnt++;
		}
	}

	file = fopen("pi_45.csv", "w");

	for (int i = 0; i < dim; i++) {
		fprint_double(file, -cutoff + dz * i);
		for (int e = 0; e < E_dim; e++) {
			pi[e][i] /= cnt;
			fprint_double(file, pi[e][i]);
		}
		fprintf(file, "\n");
	}
	fclose(file);


	/* Plot */
	gnuplot_ctrl *h, *g;
	h = gnuplot_init();
	// gnuplot_cmd(h, "load \'plot_trajectory.gp\'");
	gnuplot_cmd(h, "load \'plot_pi.gp\'");
	g = gnuplot_init();
	gnuplot_cmd(g, "load \'plot_pi_45.gp\'");
	// system("eog trajectory.png");
	// system("eog pi.png");
}