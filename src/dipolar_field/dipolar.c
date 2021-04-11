#include <gnuplot_i.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>
#include <verlet.h>

#include "util.h"

int main() {
	double dt = 1e-3;

	double E = 1.0;
	double theta = 0.0;
	double x0 = -10.0;
	Params par = {E, theta};

	double pos[2], v[2], a[2];

	double cutoff = 25.0;
	double dz = 0.1;
	int dim = 2 * (int)rint(cutoff / dz);
	double pi[dim];

	for (int i = 0; i < dim; i++) {
		pi[i] = 0.0;
	}

	double z0 = 25.0;
	double z = -z0;
	double dz_start = 0.005;
	int cnt = 0;
	while (z <= z0) {
		pos[0] = x0;
		pos[1] = z;
		v[0] = sqrt(2.0 * E);
		v[1] = 0.0;
		execute_verlet(pos, v, a, dt, &par);

		int idx = (int)rint(pos[1] / dz) + dim / 2;
		if (idx >= 0 && idx < dim) {
			pi[idx] += 1.0;
		}
		z += dz_start;
		cnt++;
	}

	FILE *file;
	file = fopen("pi.csv", "w");
	for (int i = 0; i < dim; i++) {
		pi[i] /= cnt;
		fprintf(file, "%lf\t", -cutoff + dz * i);
        fprintf(file, "%lf\n", pi[i]);
	}
	fclose(file);

	/* Plot */
	gnuplot_ctrl *h;
	h = gnuplot_init();
	// gnuplot_cmd(h, "load \'plot_trajectory.gp\'");
	gnuplot_cmd(h, "load \'plot_pi.gp\'");
	// system("eog trajectory.png");
	system("eog pi.png");
}