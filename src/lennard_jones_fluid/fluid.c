#include <assert.h>
#include <gnuplot_i.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FCC.h"
#include "util.h"

int main() {
	/*-----------------------------------------------------------
	 * definition of variables
	 *-----------------------------------------------------------
	 */

	/*  box length */
	L = pow(N / rho, 1.0 / 3.0);

	/* durantion of equilibration and of data sampling */
	double equilibration_duration = 80.0;
	int equilibration_interval = 1000;
	double sampling_duration = 80.0;
	int sampling_interval = 1000;

	/* variables for data sampling */
	int M = 0;
	double U_sample = 0.0;
	double K_sample = 0.0;
	double P_sample = 0.0;
	double dr = cutoff / S;

	/*-----------------------------------------------------------
	 * set initial conditions
	 *-----------------------------------------------------------
	 */

	/* set initial positions */
	generate_FCC(N, L, x);

	/* set velocities to zero */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			v[i][j] = 0;
		}
	}

	/* empty g(r) */
	for (int i = 0; i < S; i++) {
		g[i] = 0.0;
	}

	/* calculate  initial forces */
	calculate_forces();

	/* open output files */
	FILE *f_energy, *f_t;
	f_energy = fopen("energies.csv", "w");
	f_t = fopen("T.csv", "w");

	/* counter */
	int count = 0;

	/*-----------------------------------------------------------
	 * simulation
	 *-----------------------------------------------------------
	 */
	do {
		/* advance in trajectory */
		verlet_step();

		/* equilibration */
		if (count % equilibration_interval == 0 && count != 0 && t < equilibration_duration) {
			rescale_velocities();
		}

		/* current energies */
		double K = calculate_K();
		double U = calculate_U();
		double P = calculate_P();

		/* sampling data */
		if (t > equilibration_duration && count % sampling_interval == 0) {
			M++;
			U_sample += U;
			K_sample += K;
			P_sample += P;

			/* g(r) */
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					if (i != j) {
						double dist = r_polari(dx[i][j][0], dx[i][j][1], dx[i][j][2]);
						int pos = (int)rint(dist / dr);
						if (pos < S) {
							g[pos] += 1.0 / dr;
						}
					}
				}
			}
		}

		/* output stuff */
		/* energy */
		print_double(t, f_energy);

		print_double(K, f_energy);
		print_double(U, f_energy);
		print_double(K + U, f_energy);
		fprintf(f_energy, "\n");

		/* temperature */
		print_double(t, f_t);
		print_double(calculate_T(), f_t);
		fprintf(f_t, "\n");

		count++;
	} while (t < equilibration_duration + sampling_duration);

	/* calculate averages */
	U_sample /= M;
	K_sample /= M;
	P_sample = P_sample / (M * 6.0 * L * L * L);
	double T_sample = 2.0 / (3.0 * N) * K_sample;

	/* g(r) */
	FILE *file_g;
	double P_int = 0.0;
	file_g = fopen("g.csv", "w");
	for (int i = 0; i < S; i++) {
		double r = (double)(i * dr);
		g[i] /= 4 * M_PI * r * r * N * M * rho;
		if (i > 0) {
			P_int += 2.0 / 3.0 * M_PI * rho * rho * pow(r, 4) * lj_part(r) * g[i] * dr;
		}

		fprint_double(file_g, r);
		fprint_double(file_g, g[i]);
		fprintf(file_g, "\n");
	}
	fclose(file_g);

	/* output results */
	printf("\n----------------------------------------------\n");
	printf("RESULTS OF DATA SAMPLING\n");
	printf("----------------------------------------------\n");
	printf("Average U = %lf\nAverage K = %lf\nAverage T = %lf\nAverage P = %lf\nP with integral = %lf\n\n\n", U_sample, K_sample, T_sample, P_sample, P_int);

	/* close output files */
	fclose(f_energy);
	fclose(f_t);

	gnuplot_ctrl *h;
	h = gnuplot_init();
	gnuplot_cmd(h, "load \'plot_energies.gp\'");
	gnuplot_cmd(h, "load \'plot_T.gp\'");
	gnuplot_cmd(h, "load \'plot_g.gp\'");
}