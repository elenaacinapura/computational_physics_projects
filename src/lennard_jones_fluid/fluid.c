#include <math.h>
#include <stdio.h>
#include <string.h>

#include "FCC.h"
#include "util.h"

int main() {

	/*-----------------------------------------------------------
	 * definition of variables
	 *-----------------------------------------------------------
	 */

	/*  box length */
	L = pow(N/rho,1.0/3.0);
	printf("L/2 = %lf\n",L/2);

	/* durantion of equilibration and of data sampling */
	double equilibration_duration = 40.0;
	int equilibratio_interval = 1000;
	double sampling_duration = 40.0;
	int sampling_interval = 1000;

	/* variables for data sampling */
	int M = 0; 
	double U_sample = 0.0;
	double K_sample = 0.0;

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

	/* calculate  initial forces */
	calculate_forces();

	/* open output files */
	FILE *f_energy, *f_t, *f_traj;
	f_energy = fopen("Output/energies.csv", "w");
	f_t = fopen("Output/T.csv", "w");
	f_traj = fopen("Output/traj.csv", "w");

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
		if(count%equilibratio_interval == 0 && count != 0 && t < equilibration_duration){
			rescale_velocities();
		}

		/* current energies */
		double K = calculate_K();
		double U = calculate_U(); 

		/* sampling data */
		if (t > equilibration_duration && count%sampling_interval == 0){
			M++;
			U_sample += U;
			K_sample += K;
		}

		/* output stuff */
		/* energy */
		print_double(t, f_energy);
		
		print_double(K, f_energy);
		print_double(U, f_energy);
		print_double(K+U, f_energy);
		fprintf(f_energy, "\n");
		
		/* temperature */
		print_double(t, f_t);
		print_double(calculate_T(),f_t);
		fprintf(f_t, "\n");

		/* trajectory of particle 100 */
		print_double(t, f_traj);
		for (int i = 0; i < 3; i++) {
			print_double(x[100][i], f_traj);
		}
		fprintf(f_traj, "\n");

		count++;
	} while (t < equilibration_duration + sampling_duration);

	/* calculate averages */ 
	U_sample /= M;
	K_sample /= M;
	double T_sample = 2.0/(3.0*N)*K_sample;

	/* output results */
	printf("\nRESULTS OF DATA SAMPLING\n");
	printf("Average U = %lf\nAverage K = %lf\nAverage T = %lf", U_sample, K_sample, T_sample);

	/* close output files */
	fclose(f_energy);
	fclose(f_t);
	fclose(f_traj);
	
}