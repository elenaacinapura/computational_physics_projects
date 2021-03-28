#include <math.h>
#include <stdio.h>
#include <string.h>

#include "FCC.h"
#include "util.h"

int main() {

	// calculate box length
	L = pow(N/rho,1.0/3.0);
	printf("L/2 = %lf\n",L/2);

	// set initial positions
	generate_FCC(N, L, x);

	// set velocities to zero
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			v[i][j] = 0;
		}
	}

	// open output files
	FILE *f_energy, *f_t, *f_traj;
	f_energy = fopen("Output/energies.csv", "w");
	f_t = fopen("Output/T.csv", "w");
	f_traj = fopen("Output/traj.csv", "w");

	calculate_forces();
	int count = 0;
	do {
		
		if( count!=0 && count%1000 == 0 ){
			rescale_velocities();
		}

		verlet_step();

		// energy
		print_double(t, f_energy);
		double K = calculate_K();
		double U = calculate_U();
		print_double(K, f_energy);
		print_double(U, f_energy);
		print_double(K+U, f_energy);
		fprintf(f_energy, "\n");
		
		// temperature
		print_double(t, f_t);
		print_double(calculate_T(),f_t);
		fprintf(f_t, "\n");

		// treiettoria
		print_double(t, f_traj);
		for (int i = 0; i < 3; i++) {	// print coord of particle 100
			print_double(x[100][i], f_traj);
		}
		fprintf(f_traj, "\n");

		count++;
	} while (t < duration);

	//close output files
	fclose(f_energy);
	fclose(f_t);
	fclose(f_traj);
	
}