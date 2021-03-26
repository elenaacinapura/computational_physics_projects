#include <math.h>
#include <stdio.h>
#include <string.h>

#include "FCC.h"
#include "util.h"

int main() {
	// set initial positions
	generate_FCC(N, L, x);

	// set velocities to zero
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			v[i][j] = 0;
		}
	}

	// open output files
	FILE *f_energy, *f_traj;
	f_energy = fopen("Output/energies.csv", "w");
	f_traj = fopen("Output/trajectory.csv", "w");

	calculate_acc();
	do {
		
		verlet_step();


		// output stuff
		printf("%lf\n", t);
		// energy
		print_double(t, f_energy);
		print_double(calculate_K(), f_energy);
		print_double(calculate_U(), f_energy);
		print_double(calculate_K()+calculate_U(), f_energy);
		fprintf(f_energy, "\n");
		// trajectory
		print_double(t, f_traj);
		for (int i = 0; i < 3; i++) {	// print coord of particle 100
			print_double(x[100][i], f_traj);
		}
		fprintf(f_traj, "\n");
	} while (t < duration);

	//close output files
	fclose(f_energy);
}