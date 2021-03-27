#include <math.h>
#include <stdio.h>
#include <string.h>

#include "FCC.h"
#include "util.h"

int main() {

	// calculate box length
	L = pow(N/rho,1.0/3.0);

	// set initial positions
	generate_FCC(N, L, x);

	// set velocities to zero
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			v[i][j] = 0;
		}
	}



	// open output files
	FILE *f_energy, *f_t;
	f_energy = fopen("Output/energies.csv", "w");
	f_t = fopen("Output/T.csv", "w");

	calculate_forces();
	do {
		
		verlet_step();

		// energy
		print_double(t, f_energy);
		print_double(calculate_K(), f_energy);
		print_double(calculate_U(), f_energy);
		print_double(calculate_K()+calculate_U(), f_energy);
		fprintf(f_energy, "\n");
		
		// temperature
		print_double(t, f_t);
		print_double(calculate_T(),f_t);
		fprintf(f_t, "\n");

	} while (t < duration);

	//close output files
	fclose(f_energy);
	fclose(f_t);
	
	print_mat(x);
}