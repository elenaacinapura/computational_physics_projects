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

	calculate_acc();

	do {
		verlet_step();

        print_double(t, "Output/T.csv");
		print_double(calculate_T(), "Output/T.csv");
        FILE *outputfile;
        outputfile = fopen("Output/T.csv", "a");
        fprintf(outputfile, "\n");
        fclose(outputfile);

        char file_energies [] = "Output/energies.csv";
        print_double(t, file_energies);
        print_double(calculate_U() + calculate_K(), file_energies);
        outputfile = fopen(file_energies, "a");
        fprintf(outputfile, "\n");
        fclose(outputfile);

	} while (t < duration);
}