#include <math.h>
#include "FCC.h"

void generate_FCC(unsigned int N, double L, double X[][3]) {
	/*
    This function fills the matrix X[N][3] with the coordinates of
    particles in a FCC lattice within a box of side L whose center is the
    center of the coordinate system.
    The box in filled uniformly if N = 4 n^3 with n integer hence for
    N = 4     32    108    256    500    864   1372   2048   2916
    4000   5324 6912   8788  10976  13500  16384 ...
  */

	int i, j, k, m, n, p, c;

	/* position within the primary cell */
	const double rFCC[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
	double b, rCell[3];

	for (c = 1;; c++)
		if (4 * c * c * c >= N)
			break;

	b = L / (double)c; /* side length of the primary cell */
	p = 0;			   /* particles placed so far */
	for (i = 0; i < c; i++) {
		rCell[0] = i;
		for (j = 0; j < c; j++) {
			rCell[1] = j;
			for (k = 0; k < c; k++) {
				rCell[2] = k;
				for (m = 0; m < 4; m++) {/* 4 particles in cell */
					if (p < N) {
						/* add the com to each bead, and project to the real cell */
						for (n = 0; n < 3; n++) {
							X[p][n] = b * (rCell[n] + rFCC[m][n]);
							X[p][n] -= L * rint(X[p][n] / L);
						}
						++p;
					}
                }
			}
		}
	}
}