#include "util.h"

#include <assert.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>
#include <differential_eq/verlet.h>

double r(double x, double z) {
	return sqrt(x * x + z * z);
}

double Fx(double x, double z, void *p) {
	Params *param = (Params *)p;
	double theta = param->theta;
	assert(r(x, z) != 0.0);
	return -sin(theta) / pow(r(x, z), 3) + 3.0 * x * (x * sin(theta) + z * cos(theta)) / pow(r(x, z), 5);
}

double Fz(double x, double z, void *p) {
	Params *param = (Params *)p;
	double theta = param->theta;
	assert(r(x, z) != 0.0);
	return -cos(theta) / pow(r(x, z), 3) + 3.0 * z * (x * sin(theta) + z * cos(theta)) / pow(r(x, z), 5);
}

void calculate_acc(double pos[], double a[], void *p) {
	a[0] = Fx(pos[0], pos[1], p);
	a[1] = Fz(pos[0], pos[1], p);
}

void execute_verlet(double pos[], double v[], double a[], double dt, void *p) {
	calculate_acc(pos, a, p);
	FILE *file;
	file = fopen("trajectory.csv", "w");
	while (pos[0] < 5.0 && pos[0] > -60.0) {
		verlet_2D(pos, v, a, dt, calculate_acc, p);
		// fprint_double(file, pos[0]);
		// fprint_double(file, pos[1]);
		// fprintf(file, "\n");
	}
	fclose(file);
}
