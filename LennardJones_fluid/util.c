#include "util.h"

#include <math.h>
#include <stdio.h>

double x[N][3];
double v[N][3];
double a[N][3];

double rho = 0.7;
double cutoff = 3.0;
double T = 0.7;
double dt = 0.001;
double t = 0.0;
double duration = 1.0;


void print_mat(double m[N][3]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%lf\t", m[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_double(double d, FILE *f) {
	fprintf(f, "%lf\t", d);
}

double r_polari(double x, double y, double z) {
	return sqrt(x * x + y * y + z * z);
}

double lj_part(double r) {
	return 24 * (2 * pow(r, -14) - pow(r, -8));
}

double lj_u(double r) {
	return 4 * (pow(r, -12) - pow(r, -6));
}

void calculate_acc() {
	double L = pow(N/rho,1/3);
	// reset accelerations
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			a[i][j] = 0.0;
		}
	}
	// update accelerations
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) {
			// distances
			double dx = x[i][0] - x[j][0];
			double dy = x[i][1] - x[j][1];
			double dz = x[i][2] - x[j][2];

			dx -= L * rint(dx / L);
			dy -= L * rint(dy / L);
			dz -= L * rint(dz / L);

			double r = r_polari(dx, dy, dz);

			if (r < cutoff) {
				// common part in force
				double acc_part = lj_part(r);
				// acc on particle i due to j
				a[i][0] += acc_part * dx;
				a[i][1] += acc_part * dy;
				a[i][2] += acc_part * dz;
				// acc on particle j due to i
				a[j][0] += -a[i][0];
				a[j][1] += -a[i][1];
				a[j][2] += -a[i][2];
			}
		}
	}
}

void verlet_step() {
	double L = pow(N/rho,1/3);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			v[i][j] += dt / 2.0 * a[i][j];
			x[i][j] += dt * v[i][j];
			x[i][j] -= L * rint(x[i][j] / L);  //pacman
		}
		calculate_acc();
		for (int j = 0; j < 3; j++) {
			v[i][j] += dt / 2.0 * a[i][j];
		}
	}
	t += dt;
}

double calculate_U() {
	double L = pow(N/rho,1/3);
	double U = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				continue;
			}
			// distances
			double dx = x[i][0] - x[j][0];
			double dy = x[i][1] - x[j][1];
			double dz = x[i][2] - x[j][2];
			dx -= L * rint(dx / L);
			dy -= L * rint(dy / L);
			dz -= L * rint(dz / L);
			double r = r_polari(dx, dy, dz);

			U += lj_u(r);
		}
	}
	return U;
}

double calculate_K() {
	double K_curr = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			K_curr += 0.5 * pow(v[i][j], 2);
		}
	}
	return K_curr;
}

double calculate_T() {
	double K_curr = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			K_curr += 0.5 * pow(v[i][j], 2);
		}
	}
	double T_curr = 2.0 / (3.0 * N) * K_curr;
	return T_curr;
}

void rescale_velocities() {
	double Ko = 3.0 / 2.0 * N * T;
	double T_curr = calculate_T();
	double alpha = sqrt(T / T_curr);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 3; j++) {
			v[i][j] *= alpha;
		}
	}
}