#include <assert.h>
#include <differential_eq/verlet.h>
#include <math.h>
#include <numerical_methods/integrate_notypesafe.h>
#include <numerical_methods/zero_newton.h>
#include <print_routines.h>

/*=========================== PARAMETERS ======================*/
double E = 1.0;
double x_0 = 20.0;
double dt = 1e-3;
double dtheta = 0.5;
double b_start = 0.0;
double db = 0.01;
double b_end = 5.0;
double b;
/*=========================== STRUCT ======================*/
typedef struct Param_integral {
	double b, E;
} Param_integral;
/*=========================== FUNCTIONS ======================*/
void calculate_acc(double pos[2], double a[2], void *p) {
	double x = pos[0];
	double y = pos[1];
	a[0] = 24.0 * (2.0 * pow(x * x + y * y, -7) - pow(x * x + y * y, -4)) * x;
	a[1] = 24.0 * (2.0 * pow(x * x + y * y, -7) - pow(x * x + y * y, -4)) * y;
}
void execute_verlet(double pos[], double v[], double a[]) {
	struct empty {
	} useless_struct;
	while (fabs(pos[0]) < 2 * x_0) {
		verlet_2D(pos, v, a, dt, calculate_acc, &useless_struct);
	}
}
double integrand(double x, void *p) {
	return b / (x * x * sqrt(1.0 - b * b / (x * x) - 4.0 / E * (pow(x, -12) - pow(x, -6))));
}
double f(double x, void *p) {
	return pow(x, 12) - b * b * pow(x, 10) + 4.0 / E * pow(x, 6) - 4.0 / E;
}
double find_theta_analytically(double b) {
	struct Empty {
	} useless;
	double rmin = zero_newton(f, 10.0, &useless);
	double theta = M_PI - 2.0 * integrate(integrand, rmin, 100.0, 1000, &useless);
	return theta / M_PI * 180.0;
}
/*=========================== MAIN ======================*/
int main() {
	/* Welcome */
	printf("===================================================\n");
	printf("CLASSICAL SCATTERING\n");
	printf("===================================================\n");
	printf("Parameters:\n");
	printf("	E = %.1lf\n", E);
	printf("	db = %lf\n", db);
	printf("	b = [%.1lf, %.1lf]\n", b_start, b_end);
	printf("	dt = %lf\n", dt);
	printf("	dtheta = %.2lf\n", dtheta);
	printf("===================================================\n");

	double pos[2], v[2], a[2];
	int n_bins = (int)(180 / dtheta);
	int dim = (int)((b_end - b_start) / db);
	double b_vec[dim];
	double theta_vec[dim];
	double dtheta_db[dim];
	double cross_section[n_bins];

	for (int i = 0; i < n_bins; i++) {
		cross_section[i] = 0.0;
	}

	FILE *f_theta, *f_cs;
	f_theta = fopen("theta.csv", "w");
	fprintf(f_theta, "b\ttheta\ttheta_theo\n");

	b = b_start;
	for (int cnt = 0; cnt < dim; cnt++) {
		fflush(stdout);
		printf("\rcurrent b = %.3lf", b);
		b_vec[cnt] = b;
		/* set initial conditions */
		pos[0] = -x_0;
		pos[1] = b;
		v[0] = sqrt(2.0 * E);
		v[1] = 0.0;
		/* execute verlet */
		execute_verlet(pos, v, a);
		/* extract the angle */
		double theta = atan2(v[1], v[0])   / M_PI * 180.0; /* in degrees */
		theta_vec[cnt] = theta;
		double theta_theo = find_theta_analytically(b);
		/* print theta */
		fprint_double(f_theta, b);
		fprint_double(f_theta, theta);
		fprint_double_newline(f_theta, theta_theo);

		b += db;
	}
	fclose(f_theta);
	/* ===================== REDUCED CROSS SECTION ======================= */
	/* calculate dtheta/db */
	dtheta_db[0] = (theta_vec[1] - theta_vec[0]) / db;
	for (int i = 1; i < dim; i++) {
		dtheta_db[i] = fabs((theta_vec[i] - theta_vec[i - 1]) / db);
		// fprint_double_newline(stdout, dtheta_db[i]);
	}
	/* fill the bin */
	for (int i = 0; i < dim; i++) {
		double theta = theta_vec[i];
		if (theta < 0) {
			theta = -theta;
		}
		int idx = (int)((theta) / dtheta);
		double weight = 2 * M_PI * b_vec[i] / dtheta_db[i];
		if (idx > 0 && idx < n_bins) {
			cross_section[idx] += weight;
		}
	}
	/* print cross section */
	f_cs = fopen("cross_section.csv", "w");
	fprintf(f_cs, "theta\tcnt\n");
	for (int i = 0; i < n_bins; i++) {
		fprint_double(f_cs, i * dtheta);
		fprint_double_newline(f_cs, cross_section[i]);
	}
	fclose(f_cs);
}