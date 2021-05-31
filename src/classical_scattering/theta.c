#include <assert.h>
#include <differential_eq/verlet.h>
#include <math.h>
#include <print_routines.h>
#include <numerical_methods/integrate_notypesafe.h>
#include <numerical_methods/zero_newton.h>

/*=========================== PARAMETERS ======================*/
double E = 1.0;
double x_0 = 20.0;
double dt = 1e-3;
double dtheta = 1.0;
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
	return b/(x*x * sqrt(1.0 - b*b/(x*x) - 4.0/E*(pow(x, -12) - pow(x, -6))));
}
double f(double x, void *p) {
	return pow(x, 12) - b*b*pow(x, 10) + 4.0/E*pow(x, 6) - 4.0/E;
}
double find_theta_analytically(double b, double E) {
    struct Empty {} useless;
	double rmin = zero_newton(f, 5.0, &useless);
	double theta = M_PI - 2.0 * integrate(integrand, rmin, 100.0, 1000, &useless);
	return theta / M_PI * 180.0;
}
/*=========================== MAIN ======================*/
int main() {
	double pos[2], v[2], a[2];
	int n_bins = (int)(180 / dtheta);
	double cross_section[n_bins];

	for (int i = 0; i < n_bins; i++) {
		cross_section[i] = 0.0;
	}

	FILE *f_theta, *f_cs;
	f_theta = fopen("theta.csv", "w");
	fprintf(f_theta, "b\ttheta\ttheta_theo\n");

	b = b_start;
	while (b <= b_end) {
		fflush(stdout);
		printf("\rcurrent b = %.2lf", b);
		/* set initial conditions */
		pos[0] = -x_0;
		pos[1] = b;
		v[0] = sqrt(E / 2.0);
		v[1] = 0.0;
		/* execute verlet */
		execute_verlet(pos, v, a);
		/* extract the angle */
		double theta = atan2(v[1], v[0]) / M_PI * 180.0;
        double theta_theo = find_theta_analytically(b, E);
		/* print theta */
		fprint_double(f_theta, b);
		fprint_double(f_theta, theta);
        fprint_double_newline(f_theta, theta_theo);
		/* fill the bin */
		int idx = (int)theta / dtheta;
		double weight = 2 * M_PI * b * db / dtheta;
		if (idx > 0 && idx < n_bins) {
			cross_section[idx] += weight;
		}
		b += db;
	}
	fclose(f_theta);
	/* print cross section */
	f_cs = fopen("cross_section.csv", "w");
	fprintf(f_cs, "theta\tcnt\n");
	for (int i = 0; i < n_bins; i++) {
		fprint_double(f_cs, i * dtheta);
		fprint_double_newline(f_cs, cross_section[i]);
	}
	fclose(f_cs);
}