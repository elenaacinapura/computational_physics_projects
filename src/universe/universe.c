#include <assert.h>
#include <differential_eq/runge_kutta.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

#include "util.h"


int main() {
	printf("=============================================\n");
	printf("THE AGE OF THE UNIVERSE\n");
	printf("=============================================\n");
	/* Parameters */
	double omega_0 = 0.27;
	double omega_l = 0.73;
    double inv_H0 = 3.08/70.4 * 1e19;

	Param_universe p;
	p.omega_0 = omega_0;
	p.omega_l = omega_l;
    p.sign = 0;

	FILE *f;
	f = fopen("universe.csv", "w");
	assert(f != NULL);
	fprintf(f, "t\ta\n");

	double t = 0.0;
	double dt = 0.001;
	double t_end = -3.0;
	double a0 = 1.0;
	double a = a0;
    // double t_a0;

	while (t > t_end) {
		fprintf(f, "%lf\t", t);
		fprintf(f, "%lf\n", a);
        double a_old = a;
		a = RK4(a, t, -dt, F_universe, &p);
        if (a_old >= 0.5 && a <= 0.5) {
            printf("t_1/2 = %lf\n", t);
        }
        if (a == 0.0 || !isfinite(a)) {
            printf("a(t) = 0 for reduced t = %lf\n", t);
            // t_a0 = t;
            double age = -t*inv_H0 / (3600*24*365);
            printf("The age of the universe is %e years\n", age);
            break;
        }
		t -= dt;
	}
	fclose(f);
	printf("=============================================\n");
}