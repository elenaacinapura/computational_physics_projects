#include <assert.h>
#include <differential_eq/runge_kutta.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

#include "util.h"

int main () {
    /* Parameters */
    double omega_0 = 0.27;
    double omega_l = 0.73;
    double t = 0.0;
    double a0 = 1.0;
    double dt = 0.001;

    Param_universe p;
    p.omega_0 = omega_0;
    p.omega_l = omega_l;

    FILE *f;
    f = fopen("universe.csv", "w");
    assert(f != NULL);
    fprintf(f, "t\ta\n");

    double t_end = 80;
    double a = a0;
    while (t < t_end) {
        a = RK4(a, t, dt, F_universe, &p);
        fprintf(f, "%lf\t",t);
        fprintf(f, "%lf", a);
        fprintf(f, "\n");

        t += dt;
    }
    fclose(f);
}