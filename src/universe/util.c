#include "util.h"

#include <assert.h>
#include <differential_eq/runge_kutta.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

double F_universe (double a, double t, void *param) {
    Param_universe *p = (Param_universe *) param;
    double omega_0 = p->omega_0;
    double omega_l = p->omega_l;
    int sign = p->sign;
    return pow(-1, sign) * sqrt(omega_0/a + omega_l*a*a);
}
