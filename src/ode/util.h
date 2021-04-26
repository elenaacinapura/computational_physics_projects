#ifndef __ODE_UTIL_H__
#define __ODE_UTIL_H__

typedef struct Param_universe {
    double omega_0, omega_l;
} Param_universe;

double F_universe (double a, double t, void *param);

#endif