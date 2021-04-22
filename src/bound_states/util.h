#ifndef __BOUND_STATES_UTIL_H__
#define __BOUND_STATES_UTIL_H__

typedef struct Params_cosh {
    double a, E;    /* a = 2 m s^2 V0 / hbar^2 */
} Params_cosh;

typedef struct Params_lj {
    double a, E;
    int l;
} Params_lj;

typedef struct Params_delta {
    double L, dx, x0, a, A, B;
    int l;
} Params_delta;

double V_cosh(double x, void *param);

double F_cosh (double x, void *param);

double F_lj(double r, void *param);

void execute_numerov(double x[], double phi [], double dx, int dim, double F (double, void*), void *p);

double calculate_delta (double x0, double phiF[], double phiB[], double dx, int dimF, int dimB, double F (double, void *), void *p);

double Delta_E_cosh(double E, void *param);

double Delta_E_lj(double E, void *param);

#endif