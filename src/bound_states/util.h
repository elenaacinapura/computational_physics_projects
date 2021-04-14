#ifndef __BOUND_STATES_UTIL_H__
#define __BOUND_STATES_UTIL_H__

typedef struct Params {
    double a, E;
} Params;

double F_cosh (double x, void *param);

void execute_numerov(double x[], double phi [], double dx, int dim, double F (double, void*), void *p);

double calculate_delta (double x0, double phiF[], double phiB[], double dx, int dimF, int dimB, double F (double, void *), void *p);

#endif