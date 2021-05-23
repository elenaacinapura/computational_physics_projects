#ifndef __BOUND_STATES_NUMEROV_ROUTINES_H__
#define __BOUND_STATES_NUMEROV_ROUTINES_H__

void execute_numerov(double x[], double phi [], double dx, int dim, double F (double, void*), void *p);

double calculate_delta (double x0, double phiF[], double phiB[], double dx, int dimF, int dimB, double F (double, void *), void *p);

#endif