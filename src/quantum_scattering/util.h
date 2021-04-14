#ifndef __QUANTUM_SCATTERING_H__
#define __QUANTUM_SCATTERING_H__

typedef struct Params {
    double xi, E;
    int l;
} Params;

double F_hard(double r, void *param);

void execute_numerov (double x[], double u[], int dim, double dx, double F (double, void *), void *p);

double find_delta (double r1, double r2, double u1, double u2, int l, double k);

double sin_2_delta_th(double k, int l);


#endif