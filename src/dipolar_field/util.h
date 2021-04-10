#ifndef __DIPOLAR_UTIL_H__
#define __DIPOLAR_UTIL_H__

typedef struct Params {
    double E, theta;
} Params;

double Fx (double x, double z, void *p);

double Fz (double x, double z, void *p);

void calculate_acc (double pos[], double a[], void *p);

void execute_verlet (double pos[], double v[], double a[], double dt, void *p);

#endif