#ifndef __SUN_H__
#define __SUN_H__

#define START 0.0
#define EPS 1e-4
#define DIM_MAX (int)1e5

typedef struct Vec{
    double v[DIM_MAX];
    int dim;
}Vec;

typedef struct Param{
    double n;
    double theta,eta;
}Param;

double F_theta(double useless, double xi, void *p);
double F_eta(double useless, double xi, void *p);
void set_initial_condition(Vec *xi, Vec *theta, Vec *eta);
void increment_dimension(Vec *xi, Vec *theta, Vec *eta);
double integral(double xi_0, double dxi, Vec *xi, Vec *theta, Param *p);

#endif