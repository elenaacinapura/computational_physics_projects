#ifndef __VARIATIONAL_UTIL_H__
#define __VARIATIONAL_UTIL_H__

#include <stdbool.h>
#include <stdio.h>

/* STRUCTURES */
typedef struct Param_E{
    double xi;
    double alpha;
}Param_E;

/* VARIABLES */
extern double alpha_max, dalpha;

/* FUNCTIONS */
double integrand(double x,void *p);
double integrate_trap (double f (double, void *), double low, double high, int density, void *param);
void fill_alpha(double alpha[],int dim);


#endif