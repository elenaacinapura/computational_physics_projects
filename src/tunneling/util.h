#ifndef __UTIL_H__
#define __UTIL_H__

#include <complex.h>
#include <stdbool.h>
#include <stdio.h>

/* STRUCTURES */
typedef struct Param_F {
    double xi; /* xi = 2*m*a*a*|V_0| / h*h */
} Param_F;

/* VARIABLES */
extern double E; /* in reduced units */

/* FUNCTIONS */
double F_square (double x, void *p);

double F_gauss (double x, void *p);

double F_asymm_L (double x, void *p);

double F_asymm_R (double x, void *p);

double F_cosh (double x, void *p);

void solve_numerov (double x[], complex double phi[], int dim, double dx, double F (double, void *p), void *p, bool printoutput, FILE *outfile);

complex double T_coeffcient(double x[], complex double phi[], int dim, double k);


#endif