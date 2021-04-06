#ifndef __UTIL_H__
#define __UTIL_H__

#include <complex.h>
#include <stdbool.h>
#include <stdio.h>

/* STRUCTURES */
typedef struct Param_F_square {
    double sqrt_xi;
} Param_F_square;

/* VARIABLES */
extern double E; /* in reduced units */

/* FUNCTIONS */
double F_square (double x, void *p);

void solve_numerov (double x[], complex double phi[], int dim, double dx, double F (double, void *p), void *p, bool printoutput, FILE *outfile);

complex double T_coeffcient(double x[], complex double phi[], int dim, double k);


#endif