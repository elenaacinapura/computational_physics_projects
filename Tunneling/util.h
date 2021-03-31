#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdio.h>

#define MAXDIM (int)1e6


typedef struct Array {
    int dim;
    double v [MAXDIM];
} Array;

typedef struct Param_V_square {
    double V0, a;
} Param_V_square;


extern Array phi, X;
extern int dim;
extern double dx, E, L;

void fill_X ();

double V_square (double x, void *p);

double F (double f (double, void *), double x, void *p);

void numerov_step (int idx, double f (double, void *), void *p);

double calculate_T (int idx1, int idx2);

void fprint_vec (FILE *file, double v [], int dim);

void fprint_double (FILE *file, double d);

#endif