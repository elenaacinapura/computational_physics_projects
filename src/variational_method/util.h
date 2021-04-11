#ifndef __VARIATIONAL_UTIL_H__
#define __VARIATIONAL_UTIL_H__

#include <stdio.h>

typedef struct Params {
    double xi, alpha;
} Params;

double integrand (double x, void *p);

double energy (double a, void *p);

#endif