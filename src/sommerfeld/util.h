#ifndef __SOMMERFELD_UTIL_H__
#define __SOMMERFELD_UTIL_H__

#include <setjmp.h>
extern jmp_buf JumpBuffer;

typedef struct Params {
    double E, xi;
    int n;
} Params;

double integrand (double x, void *p);

double f (double E, void *p);

#endif