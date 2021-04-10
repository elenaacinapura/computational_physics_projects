#include "util.h"

#include <gnuplot_i.h>
#include <integrate_notypesafe.h>
#include <print_routines.h>
#include <zeros_newton.h>

#include <assert.h>
#include <math.h>
#include <setjmp.h>
#include <stdio.h>

jmp_buf JumpBuffer;

double integrand (double x, void *p) {
    Params *params = (Params *) p;
    double E = params->E;
    assert(E + pow(cosh(x), -4) >= 0.0);
    return sqrt(E + pow(cosh(x), -4));
}

/* Functions of which I will find the zeros */
double f (double E, void *p) {
    Params *params = (Params *) p;
    params->E = E;
    double xi = params-> xi;
    int n = params->n;

    if (E >= 0.0) {
        longjmp(JumpBuffer, 1);
    }

    double high = acosh(pow(-E, -0.25));
    double low = -high;
    return integrate(integrand, low, high, 100, p) - M_PI* (n + 0.5) * sqrt(xi);
}

