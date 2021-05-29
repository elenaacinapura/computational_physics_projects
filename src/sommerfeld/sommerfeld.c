#include "util.h"

#include <gnuplot_i.h>
#include <numerical_methods/integrate_notypesafe.h>
#include <print_routines.h>
#include <numerical_methods/zero_newton.h>

#include <assert.h>
#include <math.h>
#include <setjmp.h>
#include <stdbool.h>
#include <stdio.h>

int main () {
    Params p;
    int num_xi = 3;
    double xi[] = {0.025, 0.005, 0.0025};

    /* Welcome */
    printf("========================================================\n");
    printf("BOHR-SOMMERFELD METHOD FOR FINDING BOUND STATES\n");
    printf("========================================================\n");
    printf("Values of xi = hbar^2 / (m a^2 V0) considered:\n");
    for (int i = 0; i < num_xi; i++) {
        printf("%.4lf\t", xi[i]);
    }
    printf("\n");
    printf("========================================================\n");
    printf("Calculating...\n");
    bool stop = 0;
    for (int i = 0; i < num_xi; i++) {
        p.xi = xi[i];
        printf("---------------------------------------------\nxi = %.04lf\n", xi[i]);
        for (int n = 0; n < 15; n++) {
            if (stop) {
                stop = 0;
                break;
            }
            p.n = n;
            if (setjmp(JumpBuffer)) {
                printf("End of bound states reached for n = %d\n", n);
                stop = 1;
                continue;
            }
            double E = zero_newton(f, -0.0001, &p);
            printf("n = %d\t E = %.4lf\n", n, E);
        }
    }
    printf("========================================================\n");
}